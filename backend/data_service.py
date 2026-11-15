"""SPHEREx cutout data service.

Queries IRSA SIA2 for SPHEREx spectral-image MEFs overlapping a sky
position, downloads small cutouts via IRSA's cutout service, caches them
locally, renders PNG previews, and returns chronological frame metadata.

Public entry point:
    get_spherex_cutouts(ra, dec, fov_size, ...)

Verified against official IRSA docs (see phase-0/research/).
"""

from __future__ import annotations

import hashlib
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import astropy.units as u
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa import Irsa

# Allow `python backend/data_service.py` smoke tests and imports from repo root.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from processing.render import render_fits_png  # noqa: E402

CACHE_DIR = Path(__file__).resolve().parent / "cache"
MANIFEST_PATH = CACHE_DIR / "session_manifest.json"
OFFLINE_COORD_TOL_DEG = 0.05   # ~3 arcmin: treat as same target
COLLECTIONS = ["spherex_dr1", "spherex_dr1_deep", "spherex_qr3", "spherex_qr3_deep"]
DEFAULT_FOV_DEG = 0.05          # ~30 arcsec cutout (verified in Phase 0)
DEFAULT_MAX_FRAMES = 40         # cap downloads per query for responsiveness
HTTP_HEADERS = {"User-Agent": "planet-x-spherex/0.1 (NASA Space Apps 2026)"}
HTTP_TIMEOUT = 120


def _cache_key(url: str) -> str:
    return hashlib.sha256(url.encode()).hexdigest()[:16]


def _mjd_to_iso(mjd) -> str | None:
    """Convert IRSA t_min/t_max (MJD) to ISO-8601 UTC."""
    try:
        from astropy.time import Time
        return Time(float(mjd), format="mjd").iso
    except Exception:  # noqa: BLE001
        return None


def _fetch(url: str, dest: Path) -> None:
    """Download URL to dest with retries. Skips if dest already exists."""
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_exc: Exception | None = None
    for attempt in range(4):
        try:
            with requests.get(url, headers=HTTP_HEADERS,
                              timeout=HTTP_TIMEOUT, stream=True) as r:
                r.raise_for_status()
                with tmp.open("wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            tmp.replace(dest)
            return
        except Exception as exc:  # noqa: BLE001 - retry transient errors
            last_exc = exc
            time.sleep(2 * (attempt + 1))
    if tmp.exists():
        tmp.unlink(missing_ok=True)
    raise RuntimeError(f"download failed for {url}: {last_exc}")


def _process_row(row, ra: float, dec: float, fov_size: float,
                 index: int) -> dict | None:
    """Download + render one MEF cutout. Returns frame metadata dict."""
    base_url = row["access_url"]
    cutout_url = f"{base_url}?center={ra},{dec}&size={fov_size}"
    key = _cache_key(cutout_url)
    fits_path = CACHE_DIR / f"{key}.fits"
    png_path = CACHE_DIR / f"{key}.png"

    _fetch(cutout_url, fits_path)
    if not png_path.exists():
        render_fits_png(fits_path, png_path)

    obs_id = str(row.get("obs_id", ""))
    bandpass = str(row.get("energy_bandpassname", ""))
    t_min = float(row["t_min"]) if row.get("t_min") is not None else None
    t_max = float(row["t_max"]) if row.get("t_max") is not None else None

    return {
        "index": index,
        "obs_id": obs_id,
        "bandpass": bandpass,
        "detector": bandpass.replace("SPHEREx-", ""),
        "t_min": t_min,
        "t_max": t_max,
        "t_min_iso": _mjd_to_iso(t_min),
        "t_max_iso": _mjd_to_iso(t_max),
        "fits_path": str(fits_path),
        "png_path": str(png_path),
        "access_url": base_url,
        "cutout_url": cutout_url,
    }


def get_spherex_cutouts(ra: float, dec: float,
                        fov_size: float = DEFAULT_FOV_DEG,
                        collections: list[str] | None = None,
                        max_frames: int | None = DEFAULT_MAX_FRAMES,
                        workers: int = 6,
                        search_radius_arcmin: float = 2.0) -> dict:
    """Query IRSA for SPHEREx MEFs at (ra, dec) and prepare cutout frames.

    Parameters
    ----------
    ra, dec : float           ICRS degrees.
    fov_size : float          Cutout size in degrees (IRSA cutout service).
    collections : list[str]   SIA2 collections; defaults to DR1+QR3 spectral
                              images (wide + deep).
    max_frames : int | None   Cap on downloaded frames; None = all matches.
    workers : int             Parallel download workers.
    search_radius_arcmin : float  SIA2 cone-search radius.

    Returns
    -------
    dict with keys:
        ra, dec, fov_size, query_time, n_results (SIA rows),
        n_epochs (unique t_min), frames (chronological list of dicts),
        error (str or None).
    """
    out = {
        "ra": ra, "dec": dec, "fov_size": fov_size,
        "query_time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_results": 0, "n_epochs": 0, "frames": [], "error": None,
        "offline": False,
    }

    results = None
    last_exc = None
    for attempt in range(3):
        try:
            coord = SkyCoord(ra, dec, unit="deg", frame="icrs")
            radius = search_radius_arcmin * u.arcmin
            cols = collections or COLLECTIONS
            results = Irsa.query_sia(pos=(coord, radius), collection=cols)
            break
        except Exception as exc:  # noqa: BLE001 - retry transient NASA errors
            last_exc = exc
            print(f"[data_service] SIA2 attempt {attempt + 1} failed: {exc}",
                  file=sys.stderr)
            time.sleep(3 * (attempt + 1))

    if results is None or len(results) == 0:
        # Offline/cached fallback: serve the last good session if the
        # requested position matches it (e.g., the WISE 0855 demo target).
        cached = _load_offline_session(ra, dec, fov_size)
        if cached is not None:
            print("[data_service] serving OFFLINE cached session "
                  f"({cached['n_epochs']} epochs)", file=sys.stderr)
            cached["offline"] = True
            cached["error"] = (f"IRSA unreachable ({last_exc}); "
                               "serving cached data")
            return cached
        out["error"] = (f"SIA2 query failed after retries: {last_exc}"
                        if last_exc else "SIA2 query returned no results")
        return out

    out["n_results"] = len(results)
    tmins = results["t_min"] if "t_min" in results.colnames else []
    out["n_epochs"] = len(set(float(t) for t in tmins))

    # Sort by time, then deduplicate to one row per (t_min, detector/band)
    # is unnecessary — keep all rows; the UI can pick epochs.
    order = np.argsort(results["t_min"]) if "t_min" in results.colnames else range(len(results))
    rows = [results[i] for i in order]
    if max_frames is not None:
        rows = rows[:max_frames]

    frames: list[dict] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_process_row, row, ra, dec, fov_size, i): i
            for i, row in enumerate(rows)
        }
        for fut in as_completed(futures):
            try:
                meta = fut.result()
                if meta:
                    frames.append(meta)
            except Exception as exc:  # noqa: BLE001 - skip bad frames
                print(f"[data_service] frame failed: {exc}", file=sys.stderr)

    frames.sort(key=lambda f: (f["t_min"] is None, f["t_min"] or 0.0))
    for i, f in enumerate(frames):
        f["index"] = i
    out["frames"] = frames
    _save_manifest(out)
    return out


def _save_manifest(result: dict) -> None:
    """Persist the latest successful session for offline fallback."""
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        slim = {
            "ra": result["ra"], "dec": result["dec"],
            "fov_size": result["fov_size"],
            "n_results": result["n_results"],
            "n_epochs": result["n_epochs"],
            "query_time": result["query_time"],
            "frames": result["frames"],
        }
        MANIFEST_PATH.write_text(__import__("json").dumps(slim))
    except Exception as exc:  # noqa: BLE001
        print(f"[data_service] manifest save failed: {exc}", file=sys.stderr)


def _load_offline_session(ra: float, dec: float,
                          fov_size: float) -> dict | None:
    """Return cached session if it matches the requested position."""
    try:
        if not MANIFEST_PATH.exists():
            return None
        import json
        m = json.loads(MANIFEST_PATH.read_text())
        if (abs(m["ra"] - ra) > OFFLINE_COORD_TOL_DEG or
                abs(m["dec"] - dec) > OFFLINE_COORD_TOL_DEG):
            return None
        # Verify cached artifacts still exist
        frames = [f for f in m["frames"]
                  if Path(f["png_path"]).exists()]
        if not frames:
            return None
        for i, f in enumerate(frames):
            f["index"] = i
        return {
            "ra": ra, "dec": dec, "fov_size": fov_size,
            "query_time": m.get("query_time"),
            "n_results": m.get("n_results", len(frames)),
            "n_epochs": m.get("n_epochs", 0),
            "frames": frames, "error": None, "offline": True,
        }
    except Exception as exc:  # noqa: BLE001
        print(f"[data_service] offline load failed: {exc}", file=sys.stderr)
        return None


if __name__ == "__main__":
    import json
    res = get_spherex_cutouts(133.795129, -7.245147, max_frames=5)
    res["frames"] = [{k: v for k, v in f.items()
                      if k not in ("fits_path", "png_path", "access_url",
                                   "cutout_url")} for f in res["frames"]]
    print(json.dumps(res, indent=2, default=str))
