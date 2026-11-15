"""Catalog cross-matching service (Phase 4).

Queries Gaia DR3 (via astroquery.gaia) and AllWISE (via IRSA TAP region
query) around a target position, propagates proper motions to a reference
epoch, and classifies the target:

  - "known_object"        : a catalog source lies within tolerance
  - "discovery_candidate" : no catalog match (UI must present this as a
                            *candidate*, never a confirmed discovery)

All catalog calls are best-effort: failures are captured in `warnings`
and never break the endpoint.
"""

from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutTimeout
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

MATCH_TOLERANCE_ARCSEC = 3.0     # after proper-motion propagation
DEFAULT_SEARCH_RADIUS_ARCSEC = 10.0
GAIA_EPOCH_JYEAR = 2016.0        # Gaia DR3 reference epoch
ALLWISE_EPOCH_JYEAR = 2010.0     # approximate AllWISE mean epoch


def _propagate(ra_deg, dec_deg, pmra_masyr, pmdec_masyr,
               from_epoch, to_epoch):
    """Propagate ICRS position by proper motion (linear, small-angle)."""
    dt = to_epoch - from_epoch
    if pmra_masyr is None or pmdec_masyr is None or not np.isfinite(dt):
        return ra_deg, dec_deg
    try:
        pmra = float(pmra_masyr)
        pmdec = float(pmdec_masyr)
    except (TypeError, ValueError):
        return ra_deg, dec_deg
    if not (np.isfinite(pmra) and np.isfinite(pmdec)):
        return ra_deg, dec_deg
    dec = np.deg2rad(dec_deg)
    dra = pmra * dt / 3.6e6 / np.cos(dec)   # pmra is mu_alpha* (incl. cos dec)
    ddec = pmdec * dt / 3.6e6
    return ra_deg + dra, dec_deg + ddec


CATALOG_TIMEOUT_SEC = 150


def _with_timeout(fn, label: str, warnings: list[str]):
    """Run fn() in a thread; return result or None on timeout/error.

    Uses shutdown(wait=False) so a hung upstream query cannot block the
    endpoint forever (the orphaned thread finishes whenever it can).
    """
    pool = ThreadPoolExecutor(max_workers=1)
    fut = pool.submit(fn)
    try:
        result = fut.result(timeout=CATALOG_TIMEOUT_SEC)
    except FutTimeout:
        pool.shutdown(wait=False, cancel_futures=True)
        warnings.append(f"{label} query timed out after "
                        f"{CATALOG_TIMEOUT_SEC}s")
        return None
    except Exception as exc:  # noqa: BLE001
        pool.shutdown(wait=False)
        warnings.append(f"{label} query failed: {exc}")
        return None
    pool.shutdown(wait=False)
    return result


def _ang_sep_arcsec(ra1, dec1, ra2, dec2):
    """Small-angle separation in arcsec (positions in deg)."""
    dra = (ra1 - ra2) * np.cos(np.deg2rad((dec1 + dec2) / 2))
    ddec = dec1 - dec2
    return float(np.hypot(dra, ddec) * 3600.0)


def _query_gaia(ra, dec, radius_arcsec, epoch, warnings):
    """Return list of Gaia DR3 matches propagated to `epoch`."""
    out = []
    try:
        from astroquery.gaia import Gaia
        coord = SkyCoord(ra, dec, unit="deg", frame="icrs")
        rad = radius_arcsec * u.arcsec
        job = Gaia.cone_search_async(coord, radius=rad)
        tbl = job.get_results()
        def g(row, name):
            return row[name] if name in tbl.colnames else None
        for row in tbl:
            try:
                pra, pdec = _propagate(
                    float(row["ra"]), float(row["dec"]),
                    g(row, "pmra"), g(row, "pmdec"), GAIA_EPOCH_JYEAR, epoch)
            except Exception:  # noqa: BLE001
                pra, pdec = float(row["ra"]), float(row["dec"])
            sep = _ang_sep_arcsec(ra, dec, pra, pdec)
            flags = []
            try:
                ruwe = _f(g(row, "ruwe"))
                if ruwe is not None and ruwe > 1.4:
                    flags.append("RUWE>1.4 (possible unresolved binary)")
                aen = _f(g(row, "astrometric_excess_noise"))
                if aen is not None and aen > 1.0:
                    flags.append("astrometric_excess_noise>1")
                if _f(g(row, "duplicated_source")) == 1.0:
                    flags.append("duplicated_source")
            except Exception:  # noqa: BLE001
                pass
            out.append({
                "catalog": "Gaia DR3",
                "source_id": str(row["source_id"]),
                "ra_epoch": pra, "dec_epoch": pdec,
                "separation_arcsec": sep,
                "g_mean_mag": _f(g(row, "phot_g_mean_mag")),
                "parallax_mas": _f(g(row, "parallax")),
                "pmra_masyr": _f(g(row, "pmra")),
                "pmdec_masyr": _f(g(row, "pmdec")),
                "flags": flags,
            })
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"Gaia DR3 query failed: {exc}")
    out.sort(key=lambda m: m["separation_arcsec"])
    return out


def _query_allwise(ra, dec, radius_arcsec, epoch, warnings):
    """Return list of AllWISE matches propagated to `epoch`."""
    from astroquery.ipac.irsa import Irsa
    coord = SkyCoord(ra, dec, unit="deg", frame="icrs")
    out = []
    tbl = None
    for attempt in range(2):
        try:
            tbl = Irsa.query_region(coord, catalog="allwise_p3as_psd",
                                    radius=radius_arcsec * u.arcsec)
            break
        except Exception as exc:  # noqa: BLE001
            if attempt == 1:
                raise
            warnings.append(f"AllWISE attempt {attempt + 1} failed: {exc}")
            import time
            time.sleep(5)
    try:
        for row in tbl:
            try:
                pra, pdec = _propagate(
                    float(row["ra"]), float(row["dec"]),
                    row["pmra"], row["pmdec"], ALLWISE_EPOCH_JYEAR, epoch)
            except Exception:  # noqa: BLE001
                pra, pdec = float(row["ra"]), float(row["dec"])
            sep = _ang_sep_arcsec(ra, dec, pra, pdec)
            flags = []
            try:
                cc = str(row["cc_flags"]) if "cc_flags" in tbl.colnames else ""
                if any(c not in ("0", " ", "-") for c in cc):
                    flags.append(f"cc_flags={cc} (contamination/artifact)")
                ph = str(row["ph_qual"]) if "ph_qual" in tbl.colnames else ""
                if "U" in ph or "X" in ph:
                    flags.append(f"ph_qual={ph} (poor WISE photometry)")
            except Exception:  # noqa: BLE001
                pass
            out.append({
                "catalog": "AllWISE",
                "source_id": str(row["source_id"]) if "source_id" in tbl.colnames else str(row["designation"]) if "designation" in tbl.colnames else "?",
                "ra_epoch": pra, "dec_epoch": pdec,
                "separation_arcsec": sep,
                "w1mpro": _f(row["w1mpro"]) if "w1mpro" in tbl.colnames else None,
                "w2mpro": _f(row["w2mpro"]) if "w2mpro" in tbl.colnames else None,
                "pmra_masyr": _f(row["pmra"]) if "pmra" in tbl.colnames else None,
                "pmdec_masyr": _f(row["pmdec"]) if "pmdec" in tbl.colnames else None,
                "flags": flags,
            })
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"AllWISE query failed: {exc}")
    out.sort(key=lambda m: m["separation_arcsec"])
    return out


def _await(fut, label: str, warnings: list[str]):
    """Get result of a catalog future; [] + warning on timeout."""
    try:
        return fut.result(timeout=CATALOG_TIMEOUT_SEC)
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"{label} query timed out/failed: {exc}")
        return []


def _col(tbl, row, name):
    return row[name] if name in tbl.colnames else None


def _query_gaia_vizier(ra, dec, radius_arcsec, epoch, warnings):
    """Fallback: Gaia DR3 via VizieR (I/355/gaiadr3)."""
    out = []
    try:
        from astroquery.vizier import Vizier
        coord = SkyCoord(ra, dec, unit="deg", frame="icrs")
        res = Vizier.query_region(coord,
                                  radius=radius_arcsec * u.arcsec,
                                  catalog="I/355/gaiadr3")
        if not res:
            return out
        tbl = res[0]
        for row in tbl:
            try:
                pra, pdec = _propagate(
                    float(_col(tbl, row, "RA_ICRS")),
                    float(_col(tbl, row, "DE_ICRS")),
                    _col(tbl, row, "pmRA"), _col(tbl, row, "pmDE"),
                    GAIA_EPOCH_JYEAR, epoch)
            except Exception:  # noqa: BLE001
                pra = pdec = None
            if pra is None:
                continue
            sep = _ang_sep_arcsec(ra, dec, pra, pdec)
            flags = []
            ruwe = _f(_col(tbl, row, "RUWE"))
            if ruwe is not None and ruwe > 1.4:
                flags.append("RUWE>1.4 (possible unresolved binary)")
            out.append({
                "catalog": "Gaia DR3 (VizieR)",
                "source_id": str(_col(tbl, row, "Source") or "?"),
                "ra_epoch": pra, "dec_epoch": pdec,
                "separation_arcsec": sep,
                "g_mean_mag": _f(_col(tbl, row, "Gmag")),
                "parallax_mas": _f(_col(tbl, row, "Plx")),
                "pmra_masyr": _f(_col(tbl, row, "pmRA")),
                "pmdec_masyr": _f(_col(tbl, row, "pmDE")),
                "flags": flags,
            })
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"Gaia DR3 VizieR fallback failed: {exc}")
    out.sort(key=lambda m: m["separation_arcsec"])
    return out


def _query_allwise_vizier(ra, dec, radius_arcsec, epoch, warnings):
    """Fallback: AllWISE via VizieR (II/328/allwise)."""
    out = []
    try:
        from astroquery.vizier import Vizier
        coord = SkyCoord(ra, dec, unit="deg", frame="icrs")
        res = Vizier.query_region(coord,
                                  radius=radius_arcsec * u.arcsec,
                                  catalog="II/328/allwise")
        if not res:
            return out
        tbl = res[0]
        for row in tbl:
            try:
                pra, pdec = _propagate(
                    float(_col(tbl, row, "RAJ2000")),
                    float(_col(tbl, row, "DEJ2000")),
                    _col(tbl, row, "pmRA"), _col(tbl, row, "pmDE"),
                    ALLWISE_EPOCH_JYEAR, epoch)
            except Exception:  # noqa: BLE001
                pra = pdec = None
            if pra is None:
                continue
            sep = _ang_sep_arcsec(ra, dec, pra, pdec)
            flags = []
            cc = str(_col(tbl, row, "ccflg") or "")
            if cc and cc != "0":
                flags.append(f"ccflg={cc} (contamination/artifact)")
            ph = str(_col(tbl, row, "phqual") or "")
            if "U" in ph or "X" in ph:
                flags.append(f"phqual={ph} (poor WISE photometry)")
            out.append({
                "catalog": "AllWISE (VizieR)",
                "source_id": str(_col(tbl, row, "AllWISE") or "?"),
                "ra_epoch": pra, "dec_epoch": pdec,
                "separation_arcsec": sep,
                "w1mpro": _f(_col(tbl, row, "W1mag")),
                "w2mpro": _f(_col(tbl, row, "W2mag")),
                "pmra_masyr": _f(_col(tbl, row, "pmRA")),
                "pmdec_masyr": _f(_col(tbl, row, "pmDE")),
                "flags": flags,
            })
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"AllWISE VizieR fallback failed: {exc}")
    out.sort(key=lambda m: m["separation_arcsec"])
    return out


def _f(v):
    try:
        if v is None:
            return None
        f = float(v)
        return f if np.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def crossmatch(ra: float, dec: float,
               radius_arcsec: float = DEFAULT_SEARCH_RADIUS_ARCSEC,
               epoch: float = 2025.5) -> dict:
    """Query Gaia DR3 + AllWISE (parallel, with VizieR fallbacks)."""
    warnings: list[str] = []
    pool = ThreadPoolExecutor(max_workers=2)
    fg = pool.submit(_query_gaia, ra, dec, radius_arcsec, epoch, warnings)
    fw = pool.submit(_query_allwise, ra, dec, radius_arcsec, epoch, warnings)
    gaia = _await(fg, "Gaia DR3", warnings)
    wise = _await(fw, "AllWISE", warnings)
    pool.shutdown(wait=False)

    # VizieR fallbacks when the primary service returned nothing
    if not gaia:
        gv = _with_timeout(
            lambda: _query_gaia_vizier(ra, dec, radius_arcsec, epoch,
                                       warnings),
            "Gaia DR3 (VizieR)", warnings)
        gaia = gv or []
    if not wise:
        wv = _with_timeout(
            lambda: _query_allwise_vizier(ra, dec, radius_arcsec, epoch,
                                          warnings),
            "AllWISE (VizieR)", warnings)
        wise = wv or []

    matches = gaia + wise
    nearest = matches[0] if matches else None
    within = [m for m in matches
              if m["separation_arcsec"] <= MATCH_TOLERANCE_ARCSEC]

    if within:
        classification = "known_object"
        best = min(within, key=lambda m: m["separation_arcsec"])
        summary = (f"Matched {best['catalog']} source "
                   f"{best['source_id']} at {best['separation_arcsec']:.2f}\"")
    elif matches:
        classification = "unmatched_within_tolerance"
        summary = (f"{len(matches)} catalog source(s) in search radius, "
                   f"but none within {MATCH_TOLERANCE_ARCSEC}\" — "
                   f"treat as discovery CANDIDATE, not a confirmed discovery")
    else:
        classification = "discovery_candidate"
        summary = ("No catalog counterparts found — DISCOVERY CANDIDATE "
                   "only; requires follow-up confirmation")

    all_flags = sorted({fl for m in matches for fl in m["flags"]})
    result = {
        "ra": ra, "dec": dec,
        "search_radius_arcsec": radius_arcsec,
        "epoch": epoch,
        "match_tolerance_arcsec": MATCH_TOLERANCE_ARCSEC,
        "classification": classification,
        "summary": summary,
        "n_matches": len(matches),
        "matches": matches[:20],
        "flags": all_flags,
        "warnings": warnings,
        "disclaimer": ("Catalog matches indicate a likely known object; "
                       "absence of a match does NOT confirm a discovery — "
                       "it only means no counterpart in queried catalogs."),
    }
    return _catalog_cache_store_or_load(ra, dec, result, warnings)


def classify_candidate_target(candidate: dict[str, Any],
                              catalog_result: dict[str, Any]) -> str:
    """Classify a candidate using both pipeline detection and catalog cross-matching.

    Categories:
      - 'known_object': Source matches Gaia DR3 or AllWISE within tolerance.
      - 'solar_system_target': Fast-moving candidate (displacement > 0.5") with
        no optical stellar counterpart in Gaia.
      - 'high_priority_candidate': Uncataloged point source with good PSF and SNR,
        consistent with cold brown dwarf or Planet X analog.
      - 'discovery_candidate': Uncataloged source requiring human follow-up.
      - 'unmatched_within_tolerance': Nearby catalog sources exist (3-10") but none within 3".
    """
    cat_class = catalog_result.get("classification", "discovery_candidate")
    motion = candidate.get("motion", {})
    is_moving = motion.get("is_moving", False)
    disp = motion.get("displacement_arcsec", 0.0)
    psf_eval = candidate.get("psf_eval", {})
    is_point_source = psf_eval.get("is_point_source", False)
    confidence = candidate.get("confidence", 0.0)

    if cat_class == "known_object":
        return "known_object"

    if is_moving and disp >= 0.5 and not catalog_result.get("matches"):
        return "solar_system_target"

    if cat_class in ("discovery_candidate", "unmatched_within_tolerance"):
        if is_point_source and confidence >= 50.0:
            return "high_priority_candidate"
        if is_moving:
            return "solar_system_target"
        return cat_class

    return cat_class


CATALOG_CACHE_DIR = Path(__file__).resolve().parent / "cache" / "catalog_cache"
CATALOG_CACHE_MAX_AGE_SEC = 7 * 24 * 3600


def _catalog_cache_key(ra: float, dec: float) -> str:
    import hashlib
    return hashlib.sha256(f"{ra:.3f},{dec:.3f}".encode()).hexdigest()[:16]


def _catalog_cache_store_or_load(ra, dec, result, warnings):
    """Cache successful lookups; serve cache when live queries failed."""
    import json
    import time as _time
    path = CATALOG_CACHE_DIR / f"{_catalog_cache_key(ra, dec)}.json"
    services_failed = (result["n_matches"] == 0 and len(warnings) > 0)
    try:
        CATALOG_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        if services_failed and path.exists():
            cached = json.loads(path.read_text())
            age = _time.time() - path.stat().st_mtime
            if age < CATALOG_CACHE_MAX_AGE_SEC and cached.get("n_matches", 0) > 0:
                cached["warnings"] = warnings + [
                    f"Serving cached catalog result "
                    f"({int(age // 3600)}h old) — live services unavailable"]
                cached["offline_cached"] = True
                return cached
        elif result["n_matches"] > 0 and not services_failed:
            path.write_text(json.dumps(result, default=str))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"catalog cache error: {exc}")
    return result


if __name__ == "__main__":
    import json
    print(json.dumps(crossmatch(133.795129, -7.245147), indent=2, default=str))
