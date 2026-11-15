"""Rebuild backend/cache/session_manifest.json from cached cutout FITS.

Used to seed the offline fallback when IRSA is unreachable: each cutout
MEF header carries DATE-AVG/MJD-AVG, OBSID and DETECTOR, which is all the
viewer needs.
"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "backend" / "cache"
MANIFEST = CACHE / "session_manifest.json"
RA, DEC, FOV = 133.795129, -7.245147, 0.05


def mjd_to_iso(mjd):
    try:
        return Time(float(mjd), format="mjd").iso
    except Exception:
        return None


def main():
    frames = []
    for fits_path in sorted(CACHE.glob("*.fits")):
        png_path = fits_path.with_suffix(".png")
        if not png_path.exists():
            continue
        try:
            with fits.open(fits_path) as hdul:
                hdr = hdul[1].header
            mjd_avg = float(hdr["MJD-AVG"])
            obsid = str(hdr.get("OBSID", fits_path.stem))
            det = int(hdr.get("DETECTOR", 0))
        except Exception as exc:
            print(f"skip {fits_path.name}: {exc}")
            continue
        frames.append({
            "obs_id": obsid,
            "bandpass": f"SPHEREx-D{det}" if det else "",
            "detector": f"D{det}" if det else "",
            "t_min": mjd_avg,
            "t_max": mjd_avg,
            "t_min_iso": mjd_to_iso(mjd_avg),
            "fits_path": str(fits_path),
            "png_path": str(png_path),
            "access_url": "",
            "cutout_url": "",
        })

    frames.sort(key=lambda f: f["t_min"])
    for i, f in enumerate(frames):
        f["index"] = i

    manifest = {
        "ra": RA, "dec": DEC, "fov_size": FOV,
        "n_results": len(frames) * 2,
        "n_epochs": len({f["obs_id"] for f in frames}),
        "query_time": "reconstructed-from-cache",
        "frames": frames,
    }
    MANIFEST.write_text(json.dumps(manifest))
    print(f"wrote {MANIFEST} with {len(frames)} frames, "
          f"{manifest['n_epochs']} unique obs_ids")
    print("range:", frames[0]["t_min_iso"], "->", frames[-1]["t_min_iso"])


if __name__ == "__main__":
    sys.exit(main())
