"""Phase 0 Experiment: WISE 0855 SPHEREx feasibility query.

Target: WISE 0855 (Y-dwarf), ICRS RA = 133.795129 deg, Dec = -7.245147 deg.

Does:
  1. Lists available IRSA SIA2 collections matching 'spherex'.
  2. Queries SIA2 for spectral-image MEFs overlapping the target position
     across DR1 + QR3 (wide + deep).
  3. Records number of hits, observation time range, access URLs.
  4. Optionally downloads ONE cutout FITS for inspection
     (--download flag), then reports HDU structure, shapes, and header keys.
  5. Writes phase-0/experiments/query_results.csv.

Run from repo root:
    .venv\\Scripts\\python.exe phase-0/experiments/wise0855_query.py --download
"""

import argparse
import csv
import sys
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.ipac.irsa import Irsa

RA = 133.795129
DEC = -7.245147
SEARCH_RADIUS = 30 * u.arcsec  # generous for a first feasibility check
COLLECTIONS = [
    "spherex_dr1",
    "spherex_dr1_deep",
    "spherex_qr3",
    "spherex_qr3_deep",
]
OUT_CSV = Path(__file__).with_name("query_results.csv")
CUTOUT_DIR = Path(__file__).with_name("cutouts")


def query_sia(coord, radius):
    print(f"Querying IRSA SIA2 at RA={coord.ra.deg:.6f} Dec={coord.dec.deg:.6f} "
          f"r={radius.to(u.arcsec):.1f}\n  collections={COLLECTIONS}")
    return Irsa.query_sia(pos=(coord, radius), collection=COLLECTIONS)


def summarize(results):
    print(f"\nRows returned: {len(results)}")
    if len(results) == 0:
        return
    # Distinct epochs by observation start time
    tmin = results["t_min"]
    tmax = results["t_max"]
    bands = sorted(set(str(b) for b in results["energy_bandpassname"]))
    obs_ids = sorted(set(str(o) for o in results["obs_id"]))
    print(f"  time range (t_min): {min(tmin)} .. {max(tmin)}")
    print(f"  distinct obs_id count: {len(obs_ids)}")
    print(f"  bandpasses: {bands}")
    print("\nFirst 5 rows (selected columns):")
    cols = [c for c in ["obs_id", "t_min", "t_max", "energy_bandpassname",
                        "obs_collection", "access_url"] if c in results.colnames]
    for row in results[cols][:5]:
        print("   ", {c: row[c] for c in cols})


def write_csv(results):
    if len(results) == 0:
        OUT_CSV.write_text("no_results\n")
        return
    keep = [c for c in results.colnames]
    with OUT_CSV.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(keep)
        for row in results:
            writer.writerow(["" if row[c] is None else row[c] for c in keep])
    print(f"\nWrote {len(results)} rows -> {OUT_CSV}")


def download_first_cutout(results):
    """Download a small cutout of the first MEF via IRSA's cutout service."""
    if len(results) == 0:
        print("No results — nothing to download.")
        return
    url = results["access_url"][0]
    cutout_url = f"{url}?center={RA},{DEC}&size=0.05"
    print(f"\nDownloading cutout:\n  {cutout_url}")
    CUTOUT_DIR.mkdir(exist_ok=True)
    out = CUTOUT_DIR / "wise0855_first_cutout.fits"
    from astropy.utils.data import download_file
    download_file(cutout_url, cache=True)
    # download_file caches; copy to our dir for a stable artifact
    import shutil
    cached = download_file(cutout_url, cache=True)
    shutil.copy(cached, out)
    size_mb = out.stat().st_size / 1e6
    print(f"  saved {out} ({size_mb:.2f} MB)")
    inspect_fits(out)


def inspect_fits(path):
    print("\nFITS structure:")
    with fits.open(path) as hdul:
        hdul.info()
        for i, hdu in enumerate(hdul):
            hdr = hdu.header
            print(f"\n--- HDU {i}: {hdu.name or 'UNNAMED'} ---")
            for key in ("EXTNAME", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3",
                        "CTYPE1", "CTYPE2", "CTYPE3", "CUNIT", "BUNIT",
                        "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
                        "BUNIT3", "CTYPE4", "WAVE", "WAVERANGE"):
                if key in hdr:
                    print(f"  {key:10s} = {hdr[key]}")
            # print a couple of alternate WCS / spectral keywords if present
            for key in hdr:
                if key.startswith("WCSAXES") or key.startswith("CTYPE") \
                        or key.startswith("CDELT"):
                    pass  # already covered partially above


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true",
                        help="download one cutout FITS for inspection")
    args = parser.parse_args()

    coord = SkyCoord(RA, DEC, unit="deg", frame="icrs")
    try:
        cols = Irsa.list_collections(servicetype="SIA", filter="spherex")
        print("SIA collections matching 'spherex':")
        for c in cols:
            print("  ", c)
    except Exception as exc:  # noqa: BLE001
        print(f"(list_collections failed: {exc})")

    results = query_sia(coord, SEARCH_RADIUS)
    summarize(results)
    write_csv(results)

    if args.download:
        download_first_cutout(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
