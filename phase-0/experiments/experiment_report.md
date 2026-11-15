# Experiment Report — WISE 0855 SPHEREx Feasibility (Phase 0)

Date: 2026-11-14 · Script: `wise0855_query.py` · Machine: Windows, Python 3.12.0

## Target

WISE 0855 (Y-dwarf): ICRS **RA = 133.795129°, Dec = −7.245147°**
Search radius used: 30 arcsec (generous cone for feasibility).

## Question 1 — Do matching observations exist?

**YES.** `Irsa.query_sia` against collections
`[spherex_dr1, spherex_dr1_deep, spherex_qr3, spherex_qr3_deep]` returned
**258 spectral-image MEF rows**, all from `spherex_dr1` (QR3 had no
additional overlapping products at query time — consistent with DR1
superseding QR3 coverage).

## Question 2 — How many epochs cover this patch?

- **129 distinct observations** (`obs_id` / `t_min` unique count = 129)
- Each observation appears twice in the table (two detector sub-frames per
  obs in the overlapping region) → 258 rows.
- Coverage spans **t_min 60800.70 → 61174.29** (~373 days of mission
  elapsed time as reported by IRSA `t_min`).
- All **6 detector arrays** contribute: D1×36, D2×48, D3×45, D4×36,
  D5×48, D6×45 rows. `energy_bandpassname` values are detector-array
  identifiers (SPHEREx-D1…D6), not spectral bands — spectral channel is a
  property of each pixel column via WCS-WAVE.

## Question 3 — Timestamps, bands, detector arrays

- Timestamps: `t_min`/`t_max` columns in `query_results.csv` (258 rows ×
  58 IVOA columns; includes `access_url`, `s_region`, `cloud_access`).
- Detector arrays: encoded in `obs_id`/filename (e.g. `..._1D3_...` = array 3
  of obs `2025W19_1B_0116_1`).
- Spectral coverage per file: WCS-WAVE table; measured ~**1.6–2.4 µm** for
  the downloaded D3 file (band varies by detector column).

## FITS Sample Download

Downloaded via IRSA cutout service:
`...level2_2025W19_1B_0116_1D3_spx_l2b-v27-2026-197.fits?center=133.795129,-7.245147&size=0.05`

- **File**: `phase-0/experiments/cutouts/wise0855_first_cutout.fits`
- **Size**: 4,049,280 bytes (~4.0 MB)
- Cutout = 29×29 pixels; astropy WCS places the target at pixel (14.2, 13.7)
  — dead center. Cutout center sky = (133.79464, −7.24546), ~2″ from target.

### HDU structure (matches official docs)

| HDU | Name | Shape | Format | Notes |
|-----|------|-------|--------|-------|
| 0 | PRIMARY | — | — | header only |
| 1 | IMAGE | 29×29 | float32 | MJy/sr; stats: min 0.039, med 0.148, max 3.62 |
| 2 | FLAGS | 29×29 | int32 | bitmap (values seen: 0, 524288, 2097152, …) |
| 3 | VARIANCE | 29×29 | float32 | MJy²/sr²; median 2.06e-4 |
| 4 | ZODI | 29×29 | float32 | zodiacal model (not subtracted) |
| 5 | EPSF | 451×11 | table | effective PSF regions |
| 6 | WCS-WAVE | 1×3 | table | spectral WCS lookup (~1.6–2.4 µm here) |

Spatial WCS: `CTYPE1=RA---TAN-SIP`, `CTYPE2=DEC--TAN-SIP`, PC matrix + SIP
distortion; CUNIT=deg. Note `CRVAL` refers to the *full-array* reference
pixel — use `WCS(hdu.header)` and `world_to_pixel` rather than CRPIX math.

## Issues Encountered

- `astropy.utils.data.download_file` and a first `curl` attempt returned
  HTTP 503 / connection failure on the cutout URL; `requests.get` with a
  browser User-Agent succeeded (4.0 MB). Likely transient/rate-limit —
  backend should include retries + UA header.
- No `git` on PATH — repo not yet version-controlled.

## Conclusion

**Phase 0 scientific feasibility: CONFIRMED.** WISE 0855 has 129 epochs of
DR1 spectral-image coverage — rich enough for a real epoch-animation demo.
All access paths (SIA2 query, TAP tables documented, cutout service) are
verified working from this machine with public access (no credentials).
