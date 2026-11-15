# ADR-0001: Difference Imaging Approach (Phase 3)

Date: 2026-11-14 · Status: accepted (with caveats)

## Decision

- Alignment: FFT phase-correlation (`processing/diff_imaging.py`),
  integer+subpixel shift, with a **2 px plausibility gate** — shifts larger
  than that are rejected and the frame is left unaligned.
- Static suppression: per-pixel **median template** subtraction over a
  **single detector array** (auto-selected as the most common detector in
  the session; `detector` param on `POST /api/process`).
- Residuals rendered with diverging `RdBu_r` colormap, symmetric
  robust (MAD-based) limits.

## Why per-detector

SPHEREx's six detector arrays have independent distortion solutions.
Measured phase-correlation offsets between arrays reach ~13 px on 29×29
cutouts — cross-detector stacks are not alignable by translation.

## Measured results (WISE 0855, D5 stack, 27 frames)

- Phase-correlation shifts: median ~8 px, max ~13 px → nearly all rejected
  by the 2 px gate. The cutouts are already astrometrically registered by
  IRSA's cutout service, and the field is faint/zodi-dominated, so the
  correlation peak is noise-driven.
- Median subtraction still suppresses static structure: residual std ≈
  77% of raw std for D5 frames (modest — single-epoch 30″ cutouts of a
  faint field have limited static-source SNR).

## Consequences / honest caveats

- "Processed" view currently ≈ unaligned median-subtracted frames. True
  mover highlighting requires either larger FOV (more stars for
  registration), PSF-matched subtraction, or image stacking — all P1+.
- WISE 0855's own motion is sub-pixel at these epochs/cutout centering, so
  it will NOT stand out as a mover in difference frames; residuals show
  noise, coverage edges, and artifacts. This is reported honestly in the
  UI/status, not hidden.
- `align` flag remains in the API for future improvement (e.g.,
  template-based iterative registration or brighter-source masks).
