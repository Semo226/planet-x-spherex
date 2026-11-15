# Tech Stack — Planet X and SPHEREx

## Recommended Stack

| Layer | Choice | Rationale |
|---|---|---|
| Backend | **Python + FastAPI** | Same ecosystem as astropy/astroquery; async-friendly; easy to wrap IRSA calls |
| Data access | **astroquery.ipac.irsa** (+ `pyvo` for raw TAP if needed) | Officially documented IRSA path; no invented endpoints |
| FITS processing | **astropy.io.fits**, **reproject** (if mosaicking), **numpy** | Standard; handles SPHEREx MEFs and alternate WCS |
| Visualization data | **matplotlib** (static), FITS cutouts served to frontend; optional **firefly-client** for rich exploration | Firefly is the IRSA-recommended viewer for SPHEREx WCS |
| Frontend | **React + Vite** (TypeScript), **Aladin Lite** or custom canvas for FITS→PNG tiles, **plotly/d3** for light curves | FITS isn't browser-native — plan a server-side rasterization step |
| Deployment | Local-first for hackathon; optional Docker later | Phase 0 has no deployment requirement |

## Key Architectural Constraint

SPHEREx spectral images are **3D data** (2D spatial × ~108 spectral
channels encoded via WCS-WAVE lookup). The browser cannot read FITS
natively. Therefore the backend must:

1. Query IRSA (SIA2/TAP) for MEFs overlapping a coordinate.
2. Download (or cutout) MEFs.
3. Extract requested band/channel images → PNG or compressed arrays.
4. Serve per-epoch tiles to the frontend for the animation viewer.

Caching is essential: cache cutouts keyed by `(url, center, size, band)` in
`sample-data/` or a local cache dir (gitignored).

## Data Flow (MVP)

```
User enters RA/Dec (or picks WISE 0855 preset)
  → POST /api/query  {ra, dec, radius, band?}
  → backend: Irsa.query_sia(...) [+ query_tap fallback]
  → returns epoch list (t_min, obs_id, access_url, bandpass)
User selects epochs / hits "animate"
  → GET /api/cutout?url=...&center=...&size=...&band=...
  → backend downloads MEF cutout (cached), extracts band → PNG
  → frontend plays epoch sequence
```

## Alternatives Considered

- **All-sky cube (`spherex_dr1_cube`)**: simpler (one file, HiPS) but coarse
  spatial resolution vs. spectral images; keep as P1 option.
- **Firefly embedded viewer**: IRSA's own web viewer understands SPHEREx
  WCS — good for P1 "expert mode" link-out rather than core MVP.
