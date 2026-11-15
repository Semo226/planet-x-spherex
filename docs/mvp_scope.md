# MVP Scope — Planet X and SPHEREx

Hackathon cut-list. "Planet X" framing: find epoch-to-epoch variability /
unexpected sources in SPHEREx time-domain data, with WISE 0855 as proof.

## P0 — Must Have (core demo)

1. **Coordinate input**: user enters RA/Dec (degrees or `hh:mm:ss`), with a
   **WISE 0855 preset button** (RA=133.795129, Dec=-7.245147).
2. **Data fetch**: backend queries IRSA SIA2 for `spherex_dr1` +
   `spherex_qr3` collections; returns list of overlapping spectral-image
   MEFs with timestamps.
3. **Animation viewer**: play through epochs of cutout images at the target
   position (band-selectable). Even grayscale JPEG/PNG tiles are fine.
4. **WISE 0855 test**: the demo must run end-to-end on WISE 0855 and show
   whatever coverage actually exists (including "N epochs found").
5. **Honest empty states**: if a position has zero coverage, say so clearly —
   that IS a scientific result for "Planet X" (absence of data).

## P1 — Should Have (if time permits)

- **Source/catalog overlay**: query IRSA catalogs (e.g., AllWISE, 2MASS) for
  known sources in the FOV to contextualize detections.
- **FOV footprint visualization**: draw `s_region` polygons on a sky map.
- **URL params / shareable links**: `?ra=..&dec=..&band=..` deep-linking.
- **Band/wavelength selector**: use WCS-WAVE to map spectral channel ↔ µm.
- **Light-curve view**: extract flux at position per epoch → matplotlib/Plotly.

## P2 — Nice to Have (stretch)

- **AI agent**: natural-language query ("show me everything near WISE 0855
  in QR3") that drives the API.
- **Mosaic / difference imaging** between epochs.
- **AWS cloud access** path (`cloud_access` column) for speed.
- **Export**: download cutout MEF bundle for a position.

## Explicitly Out of Scope (Phase 0)

- No frontend/backend feature code yet.
- No source detection / PSF photometry pipeline.
- No all-sky processing.
