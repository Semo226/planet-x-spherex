# IRSA Access Methods for SPHEREx

Verified against official IRSA docs (see `official_sources.md`).

## 1. astroquery (recommended for Python)

```python
from astroquery.ipac.irsa import Irsa

Irsa.list_collections(servicetype='SIA', filter='spherex')
results = Irsa.query_sia(pos=(coord, radius), collection=['spherex_dr1', 'spherex_qr3'])
results = Irsa.query_tap(adql_string)   # ADQL via TAP
```

- `query_sia` returns an astropy Table with IVOA SIA2 columns
  (`s_ra`, `s_dec`, `access_url`, `t_min`, `t_max`, `energy_bandpassname`,
  `obs_id`, `s_region`, `cloud_access`, ...).
- `query_tap` runs arbitrary ADQL against IRSA's TAP service.

## 2. SIA2 (Simple Image Access v2)

- Endpoint: `https://irsa.ipac.caltech.edu/SIA?`
- Query params: `COLLECTION`, `POS=circle <ra> <dec> <radius_deg>`,
  `RESPONSEFORMAT` (table/votable/json), `MAXREC`, `BAND`, `TIME`.
- Returns file-level metadata + `access_url` per image.
- Lag: new data appear in SIA2 ~1 day after weekly ingestion.

## 3. TAP (Table Access Protocol)

- Endpoint handled by `Irsa.query_tap`; schema `spherex` contains
  `artifact` (files, `uri`, `planeid`) and `plane` (per-plane metadata:
  `time_bounds_lower`, `energy_bandpassname`, `poly` sky footprint).
- Spatial indexing via `CONTAINS(POINT('ICRS', ra, dec), p.poly)`.
- Advantage: no SIA2 ingestion lag; needed for bulk cutout URL lists.

## 4. Cutout Service

- For any on-prem MEF URL: append `?center=<ra>,<dec>&size=<deg>`.
- Returns a cutout MEF preserving extensions.

## 5. Browsable Directories / AWS

- Root: `https://irsa.ipac.caltech.edu/ibe/data/spherex/`
- Mirror on AWS S3 for bulk access (`cloud_access` column in SIA2 results).

## Which to Use for Our App

| Need | Method |
|---|---|
| "Does anything cover coord X?" | `Irsa.query_sia` (simple, fast) |
| Epoch list + cutout URLs for animation | TAP (`spherex.artifact` ⨝ `spherex.plane`) |
| Fresh data (<1 day old) | TAP or browsable dirs |
| Quick interactive check | IRSA SPHEREx Data Explorer GUI |
