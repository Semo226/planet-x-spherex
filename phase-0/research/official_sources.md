# Official Sources — Phase 0 Research

All endpoints, collection names, and TAP table names below were verified
against official NASA/IPAC IRSA documentation on 2026-11-14. Nothing here is
invented; citations are inline.

## Primary Documentation

1. **SPHEREx Data Access — SPHEREx Archive at IRSA** (updated 2026-11-12)
   https://irsa.ipac.caltech.edu/spherex-archive-documentation/spherex-data-access/
   - Defines the two hosting locations (IPAC on-prem + AWS cloud mirror).
   - Documents SIA2 endpoint, collections, browsable directories, and cutout
     service query parameters (`center`, `size`).

2. **IRSA Tutorial: Introduction to SPHEREx Spectral Images**
   https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/
   - Shows `astroquery.ipac.irsa.Irsa.query_sia` usage against
     `['spherex_dr1', 'spherex_qr3']`.
   - Lists the SIA2 result column names (IVOA-standard).
   - Documents the MEF HDU structure (IMAGE / FLAGS / VARIANCE / ZODI /
     EPSF / WCS-WAVE) and units (MJy/sr, 2040×2040 per detector array).

3. **IRSA Tutorial: Download a collection of SPHEREx Spectral Image cutouts
   as a multi-extension FITS file**
   https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/
   - Shows `Irsa.query_tap` (ADQL) against `spherex.artifact` JOIN
     `spherex.plane` with `CONTAINS(POINT(...), p.poly)` spatial indexing,
     and constructs `?center=&size=` cutout URLs.
   - Notes TAP cutouts bypass the ~1-day SIA2 ingestion lag.

4. **IRSA SIA v2 Queries (IBE)**
   https://irsa.ipac.caltech.edu/ibe/sia.html
   - IRSA's generic SIA2 endpoint: `https://irsa.ipac.caltech.edu/SIA?`
   - Supported keywords include `COLLECTION`, `POS=circle RA DEC radius`,
     `RESPONSEFORMAT`, `MAXREC`, `BAND`, `TIME`.

5. **astroquery IPAC/IRSA module docs**
   https://docs.astropy.org/projects/astroquery/en/latest/ipac/irsa/irsa.html
   - `Irsa.list_collections(servicetype='SIA', filter='spherex')`
   - `Irsa.query_sia(pos=..., collection=...)`
   - `Irsa.query_tap(adql)`

6. **Using IRSA TAP services from Python**
   https://irsa.ipac.caltech.edu/docs/program_interface/astropy_TAP.html
   - Confirms `query_tap` is the recommended path for complex ADQL;
     `query_region` for simple cone/box/polygon queries.

## SPHEREx Mission Context

- SPHEREx (Spectro-Photometer for the History of the Universe, Epoch of
  Reionization, and Ices Explorer) launched March 2025; 2-year mission;
  0.75–5 µm all-sky spectroscopy. Source: spherex-intro notebook §2.
- Data Release 1 (DR1) and Quick Release 3 (QR3) are both live at IRSA as of
  2026-11; QR2 remains downloadable but is NOT searchable via SIA2. DR1 and
  QR3 supersede QR2 — always query DR1+QR3.

## Confirmed SIA2 Collections (SPHEREx)

| Collection | Contents |
|---|---|
| `spherex_dr1` | DR1 Spectral Image MEFs, Wide Survey |
| `spherex_dr1_deep` | DR1 Spectral Image MEFs, Deep Survey |
| `spherex_dr1_cube` | DR1 All-Sky Cube FITS files |
| `spherex_qr3` | QR3 Spectral Image MEFs, Wide Survey |
| `spherex_qr3_deep` | QR3 Spectral Image MEFs, Deep Survey |

## TAP Tables (from cutouts notebook)

- `spherex.artifact` — file-level metadata (includes `uri`, `planeid`)
- `spherex.plane` — plane-level metadata incl. `time_bounds_lower`,
  `energy_bandpassname`, and sky footprint polygon `p.poly` (ICRS)

## Cutout Service

- Append `?center=<RAdeg>,<Decdeg>&size=<size>` to any on-prem Spectral Image
  MEF URL, e.g.:
  `.../level2_....fits?center=337.8783221,-5.1891588&size=0.1`
  (size units per IRSA cutouts doc: degrees; verify in experiment).
- Docs: https://irsa.ipac.caltech.edu/ibe/cutouts.html

## Ingestion Cadence

- Weekly ingestion; SIA2 availability lags ~1 day behind browsable
  directories / TAP queries. Query results grow over time.
