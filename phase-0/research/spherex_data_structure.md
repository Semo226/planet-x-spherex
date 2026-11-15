# SPHEREx Level 2 Spectral Image MEF Structure

Source: IRSA tutorial "Introduction to SPHEREx Spectral Images"
(https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/) and the
SPHEREx Explanatory Supplement (linked from the IRSA user guide).

## File Format

- Level 2 Spectral Image products are **multi-extension FITS (MEF)** files,
  one per detector array exposure. Verified `hdulist.info()` from the
  official notebook:

```
No.  Name      Ver   Type        Cards  Dimensions   Format
 0  PRIMARY    1     PrimaryHDU  5      ()
 1  IMAGE      1     ImageHDU    293    (2040, 2040) float32
 2  FLAGS      1     CompImageHDU 185   (2040, 2040) int32
 3  VARIANCE   1     ImageHDU    119    (2040, 2040) float32
 4  ZODI       1     ImageHDU    122    (2040, 2040) float32
 5  EPSF       1     BinTableHDU 56     441R x 11C   [...]
 6  WCS-WAVE   1     BinTableHDU 22     1R x 3C      [...]
```

## Extensions

| HDU | Content |
|-----|---------|
| `IMAGE` | Calibrated surface-brightness flux density in **MJy/sr**, 2040×2040 px. Zodiacal light NOT subtracted (see ZODI). |
| `FLAGS` | Per-pixel bitmap of status/processing flags (compressed int32). |
| `VARIANCE` | Per-pixel variance estimate (float32). |
| `ZODI` | Model of zodiacal dust background — provided separately, not subtracted from IMAGE. |
| `EPSF` | Effective PSF table: 441 spatial regions × oversampled pixels; DR1/QR3 only. |
| `WCS-WAVE` | Spectral WCS lookup table mapping pixel → central wavelength and bandwidth per spectral channel. |

## World Coordinate Systems

- **Spatial WCS** lives in the IMAGE extension header (standard FITS WCS
  keywords; SPHEREx detector arrays have non-trivial distortion — SIP
  coefficients are distributed as separate `cal-sip` products).
- **Spectral WCS** is *not* a simple linear CDELT: `WCS-WAVE` is a lookup
  table (BinTable) because each detector column maps to a non-uniform
  wavelength/bandwidth grid. Use `astropy.wcs.WCS` with the appropriate
  alternate WCS (the Firefly viewer understands SPHEREx's alternate WCS
  axes — wavelength and bandwidth axes are exposed via alternate WCS
  descriptions in the headers).

## Filename Convention (from browsable directories)

`level2_<PlanningPeriod>_<ObsIndex>_<Detector?>_spx_l2b-v<ver>-<date>.fits`

Example from docs:
`level2_2025W18_1B_0074_2D1_spx_l2b-v27-2026-196.fits`
Directory layout:
`<release>/level2/<PlanningPeriod>/l2b-v<ver>-<procdate>/<Detector>/`

## Detector Arrays

- 6 detector arrays (2×3 focal plane); filenames encode array IDs like
  `1B`, `2D1` etc. Each MEF = one array, one exposure (spectral image =
  2040 spatial × ~108 spectral channels encoded in the WCS-WAVE table).

## Practical Notes for Our Pipeline

- For a sky position, a single SIA2 hit gives one array's image; the same
  epoch appears across multiple arrays/detectors — deduplicate by
  `t_min`/`obs_id` when building epoch lists.
- Cutouts preserve all extensions (MEF cutout), so FLAGS/VARIANCE/EPSF are
  available for small regions — ideal for our animation viewer.
- Flux units MJy/sr; for point-source work use EPSF + VARIANCE.
