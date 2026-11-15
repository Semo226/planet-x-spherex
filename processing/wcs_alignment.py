"""WCS-based astronomical image alignment and celestial registration for SPHEREx.

Phase 5 module: Replaces noise-sensitive pixel cross-correlation with rigorous
astrometric coordinate transformation.

Pipeline:
  1. Reads multi-epoch SPHEREx FITS cutouts from local cache or storage.
  2. Extracts WCS information (including SIP polynomial distortion terms) from
     the IMAGE HDU headers.
  3. Establishes a common celestial reference grid from a chosen reference epoch
     (or auto-selected primary frame).
  4. Reprojects secondary epoch frames onto the reference frame's celestial grid
     using `reproject.reproject_interp` so that pixel coordinates map to identical
     sky positions (RA/Dec).
  5. Computes overlap footprints, validity masks, and astrometric metadata.
  6. Returns aligned image arrays ready for immediate consumption by the
     frame-differencing pipeline (`processing/diff_imaging.py`).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

try:
    from reproject import reproject_exact, reproject_interp
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "reproject is required for WCS alignment: pip install reproject"
    ) from exc

logger = logging.getLogger(__name__)


def read_fits_wcs(fits_path: str | Path,
                  hdu_name: str = "IMAGE") -> tuple[np.ndarray, WCS, fits.Header]:
    """Read data array, WCS, and header from a SPHEREx MEF cutout.

    Parameters
    ----------
    fits_path : str or Path
        Path to the FITS file.
    hdu_name : str, default "IMAGE"
        Name of the extension containing the 2D image and spatial WCS.

    Returns
    -------
    data : np.ndarray
        2D float64 image array.
    wcs : astropy.wcs.WCS
        Celestial World Coordinate System object (including SIP distortion).
    header : astropy.io.fits.Header
        Header cards of the image extension.

    Raises
    ------
    FileNotFoundError
        If `fits_path` does not exist.
    ValueError
        If no matching HDU is found or data is not 2-dimensional.
    """
    fits_path = Path(fits_path)
    if not fits_path.exists():
        raise FileNotFoundError(f"FITS file not found: {fits_path}")

    with fits.open(fits_path) as hdul:
        target_hdu = None
        for hdu in hdul:
            if (hdu.name or "").upper() == hdu_name.upper() and hdu.data is not None:
                target_hdu = hdu
                break
        if target_hdu is None:
            raise ValueError(f"No HDU named '{hdu_name}' with data found in {fits_path}")

        data = np.asarray(target_hdu.data, dtype=np.float64)
        if data.ndim != 2:
            raise ValueError(
                f"Expected 2D image in HDU '{hdu_name}', got shape {data.shape} in {fits_path}"
            )
        header = target_hdu.header.copy()
        wcs = WCS(header)

    return data, wcs, header


def read_fits_bundle(fits_path: str | Path,
                     hdu_name: str = "IMAGE") -> dict[str, Any]:
    """Read all relevant Level 2 extensions (IMAGE, FLAGS, VARIANCE, EPSF).

    Returns a dict with 'data', 'wcs', 'header', 'flags', 'variance', 'epsf_info'.
    """
    fits_path = Path(fits_path)
    if not fits_path.exists():
        raise FileNotFoundError(f"FITS file not found: {fits_path}")

    with fits.open(fits_path) as hdul:
        target_hdu = None
        flags_data = None
        variance_data = None
        epsf_info = None

        for hdu in hdul:
            name = (hdu.name or "").upper()
            if name == hdu_name.upper() and hdu.data is not None:
                target_hdu = hdu
            elif name == "FLAGS" and hdu.data is not None:
                flags_data = np.asarray(hdu.data, dtype=np.int32)
            elif name == "VARIANCE" and hdu.data is not None:
                variance_data = np.asarray(hdu.data, dtype=np.float64)
            elif name == "EPSF" and hdu.data is not None:
                try:
                    epsf_info = {
                        "neff_mean": float(hdu.data["NEFF_MEAN"][0]),
                        "cwave": float(hdu.data["CWAVE"][0]),
                    }
                except Exception:  # noqa: BLE001
                    epsf_info = None

        if target_hdu is None:
            raise ValueError(f"No HDU named '{hdu_name}' with data found in {fits_path}")

        data = np.asarray(target_hdu.data, dtype=np.float64)
        header = target_hdu.header.copy()
        wcs = WCS(header)

    return {
        "data": data,
        "wcs": wcs,
        "header": header,
        "flags": flags_data,
        "variance": variance_data,
        "epsf_info": epsf_info,
    }


def get_wcs_metadata(wcs: WCS, shape: tuple[int, int]) -> dict[str, Any]:
    """Extract diagnostic celestial metadata from a WCS object.

    Parameters
    ----------
    wcs : WCS
        Celestial WCS object.
    shape : tuple of (int, int)
        Image dimensions (ny, nx).

    Returns
    -------
    dict
        Metadata including center (RA/Dec), pixel scales, projection type.
    """
    ny, nx = shape
    center_x = (nx - 1) / 2.0
    center_y = (ny - 1) / 2.0

    center_coord = None
    try:
        center_coord = wcs.pixel_to_world(center_x, center_y)
        center_ra = float(center_coord.ra.deg)
        center_dec = float(center_coord.dec.deg)
    except Exception:  # noqa: BLE001
        center_ra, center_dec = None, None

    pixel_scale_arcsec = None
    try:
        scales = proj_plane_pixel_scales(wcs)
        pixel_scale_arcsec = float(np.mean(scales) * 3600.0)
    except Exception:  # noqa: BLE001
        pass

    return {
        "center_ra_deg": center_ra,
        "center_dec_deg": center_dec,
        "pixel_scale_arcsec": pixel_scale_arcsec,
        "ctype": list(wcs.wcs.ctype) if hasattr(wcs, "wcs") else [],
        "shape": (ny, nx),
    }


def reproject_frame_to_wcs(data: np.ndarray,
                           input_wcs: WCS,
                           ref_wcs: WCS,
                           shape_out: tuple[int, int],
                           method: str = "interp") -> tuple[np.ndarray, np.ndarray]:
    """Reproject a single 2D image onto a reference celestial WCS grid.

    Parameters
    ----------
    data : np.ndarray
        2D source image data.
    input_wcs : WCS
        World Coordinate System of the source image.
    ref_wcs : WCS
        Target World Coordinate System grid.
    shape_out : tuple of (int, int)
        Target shape (ny, nx).
    method : str, default "interp"
        Reprojection algorithm: "interp" (bilinear via reproject_interp) or
        "exact" (flux-conserving spherical polygon intersection via reproject_exact).

    Returns
    -------
    reprojected : np.ndarray
        2D array of shape `shape_out` matching the reference celestial grid.
        Uncovered / non-overlapping pixels are NaN.
    footprint : np.ndarray
        2D array of shape `shape_out` containing the fractional overlap (0.0 to 1.0).
    """
    if method == "exact":
        reprojected, footprint = reproject_exact(
            (data, input_wcs),
            ref_wcs,
            shape_out=shape_out,
            return_footprint=True,
        )
    else:
        reprojected, footprint = reproject_interp(
            (data, input_wcs),
            ref_wcs,
            shape_out=shape_out,
            return_footprint=True,
        )

    # Standardize non-overlapping/masked regions as NaN
    reprojected[footprint == 0] = np.nan
    return reprojected.astype(np.float64), footprint.astype(np.float64)


def align_fits_stack_wcs(fits_paths: list[str | Path],
                         ref_index: int = 0,
                         method: str = "interp") -> dict[str, Any]:
    """Align a multi-epoch stack of SPHEREx FITS cutouts onto a common WCS grid.

    Parameters
    ----------
    fits_paths : list of (str or Path)
        Chronological list of FITS cutout paths.
    ref_index : int, default 0
        Index of the frame used to establish the celestial reference grid.
        If the specified frame cannot be read or is invalid, the first valid
        frame is chosen automatically.
    method : str, default "interp"
        Reprojection interpolation method ('interp' or 'exact').

    Returns
    -------
    dict with keys:
        - "aligned_images": list of np.ndarray (shape_out), precisely mapped to ref WCS.
        - "footprints": list of np.ndarray (shape_out), overlap coverage per pixel.
        - "ref_wcs": astropy.wcs.WCS of the celestial reference grid.
        - "ref_header": Header cards of the reference frame.
        - "ref_index": Selected reference index.
        - "shape_out": Output image dimensions (ny, nx).
        - "metadata": List of per-frame alignment metadata dicts.
        - "n_frames": Total number of input frames.
        - "n_aligned": Number of successfully registered frames.
    """
    if not fits_paths:
        raise ValueError("fits_paths list cannot be empty")

    n_total = len(fits_paths)
    loaded_data: list[np.ndarray | None] = [None] * n_total
    loaded_wcs: list[WCS | None] = [None] * n_total
    loaded_headers: list[fits.Header | None] = [None] * n_total
    errors: list[str | None] = [None] * n_total

    # 1. Read all FITS headers and WCS objects
    for i, p in enumerate(fits_paths):
        try:
            d, w, h = read_fits_wcs(p, hdu_name="IMAGE")
            if not np.isfinite(d).any():
                raise ValueError("Frame contains only NaNs or non-finite values")
            loaded_data[i] = d
            loaded_wcs[i] = w
            loaded_headers[i] = h
        except Exception as exc:  # noqa: BLE001
            errors[i] = str(exc)
            logger.warning("[wcs_alignment] Failed reading frame %d (%s): %s", i, p, exc)

    # 2. Select primary reference frame
    valid_indices = [i for i in range(n_total) if loaded_data[i] is not None]
    if not valid_indices:
        raise ValueError("No readable frames with valid 2D image data found in stack")

    chosen_ref = ref_index if ref_index in valid_indices else valid_indices[0]
    ref_wcs = loaded_wcs[chosen_ref]
    ref_header = loaded_headers[chosen_ref]
    ref_data = loaded_data[chosen_ref]
    assert ref_wcs is not None and ref_data is not None
    shape_out = ref_data.shape

    # 3. Reproject all frames onto the reference celestial grid
    aligned_images: list[np.ndarray] = []
    footprints: list[np.ndarray] = []
    metadata: list[dict[str, Any]] = []

    for i in range(n_total):
        p = Path(fits_paths[i])
        if loaded_data[i] is None or loaded_wcs[i] is None:
            # Preserve array indexing with NaN plane
            nan_plane = np.full(shape_out, np.nan, dtype=np.float64)
            zero_fp = np.zeros(shape_out, dtype=np.float64)
            aligned_images.append(nan_plane)
            footprints.append(zero_fp)
            metadata.append({
                "index": i,
                "path": str(p),
                "valid": False,
                "is_reference": (i == chosen_ref),
                "error": errors[i] or "Unknown loading error",
                "overlap_fraction": 0.0,
                "finite_pixels": 0,
            })
            continue

        if i == chosen_ref:
            # The reference frame is already on the reference grid
            aligned_img = ref_data.copy()
            fp = np.ones(shape_out, dtype=np.float64)
            overlap_frac = 1.0
        else:
            try:
                aligned_img, fp = reproject_frame_to_wcs(
                    loaded_data[i],
                    loaded_wcs[i],
                    ref_wcs=ref_wcs,
                    shape_out=shape_out,
                    method=method,
                )
                overlap_frac = float(np.sum(fp > 0) / fp.size)
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "[wcs_alignment] Reprojection failed for frame %d (%s): %s",
                    i, p, exc
                )
                aligned_img = np.full(shape_out, np.nan, dtype=np.float64)
                fp = np.zeros(shape_out, dtype=np.float64)
                overlap_frac = 0.0
                errors[i] = str(exc)

        finite_count = int(np.isfinite(aligned_img).sum())
        aligned_images.append(aligned_img)
        footprints.append(fp)
        metadata.append({
            "index": i,
            "path": str(p),
            "valid": (finite_count > 0),
            "is_reference": (i == chosen_ref),
            "error": errors[i],
            "overlap_fraction": overlap_frac,
            "finite_pixels": finite_count,
        })

    grid_meta = get_wcs_metadata(ref_wcs, shape_out)

    return {
        "aligned_images": aligned_images,
        "footprints": footprints,
        "ref_wcs": ref_wcs,
        "ref_header": ref_header,
        "ref_index": chosen_ref,
        "shape_out": shape_out,
        "metadata": metadata,
        "grid_metadata": grid_meta,
        "n_frames": n_total,
        "n_aligned": sum(1 for m in metadata if m["valid"]),
    }
