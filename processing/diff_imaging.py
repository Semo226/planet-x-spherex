"""Alignment + difference imaging for SPHEREx cutout stacks.

Pipeline:
  1. Load SPHEREx cutout MEFs and extract WCS / pixel data.
  2. Align secondary frames onto the reference frame's celestial grid
     (Phase 5: astropy WCS reprojection via `processing.wcs_alignment`;
      Phase 3 fallback: FFT phase correlation translation).
  3. Build a median "static sky" template from valid registered frames.
  4. Residual_i = aligned_i - template  -> suppresses static celestial sources,
     leaving transient events, variable sources, noise, and true movers.
  5. Render residual PNGs with a diverging RdBu_r colormap and symmetric robust limits.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits

from processing.render import _asinh_stretch  # reuse robust stretch
from processing.wcs_alignment import align_fits_stack_wcs

logger = logging.getLogger(__name__)


def load_image(fits_path: str | Path) -> np.ndarray:
    """Return the IMAGE HDU data as float64 2-D array."""
    with fits.open(fits_path) as hdul:
        for hdu in hdul:
            if (hdu.name or "").upper() == "IMAGE" and hdu.data is not None:
                return np.asarray(hdu.data, dtype=np.float64)
    raise ValueError(f"no IMAGE HDU in {fits_path}")


def _phase_shift(ref: np.ndarray, img: np.ndarray) -> tuple[float, float]:
    """Return (dy, dx) shift that aligns img onto ref, via FFT
    cross-correlation peak with quadratic subpixel refinement."""
    # Edge-taper with a Tukey-ish window to reduce boundary effects
    wy = np.hanning(ref.shape[0])[:, None]
    wx = np.hanning(ref.shape[1])[None, :]
    win = wy * wx

    with np.errstate(all="ignore"):
        a = np.nan_to_num(ref - np.nanmean(ref)) * win
        b = np.nan_to_num(img - np.nanmean(img)) * win
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return 0.0, 0.0
    A = np.fft.rfft2(a)
    B = np.fft.rfft2(b)
    # Cross-power spectrum; phase correlation
    R = B * np.conj(A)
    R /= (np.abs(R) + 1e-12)
    corr = np.fft.irfft2(R, s=ref.shape)
    corr = np.fft.fftshift(corr)

    iy, ix = np.unravel_index(np.argmax(corr), corr.shape)
    cy, cx = np.array(corr.shape) // 2

    # Quadratic subpixel refinement on the correlation peak
    def subpix(c_m, c_0, c_p):
        denom = (c_m - 2 * c_0 + c_p)
        if abs(denom) < 1e-12:
            return 0.0
        d = 0.5 * (c_m - c_p) / denom
        return float(np.clip(d, -1.0, 1.0))

    iy0, ix0 = iy, ix
    dy = subpix(corr[(iy - 1) % corr.shape[0], ix], corr[iy, ix],
                corr[(iy + 1) % corr.shape[0], ix])
    dx = subpix(corr[iy, (ix - 1) % corr.shape[1]], corr[iy, ix],
                corr[iy, (ix + 1) % corr.shape[1]])

    # corr peak at (iy,ix) means img is shifted by (iy-cy, ix-cx) rel. to ref
    shift_y = (iy - cy) + dy
    shift_x = (ix - cx) + dx
    return shift_y, shift_x


def shift_image(img: np.ndarray, dy: float, dx: float) -> np.ndarray:
    """Shift img by (dy, dx) pixels (positive dy moves data down) using
    bilinear interpolation. Fills exposed edges with NaN."""
    h, w = img.shape
    if not (np.isfinite(dy) and np.isfinite(dx)):
        return np.full_like(img, np.nan)
    y, x = np.mgrid[0:h, 0:w].astype(np.float64)
    # sample coordinates in the source image
    sy = y - dy
    sx = x - dx
    y0 = np.floor(sy).astype(int)
    x0 = np.floor(sx).astype(int)
    fy = sy - y0
    fx = sx - x0

    valid = (y0 >= 0) & (y0 < h - 1) & (x0 >= 0) & (x0 < w - 1)
    y0c = np.clip(y0, 0, h - 2)
    x0c = np.clip(x0, 0, w - 2)

    Ia = img[y0c, x0c]
    Ib = img[y0c, x0c + 1]
    Ic = img[y0c + 1, x0c]
    Id = img[y0c + 1, x0c + 1]

    out = ((1 - fy) * (1 - fx) * Ia + (1 - fy) * fx * Ib +
           fy * (1 - fx) * Ic + fy * fx * Id)
    out[~valid] = np.nan
    return out


def align_stack(images: list[np.ndarray],
                ref_index: int = 0,
                max_shift_px: float = 2.0
                ) -> tuple[list[np.ndarray], list[tuple[float, float]]]:
    """Align all images to images[ref_index].

    Frames whose measured shift exceeds `max_shift_px` are left unaligned
    (shift recorded as measured) — on faint fields phase correlation can
    lock onto noise, and a wrong shift is worse than none.

    Returns (aligned_images, shifts) where shifts[i] = (dy, dx) applied to
    frame i (0,0 for the reference).
    """
    ref = images[ref_index]
    aligned: list[np.ndarray] = [None] * len(images)  # type: ignore[list-item]
    shifts: list[tuple[float, float]] = [(0.0, 0.0)] * len(images)
    aligned[ref_index] = images[ref_index]
    for i, img in enumerate(images):
        if i == ref_index:
            continue
        dy, dx = _phase_shift(ref, img)
        if (dy * dy + dx * dx) ** 0.5 > max_shift_px:
            aligned[i] = img  # reject implausible shift
        else:
            aligned[i] = shift_image(img, dy, dx)
        shifts[i] = (dy, dx)
    return aligned, shifts  # type: ignore[return-value]


def median_template(aligned: list[np.ndarray]) -> np.ndarray:
    stack = np.stack(aligned)
    with np.errstate(all="ignore"):
        return np.nanmedian(stack, axis=0)


def residual_frames(aligned: list[np.ndarray],
                    template: np.ndarray) -> list[np.ndarray]:
    return [a - template for a in aligned]


def robust_limits(residuals: list[np.ndarray],
                  nsigma: float = 5.0) -> tuple[float, float]:
    """Symmetric display limits from robust std of all residuals."""
    allv = np.concatenate([r[np.isfinite(r)] for r in residuals])
    if allv.size == 0:
        return -1.0, 1.0
    med = np.median(allv)
    mad = np.median(np.abs(allv - med)) * 1.4826
    lim = max(nsigma * mad, 1e-12)
    return float(med - lim), float(med + lim)


def render_residual_png(resid: np.ndarray, png_path: str | Path,
                        vmin: float, vmax: float, size: int = 512) -> Path:
    """Render a residual frame with a diverging colormap."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    norm = np.clip((resid - vmin) / (vmax - vmin), 0.0, 1.0)
    norm[~np.isfinite(norm)] = 0.5
    cmap_obj = plt.get_cmap("RdBu_r")
    rgba = cmap_obj(norm)
    rgb = (rgba[..., :3] * 255).astype(np.uint8)
    img = Image.fromarray(rgb).resize((size, size), Image.BILINEAR)
    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    img.save(png_path)
    return png_path


def build_difference_stack(fits_paths: list[str | Path],
                           out_dir: str | Path,
                           align: bool = True,
                           align_method: str = "wcs",
                           ref_index: int = 0) -> dict[str, Any]:
    """Full pipeline: load -> align -> median subtract -> render PNGs.

    Parameters
    ----------
    fits_paths : list of (str or Path)
        Paths to SPHEREx cutout FITS files.
    out_dir : str or Path
        Output directory for rendered diff PNGs.
    align : bool, default True
        Whether to align images before template construction.
    align_method : str, default "wcs"
        Alignment strategy:
        - "wcs" (Phase 5): Astropy WCS celestial grid reprojection.
        - "phase" (Phase 3): FFT phase-correlation translation.
    ref_index : int, default 0
        Reference frame index used for the celestial reference grid.

    Returns
    -------
    dict with keys:
        'n_frames', 'frames' (PNG paths), 'template_png', 'limits',
        'shifts', 'alignment_method', 'alignment_metadata', 'out_dir'.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_total = len(fits_paths)
    if n_total == 0:
        raise ValueError("fits_paths cannot be empty")

    shifts: list[tuple[float, float]] = [(0.0, 0.0)] * n_total
    alignment_meta: list[dict[str, Any]] = []
    method_used = "none"

    if align and align_method == "wcs":
        try:
            wcs_res = align_fits_stack_wcs(fits_paths, ref_index=ref_index)
            images = wcs_res["aligned_images"]
            alignment_meta = wcs_res["metadata"]
            valid_idx = [i for i, m in enumerate(alignment_meta) if m["valid"]]
            if not valid_idx:
                raise ValueError("No valid frames produced by WCS alignment")
            method_used = "wcs"
        except Exception as exc:
            logger.warning(
                "[diff_imaging] WCS alignment failed (%s); falling back to phase correlation",
                exc
            )
            align_method = "phase"

    if not align or align_method != "wcs":
        # Load images raw; keep index alignment by substituting NaN planes
        raw: list[np.ndarray] = []
        ref_shape: tuple[int, ...] | None = None
        for p in fits_paths:
            try:
                img = load_image(p)
                if img.ndim != 2 or img.size < 16:
                    raise ValueError(f"bad image shape {img.shape} in {p}")
                if ref_shape is None:
                    ref_shape = img.shape
                elif img.shape != ref_shape:
                    raise ValueError(f"shape mismatch {img.shape} vs {ref_shape} in {p}")
                raw.append(img)
            except Exception as exc:  # noqa: BLE001
                logger.warning("[diff_imaging] skipping bad frame %s: %s", p, exc)
                raw.append(np.full(ref_shape or (0, 0), np.nan))

        images = raw
        valid_idx = [i for i, im in enumerate(images)
                     if np.isfinite(im).sum() > im.size * 0.5]
        if not valid_idx:
            raise ValueError("no usable frames in stack")

        if align and align_method == "phase" and len(valid_idx) > 1:
            sub = [images[i] for i in valid_idx]
            aligned_sub, sub_shifts = align_stack(sub)
            for k, i in enumerate(valid_idx):
                images[i] = aligned_sub[k]
                shifts[i] = sub_shifts[k]
            method_used = "phase"
        else:
            method_used = "none"

    valid_aligned = [images[i] for i in valid_idx]
    valid_set = set(valid_idx)
    template = median_template(valid_aligned)
    th, tw = template.shape
    residuals = [
        (images[i] - template) if i in valid_set
        else np.full((th, tw), np.nan)
        for i in range(len(images))
    ]
    vmin, vmax = robust_limits([r for r in residuals if np.isfinite(r).any()])

    png_paths = []
    for i, r in enumerate(residuals):
        png_paths.append(render_residual_png(r, out_dir / f"diff_{i:04d}.png",
                                             vmin, vmax))

    # Render the template itself for reference
    t_norm = _asinh_stretch(template)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    rgba = plt.get_cmap("magma")(np.clip(t_norm, 0, 1))
    rgb = (rgba[..., :3] * 255).astype(np.uint8)
    Image.fromarray(rgb).resize((512, 512), Image.BILINEAR).save(
        out_dir / "template.png")

    return {
        "n_frames": len(png_paths),
        "frames": [str(p) for p in png_paths],
        "template_png": str(out_dir / "template.png"),
        "limits": [vmin, vmax],
        "shifts": shifts,
        "alignment_method": method_used,
        "alignment_metadata": alignment_meta,
        "out_dir": str(out_dir),
    }
