"""Render SPHEREx FITS cutouts to PNG for the web viewer.

Reads the IMAGE extension, applies an asinh/zscale-style stretch,
and writes a PNG next to the FITS file.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits

try:
    from PIL import Image
except ImportError as exc:  # pragma: no cover
    raise ImportError("Pillow is required: pip install pillow") from exc


def _asinh_stretch(data: np.ndarray, lower_pct: float = 1.0,
                   upper_pct: float = 99.5) -> np.ndarray:
    """Simple robust stretch: clip to percentiles then asinh."""
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return np.zeros_like(data, dtype=np.float32)
    lo, hi = np.percentile(finite, [lower_pct, upper_pct])
    clipped = np.clip(data, lo, hi)
    span = max(hi - lo, 1e-12)
    norm = (clipped - lo) / span
    # asinh soft-knee, scaled to [0, 1]
    return np.arcsinh(norm * 8.0) / np.arcsinh(8.0)


def render_fits_png(fits_path: str | Path, png_path: str | Path,
                    size: int = 512, cmap: str = "magma") -> Path:
    """Render the IMAGE HDU of a (cutout) MEF to a PNG. Returns png_path."""
    fits_path = Path(fits_path)
    png_path = Path(png_path)

    with fits.open(fits_path) as hdul:
        image_hdu = None
        for hdu in hdul:
            if (hdu.name or "").upper() == "IMAGE" and hdu.data is not None:
                image_hdu = hdu
                break
        if image_hdu is None:
            raise ValueError(f"no IMAGE extension in {fits_path}")
        data = np.asarray(image_hdu.data, dtype=np.float64)

    stretched = _asinh_stretch(data)

    # Apply colormap via matplotlib (already a project dependency)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cmap_obj = plt.get_cmap(cmap)
    rgba = cmap_obj(np.clip(stretched, 0.0, 1.0))
    rgb = (rgba[..., :3] * 255).astype(np.uint8)
    img = Image.fromarray(rgb).resize((size, size), Image.BILINEAR)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    img.save(png_path)
    return png_path


if __name__ == "__main__":
    import sys
    src = Path(sys.argv[1])
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else src.with_suffix(".png")
    print(render_fits_png(src, dst))
