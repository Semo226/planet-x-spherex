"""Unit tests for Phase 5 WCS astronomical alignment and registration."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from processing.diff_imaging import build_difference_stack
from processing.wcs_alignment import (
    align_fits_stack_wcs,
    get_wcs_metadata,
    read_fits_wcs,
    reproject_frame_to_wcs,
)


def _make_dummy_fits(path: Path,
                     center_ra: float = 133.795,
                     center_dec: float = -7.245,
                     shape: tuple[int, int] = (29, 29),
                     offset_px: tuple[float, float] = (0.0, 0.0),
                     flux_val: float = 1.0) -> Path:
    """Create a synthetic SPHEREx MEF FITS file with a valid WCS IMAGE HDU."""
    ny, nx = shape
    data = np.ones((ny, nx), dtype=np.float32) * flux_val
    # Add a mock point source
    cy, cx = ny // 2, nx // 2
    data[cy, cx] += 10.0

    w = WCS(naxis=2)
    w.wcs.crpix = [nx / 2.0 + 0.5 + offset_px[1], ny / 2.0 + 0.5 + offset_px[0]]
    w.wcs.cdelt = [-0.0017, 0.0017]  # ~6.12 arcsec/px
    w.wcs.crval = [center_ra, center_dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    header = w.to_header()
    header["EXTNAME"] = "IMAGE"
    header["DETECTOR"] = 5
    header["OBSID"] = "test_obs"
    header["MJD-AVG"] = 60800.0

    primary_hdu = fits.PrimaryHDU()
    image_hdu = fits.ImageHDU(data=data, header=header, name="IMAGE")
    hdul = fits.HDUList([primary_hdu, image_hdu])
    hdul.writeto(path, overwrite=True)
    return path


class TestWCSAlignment(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.tmp_dir.name)

        # Create two slightly offset synthetic FITS cutouts
        self.f0 = _make_dummy_fits(
            self.dir_path / "frame_0.fits",
            offset_px=(0.0, 0.0),
            flux_val=2.0,
        )
        self.f1 = _make_dummy_fits(
            self.dir_path / "frame_1.fits",
            offset_px=(1.5, -1.0),
            flux_val=2.0,
        )

    def tearDown(self):
        self.tmp_dir.cleanup()

    def test_read_fits_wcs(self):
        data, wcs, hdr = read_fits_wcs(self.f0, hdu_name="IMAGE")
        self.assertEqual(data.shape, (29, 29))
        self.assertIsInstance(wcs, WCS)
        self.assertEqual(hdr["EXTNAME"], "IMAGE")

    def test_get_wcs_metadata(self):
        _, wcs, _ = read_fits_wcs(self.f0)
        meta = get_wcs_metadata(wcs, (29, 29))
        self.assertAlmostEqual(meta["center_ra_deg"], 133.795, places=3)
        self.assertAlmostEqual(meta["center_dec_deg"], -7.245, places=3)
        self.assertIsNotNone(meta["pixel_scale_arcsec"])

    def test_reproject_frame_to_wcs(self):
        d0, w0, _ = read_fits_wcs(self.f0)
        d1, w1, _ = read_fits_wcs(self.f1)
        reproj, fp = reproject_frame_to_wcs(d1, w1, w0, (29, 29))
        self.assertEqual(reproj.shape, (29, 29))
        self.assertEqual(fp.shape, (29, 29))
        self.assertTrue(np.sum(fp > 0) > 0.8 * fp.size)

    def test_align_fits_stack_wcs(self):
        res = align_fits_stack_wcs([self.f0, self.f1], ref_index=0)
        self.assertEqual(res["n_frames"], 2)
        self.assertEqual(res["n_aligned"], 2)
        self.assertEqual(len(res["aligned_images"]), 2)
        self.assertEqual(res["aligned_images"][0].shape, (29, 29))
        self.assertEqual(res["aligned_images"][1].shape, (29, 29))
        # Ref frame overlap is exactly 1.0
        self.assertEqual(res["metadata"][0]["overlap_fraction"], 1.0)
        # Second frame has large overlap
        self.assertGreater(res["metadata"][1]["overlap_fraction"], 0.8)

    def test_build_difference_stack_wcs(self):
        out_dir = self.dir_path / "diff_out"
        res = build_difference_stack(
            [self.f0, self.f1],
            out_dir=out_dir,
            align=True,
            align_method="wcs",
        )
        self.assertEqual(res["n_frames"], 2)
        self.assertEqual(res["alignment_method"], "wcs")
        self.assertTrue(Path(res["template_png"]).exists())
        for f in res["frames"]:
            self.assertTrue(Path(f).exists())


if __name__ == "__main__":
    unittest.main()
