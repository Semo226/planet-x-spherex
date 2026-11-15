"""Unit tests for Phase 5 candidate pipeline and API integration."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from backend.main import app, _session, _lock
from processing.candidate_filter import (
    calculate_candidate_confidence,
    evaluate_flags_at_pixel,
    evaluate_psf_consistency,
    evaluate_variance_noise,
)
from processing.pipeline import run_candidate_pipeline


def _make_mock_fits(path: Path,
                    offset_px: tuple[float, float] = (0.0, 0.0),
                    point_source_flux: float = 15.0) -> Path:
    """Create a mock SPHEREx cutout MEF with IMAGE, FLAGS, and VARIANCE."""
    ny, nx = 29, 29
    data = np.ones((ny, nx), dtype=np.float32) * 2.0
    cy, cx = ny // 2, nx // 2
    # Mock moving point source
    sy = int(round(cy + offset_px[0]))
    sx = int(round(cx + offset_px[1]))
    if 0 <= sy < ny and 0 <= sx < nx:
        data[sy, sx] += point_source_flux

    w = WCS(naxis=2)
    w.wcs.crpix = [nx / 2.0 + 0.5, ny / 2.0 + 0.5]
    w.wcs.cdelt = [-0.0017, 0.0017]
    w.wcs.crval = [133.795, -7.245]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    hdr = w.to_header()
    hdr["EXTNAME"] = "IMAGE"
    hdr["DETECTOR"] = 5
    hdr["OBSID"] = "mock_obs"
    hdr["MJD-AVG"] = 60800.0

    primary = fits.PrimaryHDU()
    im_hdu = fits.ImageHDU(data=data, header=hdr, name="IMAGE")
    flg_hdu = fits.ImageHDU(data=np.zeros((ny, nx), dtype=np.int32), name="FLAGS")
    var_hdu = fits.ImageHDU(data=np.full((ny, nx), 0.0004, dtype=np.float32), name="VARIANCE")

    hdul = fits.HDUList([primary, im_hdu, flg_hdu, var_hdu])
    hdul.writeto(path, overwrite=True)
    return path


class TestCandidatePipeline(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.p = Path(self.tmp.name)
        self.f0 = _make_mock_fits(self.p / "f0.fits", offset_px=(0.0, 0.0))
        self.f1 = _make_mock_fits(self.p / "f1.fits", offset_px=(1.0, 1.0))
        self.f2 = _make_mock_fits(self.p / "f2.fits", offset_px=(2.0, 2.0))

    def tearDown(self):
        self.tmp.cleanup()

    def test_filter_functions(self):
        # Clean flags
        clean_flags = np.zeros((10, 10), dtype=np.int32)
        res = evaluate_flags_at_pixel(clean_flags, 5, 5)
        self.assertTrue(res["clean"])
        self.assertEqual(res["penalty"], 0.0)

        # Dirty flags (cosmic ray)
        dirty_flags = np.zeros((10, 10), dtype=np.int32)
        dirty_flags[5, 5] = 64
        res2 = evaluate_flags_at_pixel(dirty_flags, 5, 5)
        self.assertFalse(res2["clean"])
        self.assertGreater(res2["penalty"], 0.5)

        # Variance
        var_map = np.full((10, 10), 0.01)
        res_var = evaluate_variance_noise(var_map, 5, 5, signal_flux=1.0)
        self.assertTrue(res_var["noise_consistent"])

        # Confidence
        conf = calculate_candidate_confidence(snr=5.0, psf_score=80.0, flag_penalty=0.0, n_epoch_detections=3)
        self.assertGreater(conf, 60.0)

    def test_run_candidate_pipeline(self):
        paths = [self.f0, self.f1, self.f2]
        res = run_candidate_pipeline(paths, ra=133.795, dec=-7.245, detection_sigma=1.5)
        self.assertIn("candidates", res)
        self.assertIn("alignment", res)
        self.assertIn("stack_stats", res)
        self.assertGreater(res["alignment"]["n_aligned"], 0)

    def test_api_candidates_endpoint(self):
        from backend.main import candidates
        # Populate session
        with _lock:
            _session["frames"] = [
                {"index": 0, "fits_path": str(self.f0), "detector": "D5", "obs_id": "0", "bandpass": "D5", "t_min": 60800.0, "t_min_iso": "2025-05-01T00:00:00Z"},
                {"index": 1, "fits_path": str(self.f1), "detector": "D5", "obs_id": "1", "bandpass": "D5", "t_min": 60801.0, "t_min_iso": "2025-05-02T00:00:00Z"},
            ]
            _session["request"] = {"ra": 133.795, "dec": -7.245, "fov_size": 0.05}

        data = candidates()
        self.assertIn("candidates", data)
        self.assertIn("n_candidates", data)


if __name__ == "__main__":
    unittest.main()
