"""False Positive Filtering, Noise Management & PSF Verification.

Phase 5 module: Evaluates candidate detections in SPHEREx difference images
against detector quality flags (FLAGS), noise variance maps (VARIANCE),
effective point spread function (EPSF) models, and spectral consistency.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

# SPHEREx Level 2 Flag Bit Definitions
FLAG_DO_NOT_USE = 1 << 0       # 1: Generic bad/dead pixel
FLAG_NONLINEAR = 1 << 1        # 2: Nonlinear detector response
FLAG_SATURATED = 1 << 2        # 4: Saturated pixel
FLAG_COSMIC_RAY = 1 << 6       # 64: Cosmic ray hit / jump detection
FLAG_PERSISTENCE = 1 << 17     # 131072: Latent persistence ghost
FLAG_EDGE_OR_LIMIT = 1 << 21   # 2097152: Array edge / FOV boundary
FLAG_NO_COVERAGE = 1 << 28     # 268435456: Uncovered / masked

CRITICAL_FLAGS_MASK = (
    FLAG_DO_NOT_USE
    | FLAG_SATURATED
    | FLAG_COSMIC_RAY
    | FLAG_NO_COVERAGE
)


def evaluate_flags_at_pixel(flags_map: np.ndarray | None,
                            y: float,
                            x: float,
                            radius: int = 1) -> dict[str, Any]:
    """Check detector quality flags in a small neighborhood around (y, x).

    Parameters
    ----------
    flags_map : np.ndarray or None
        2D integer bitmask from the FLAGS HDU.
    y, x : float
        Pixel coordinates of the candidate.
    radius : int, default 1
        Aperture radius in pixels to inspect (radius 1 inspects a 3x3 region).

    Returns
    -------
    dict
        Flag cleanliness evaluation: 'clean' (bool), 'flag_value' (int),
        'issues' (list of str), 'penalty' (float between 0.0 and 1.0).
    """
    if flags_map is None:
        return {"clean": True, "flag_value": 0, "issues": [], "penalty": 0.0}

    h, w = flags_map.shape
    iy, ix = int(round(y)), int(round(x))
    if iy < 0 or iy >= h or ix < 0 or ix >= w:
        return {
            "clean": False,
            "flag_value": int(FLAG_NO_COVERAGE),
            "issues": ["Candidate coordinate outside image boundaries"],
            "penalty": 1.0,
        }

    y_min = max(0, iy - radius)
    y_max = min(h, iy + radius + 1)
    x_min = max(0, ix - radius)
    x_max = min(w, ix + radius + 1)

    sub_flags = flags_map[y_min:y_max, x_min:x_max]
    center_flag = int(flags_map[iy, ix])
    combined_flags = int(np.bitwise_or.reduce(sub_flags.ravel()))

    issues: list[str] = []
    penalty = 0.0

    if combined_flags & FLAG_COSMIC_RAY:
        issues.append("Cosmic ray detected in aperture")
        penalty += 0.8
    if combined_flags & FLAG_SATURATED:
        issues.append("Saturated detector pixels")
        penalty += 0.7
    if combined_flags & FLAG_DO_NOT_USE:
        issues.append("Dead/bad pixel in aperture")
        penalty += 0.6
    if combined_flags & FLAG_NO_COVERAGE:
        issues.append("Edge/masked region without full coverage")
        penalty += 0.9
    if combined_flags & FLAG_PERSISTENCE:
        issues.append("Possible latent persistence artifact")
        penalty += 0.4
    if combined_flags & FLAG_NONLINEAR:
        issues.append("Non-linear response regime")
        penalty += 0.3

    penalty = min(1.0, penalty)
    clean = (penalty < 0.5)

    return {
        "clean": clean,
        "center_flag": center_flag,
        "combined_flags": combined_flags,
        "issues": issues,
        "penalty": penalty,
    }


def evaluate_variance_noise(variance_map: np.ndarray | None,
                            y: float,
                            x: float,
                            signal_flux: float) -> dict[str, Any]:
    """Evaluate local pixel variance against candidate signal.

    Parameters
    ----------
    variance_map : np.ndarray or None
        2D float array of per-pixel variance (MJy/sr)^2.
    y, x : float
        Candidate pixel coordinates.
    signal_flux : float
        Candidate peak or aperture flux.

    Returns
    -------
    dict
        Variance evaluation: 'local_sigma', 'pixel_snr', 'noise_consistent'.
    """
    if variance_map is None:
        return {
            "local_sigma": 1.0,
            "pixel_snr": None,
            "noise_consistent": True,
        }

    h, w = variance_map.shape
    iy, ix = int(round(y)), int(round(x))
    if iy < 0 or iy >= h or ix < 0 or ix >= w:
        return {
            "local_sigma": 1.0,
            "pixel_snr": None,
            "noise_consistent": False,
        }

    # Extract 5x5 local neighborhood around candidate
    y_min = max(0, iy - 2)
    y_max = min(h, iy + 3)
    x_min = max(0, ix - 2)
    x_max = min(w, ix + 3)

    sub_var = variance_map[y_min:y_max, x_min:x_max]
    finite_var = sub_var[np.isfinite(sub_var)]

    if finite_var.size == 0 or np.nanmedian(finite_var) <= 0:
        local_var = 1e-4
    else:
        local_var = float(np.nanmedian(finite_var))

    local_sigma = float(np.sqrt(local_var))
    pixel_snr = float(abs(signal_flux) / max(local_sigma, 1e-9))
    noise_consistent = (pixel_snr >= 2.0)

    return {
        "local_sigma": local_sigma,
        "pixel_snr": pixel_snr,
        "noise_consistent": noise_consistent,
    }


def evaluate_psf_consistency(cutout_subimage: np.ndarray,
                             expected_fwhm_px: float = 2.0) -> dict[str, Any]:
    """Test candidate spatial morphology against SPHEREx point spread function.

    Distinguishes genuine point sources (FWHM ~ 1.5-2.5 px) from single-pixel
    cosmic rays (compactness ~ 1.0, FWHM < 0.8 px) and broad scattered-light
    gradients (compactness < 0.2, FWHM > 4 px).

    Parameters
    ----------
    cutout_subimage : np.ndarray
        Small square stamp (e.g. 5x5 or 7x7) centered on candidate.
    expected_fwhm_px : float, default 2.0
        Expected SPHEREx beam FWHM in pixels.

    Returns
    -------
    dict with keys:
        - "psf_score": float (0.0 to 100.0)
        - "compactness": float (ratio of peak to 3x3 sum)
        - "correlation": float (Pearson correlation with 2D Gaussian PSF)
        - "is_point_source": bool
        - "morphology": "point_source" | "cosmic_ray" | "diffuse_artifact" | "indeterminate"
    """
    if cutout_subimage.ndim != 2 or cutout_subimage.size < 9:
        return {
            "psf_score": 50.0,
            "compactness": 0.4,
            "correlation": 0.5,
            "is_point_source": True,
            "morphology": "indeterminate",
        }

    # Normalize sub-raster so background is zero and peak is positive
    stamp = np.nan_to_num(cutout_subimage.copy(), nan=0.0)
    h, w = stamp.shape
    cy, cx = h // 2, w // 2

    peak_val = stamp[cy, cx]
    bg_median = float(np.median(stamp))
    peak_above_bg = peak_val - bg_median

    if peak_above_bg <= 0:
        # Negative peak (absorption / negative residual)
        stamp_clean = np.abs(stamp - bg_median)
        peak_above_bg = abs(peak_above_bg)
    else:
        stamp_clean = stamp - bg_median

    # 1. Compactness ratio: Peak flux / 3x3 flux
    y_min = max(0, cy - 1)
    y_max = min(h, cy + 2)
    x_min = max(0, cx - 1)
    x_max = min(w, cx + 2)
    flux_3x3 = float(np.sum(np.clip(stamp_clean[y_min:y_max, x_min:x_max], 0, None)))
    compactness = float(peak_above_bg / max(flux_3x3, 1e-9))

    # 2. Pearson correlation with 2D Gaussian model (sigma = FWHM / 2.355)
    sigma = expected_fwhm_px / 2.355
    y, x = np.mgrid[0:h, 0:w]
    model_psf = np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2.0 * sigma ** 2))

    a = stamp_clean.ravel()
    b = model_psf.ravel()
    a_norm = a - np.mean(a)
    b_norm = b - np.mean(b)
    denom = np.sqrt(np.sum(a_norm ** 2) * np.sum(b_norm ** 2))

    if denom > 1e-12:
        correlation = float(np.sum(a_norm * b_norm) / denom)
    else:
        correlation = 0.0

    # 3. Classify morphology
    if compactness > 0.85 and correlation < 0.65:
        morphology = "cosmic_ray"
        psf_score = max(5.0, 30.0 * max(0.0, correlation))
        is_point_source = False
    elif compactness < 0.15:
        morphology = "diffuse_artifact"
        psf_score = max(10.0, 40.0 * max(0.0, correlation))
        is_point_source = False
    else:
        morphology = "point_source"
        # Optimal point source has compactness ~ 0.35-0.65 and correlation > 0.7
        compact_factor = 1.0 - min(1.0, abs(compactness - 0.45) / 0.45)
        psf_score = float(np.clip(100.0 * max(0.0, correlation) * (0.5 + 0.5 * compact_factor), 10.0, 99.0))
        is_point_source = (psf_score >= 50.0)

    return {
        "psf_score": round(psf_score, 1),
        "compactness": round(compactness, 3),
        "correlation": round(correlation, 3),
        "is_point_source": is_point_source,
        "morphology": morphology,
    }


def calculate_candidate_confidence(snr: float,
                                   psf_score: float,
                                   flag_penalty: float,
                                   n_epoch_detections: int,
                                   is_moving: bool = False,
                                   is_known: bool = False) -> float:
    """Calculate an integrated, human-interpretable confidence score (0-100%).

    Parameters
    ----------
    snr : float
        Statistical signal-to-noise ratio in difference frames.
    psf_score : float
        PSF shape consistency score (0-100).
    flag_penalty : float
        Penalty from detector quality flags (0.0 to 1.0).
    n_epoch_detections : int
        Number of distinct epochs where candidate is detected.
    is_moving : bool, default False
        Whether target exhibits consistent coordinate motion across epochs.
    is_known : bool, default False
        Whether target matches a known catalog source.

    Returns
    -------
    float
        Confidence rating from 0.0 to 100.0.
    """
    # 1. Base score from SNR (0 to 45 pts)
    # SNR=3 -> 25 pts, SNR=5 -> 38 pts, SNR>=8 -> 45 pts
    snr_pts = float(np.clip(10.0 + 5.0 * snr, 0.0, 45.0))

    # 2. PSF score contribution (0 to 30 pts)
    psf_pts = float(np.clip(psf_score * 0.30, 0.0, 30.0))

    # 3. Multi-epoch persistence (0 to 20 pts)
    # 1 epoch -> 5 pts, 2 epochs -> 12 pts, >=3 epochs -> 20 pts
    epoch_pts = float(np.clip(n_epoch_detections * 7.0, 5.0, 20.0))

    # 4. Motion / astrophysical consistency bonus (0 to 10 pts)
    bonus = 8.0 if is_moving else (5.0 if is_known else 2.0)

    raw_score = (snr_pts + psf_pts + epoch_pts + bonus)

    # 5. Apply flag penalty
    clean_score = raw_score * (1.0 - 0.7 * flag_penalty)

    return float(np.clip(round(clean_score, 1), 5.0, 99.0))
