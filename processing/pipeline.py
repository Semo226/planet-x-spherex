"""End-to-end Phase 5 candidate detection pipeline for SPHEREx.

Orchestrates WCS alignment, temporal frame differencing, false positive
filtering, and motion detection to identify candidate transient or moving
sources in multi-epoch SPHEREx cutout stacks.

Pipeline stages:
  1. WCS-align the input FITS stack onto a common celestial grid.
  2. Build a median "static sky" template and compute per-frame residuals.
  3. Detect candidate peaks above a configurable sigma threshold in the
     mean absolute residual map.
  4. For each candidate pixel, read FLAGS and VARIANCE HDUs from the
     reference FITS to run false-positive filtering.
  5. Evaluate detector flags, variance noise, and PSF morphology.
  6. Compute integrated confidence scores.
  7. Track per-epoch pixel-level motion across individual residual frames.
  8. Classify each candidate as high_priority, low_priority, or artifact.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
from astropy.wcs import WCS

from processing.candidate_filter import (
    calculate_candidate_confidence,
    evaluate_flags_at_pixel,
    evaluate_psf_consistency,
    evaluate_variance_noise,
)
from processing.diff_imaging import median_template, residual_frames
from processing.wcs_alignment import (
    align_fits_stack_wcs,
    get_wcs_metadata,
    read_fits_bundle,
)

logger = logging.getLogger(__name__)


def _detect_peaks(mean_abs_residual: np.ndarray,
                  detection_sigma: float) -> list[tuple[int, int, float]]:
    """Find pixel coordinates where the mean absolute residual exceeds a
    threshold defined as ``detection_sigma * robust_std``.

    Parameters
    ----------
    mean_abs_residual : np.ndarray
        2D mean absolute residual map.
    detection_sigma : float
        Number of robust standard deviations above the background to
        qualify as a candidate peak.

    Returns
    -------
    list of (y, x, snr)
        Peak coordinates and their approximate SNR above the noise floor.
    """
    finite = mean_abs_residual[np.isfinite(mean_abs_residual)]
    if finite.size == 0:
        return []

    median_bg = float(np.median(finite))
    mad = float(np.median(np.abs(finite - median_bg)))
    robust_std = mad * 1.4826  # MAD-to-sigma conversion

    if robust_std <= 0:
        logger.warning("[pipeline] Robust std is zero — cannot detect peaks")
        return []

    threshold = median_bg + detection_sigma * robust_std

    above_mask = (mean_abs_residual > threshold) & np.isfinite(mean_abs_residual)
    peak_ys, peak_xs = np.where(above_mask)

    peaks: list[tuple[int, int, float]] = []
    for py, px in zip(peak_ys, peak_xs):
        val = float(mean_abs_residual[py, px])
        snr = (val - median_bg) / robust_std
        peaks.append((int(py), int(px), float(snr)))

    # Sort by descending SNR
    peaks.sort(key=lambda t: t[2], reverse=True)
    return peaks


def _extract_stamp(image: np.ndarray,
                   y: int, x: int,
                   radius: int) -> np.ndarray:
    """Extract a square stamp of size ``(2*radius+1)`` around ``(y, x)``.

    Parameters
    ----------
    image : np.ndarray
        2D image array.
    y, x : int
        Center pixel coordinates.
    radius : int
        Half-size of the stamp (stamp width = 2*radius + 1).

    Returns
    -------
    np.ndarray
        Stamp array, zero-padded where the extraction exceeds image bounds.
    """
    h, w = image.shape
    size = 2 * radius + 1
    stamp = np.zeros((size, size), dtype=np.float64)

    y_lo = max(0, y - radius)
    y_hi = min(h, y + radius + 1)
    x_lo = max(0, x - radius)
    x_hi = min(w, x + radius + 1)

    sy_lo = y_lo - (y - radius)
    sy_hi = sy_lo + (y_hi - y_lo)
    sx_lo = x_lo - (x - radius)
    sx_hi = sx_lo + (x_hi - x_lo)

    stamp[sy_lo:sy_hi, sx_lo:sx_hi] = image[y_lo:y_hi, x_lo:x_hi]
    return stamp


def _count_epoch_detections(residuals: list[np.ndarray],
                            y: int, x: int,
                            residual_std: float,
                            detection_sigma: float) -> int:
    """Count how many individual residual frames show a detection at (y, x).

    Parameters
    ----------
    residuals : list of np.ndarray
        Per-epoch residual frames.
    y, x : int
        Pixel coordinates of the candidate.
    residual_std : float
        Global robust standard deviation of the residual stack.
    detection_sigma : float
        Significance threshold (same as the main detection pass).

    Returns
    -------
    int
        Number of epochs where abs(residual) exceeds the threshold.
    """
    threshold = detection_sigma * residual_std
    count = 0
    for resid in residuals:
        if y < resid.shape[0] and x < resid.shape[1]:
            val = resid[y, x]
            if np.isfinite(val) and abs(val) > threshold:
                count += 1
    return count


def _measure_motion(residuals: list[np.ndarray],
                    y: int, x: int,
                    search_radius: int,
                    pixel_scale_arcsec: float | None,
                    residual_std: float,
                    detection_sigma: float) -> dict[str, Any]:
    """Track per-epoch peak-residual shifts around a candidate position.

    For each residual frame, finds the peak absolute value within a small
    search box.  If the peak position shifts between epochs, the candidate
    is flagged as potentially moving.

    Parameters
    ----------
    residuals : list of np.ndarray
        Per-epoch residual frames.
    y, x : int
        Nominal candidate pixel coordinates.
    search_radius : int
        Half-size of the search box around (y, x).
    pixel_scale_arcsec : float or None
        Pixel scale in arcseconds (used to convert displacement).
    residual_std : float
        Robust standard deviation of the residual stack.
    detection_sigma : float
        Significance threshold for considering a frame "detected".

    Returns
    -------
    dict
        Motion analysis with ``is_moving``, ``displacement_arcsec``,
        ``direction_deg``, and ``epoch_positions``.
    """
    h, w = residuals[0].shape
    threshold = detection_sigma * residual_std
    positions: list[tuple[int, int]] = []

    for resid in residuals:
        y_lo = max(0, y - search_radius)
        y_hi = min(h, y + search_radius + 1)
        x_lo = max(0, x - search_radius)
        x_hi = min(w, x + search_radius + 1)

        sub = np.abs(resid[y_lo:y_hi, x_lo:x_hi])
        if not np.isfinite(sub).any():
            continue

        sub_finite = np.where(np.isfinite(sub), sub, 0.0)
        peak_val = float(np.max(sub_finite))
        if peak_val < threshold:
            continue

        local_iy, local_ix = np.unravel_index(np.argmax(sub_finite), sub.shape)
        abs_y = int(y_lo + local_iy)
        abs_x = int(x_lo + local_ix)
        positions.append((abs_y, abs_x))

    if len(positions) < 2:
        return {
            "is_moving": False,
            "displacement_arcsec": 0.0,
            "direction_deg": 0.0,
            "epoch_positions": positions,
        }

    # Compute displacement between first and last detected position
    dy = positions[-1][0] - positions[0][0]
    dx = positions[-1][1] - positions[0][1]
    displacement_px = float(np.hypot(dy, dx))

    scale = pixel_scale_arcsec if pixel_scale_arcsec is not None else 6.2
    displacement_arcsec = displacement_px * scale

    direction_deg = float(np.degrees(np.arctan2(dx, -dy))) % 360.0

    # Require at least 1 pixel of coherent displacement to flag as moving
    is_moving = displacement_px >= 1.0

    return {
        "is_moving": is_moving,
        "displacement_arcsec": round(displacement_arcsec, 3),
        "direction_deg": round(direction_deg, 1),
        "epoch_positions": positions,
    }


def _classify_candidate(confidence: float,
                         psf_eval: dict[str, Any],
                         flag_eval: dict[str, Any],
                         motion: dict[str, Any]) -> str:
    """Assign a classification label based on filter results.

    Parameters
    ----------
    confidence : float
        Integrated confidence score (0-100).
    psf_eval : dict
        Result from ``evaluate_psf_consistency``.
    flag_eval : dict
        Result from ``evaluate_flags_at_pixel``.
    motion : dict
        Result from ``_measure_motion``.

    Returns
    -------
    str
        One of ``"high_priority"``, ``"low_priority"``, or ``"artifact"``.
    """
    morphology = psf_eval.get("morphology", "indeterminate")
    flag_penalty = flag_eval.get("penalty", 0.0)

    # Artifacts: cosmic rays, flagged regions, or very low confidence
    if morphology == "cosmic_ray":
        return "artifact"
    if flag_penalty >= 0.8:
        return "artifact"
    if confidence < 20.0:
        return "artifact"

    # High priority: good PSF, clean flags, decent confidence
    if confidence >= 50.0 and morphology == "point_source" and flag_penalty < 0.3:
        return "high_priority"

    # Boost moving candidates — even moderate confidence can be high priority
    if motion.get("is_moving", False) and confidence >= 35.0 and flag_penalty < 0.5:
        return "high_priority"

    return "low_priority"


def run_candidate_pipeline(fits_paths: list[str | Path],
                           ra: float, dec: float,
                           detection_sigma: float = 3.0,
                           stamp_radius: int = 2) -> dict:
    """Run the end-to-end Phase 5 candidate detection pipeline.

    Parameters
    ----------
    fits_paths : list of str or Path
        Chronological list of SPHEREx multi-extension FITS cutout paths.
    ra : float
        Right ascension of the query center (degrees, ICRS).
    dec : float
        Declination of the query center (degrees, ICRS).
    detection_sigma : float, default 3.0
        Number of robust standard deviations above background required for
        a pixel to be flagged as a candidate peak.
    stamp_radius : int, default 2
        Half-width (in pixels) of the postage stamp extracted around each
        candidate for PSF morphology analysis.  A radius of 2 yields a
        5×5 stamp.

    Returns
    -------
    dict
        Pipeline results with keys:

        - ``"candidates"``: list of candidate dicts (see below).
        - ``"n_candidates"``: int — number of candidates found.
        - ``"alignment"``: dict with ``method``, ``n_aligned``, ``ref_index``.
        - ``"stack_stats"``: dict with ``residual_std``, ``median_background``,
          ``n_valid_frames``.

        Each candidate dict contains:
            ``ra_deg``, ``dec_deg``, ``pixel_y``, ``pixel_x``, ``snr``,
            ``n_epoch_detections``, ``confidence``, ``classification``,
            ``psf_eval``, ``flag_eval``, ``variance_eval``, ``motion``.

    Raises
    ------
    ValueError
        If ``fits_paths`` is empty or no frames can be aligned.
    """
    if not fits_paths:
        raise ValueError("fits_paths list cannot be empty")

    logger.info(
        "[pipeline] Starting candidate pipeline: %d frames, "
        "center=(%.5f, %.5f), sigma=%.1f",
        len(fits_paths), ra, dec, detection_sigma,
    )

    # ------------------------------------------------------------------
    # Stage 1: WCS alignment
    # ------------------------------------------------------------------
    logger.info("[pipeline] Stage 1: WCS alignment")
    alignment_result = align_fits_stack_wcs(fits_paths)

    aligned_images: list[np.ndarray] = alignment_result["aligned_images"]
    ref_wcs: WCS = alignment_result["ref_wcs"]
    ref_index: int = alignment_result["ref_index"]
    shape_out: tuple[int, int] = alignment_result["shape_out"]
    n_aligned: int = alignment_result["n_aligned"]

    if n_aligned == 0:
        raise ValueError("No frames were successfully WCS-aligned")

    # Identify which frames are valid (not all-NaN)
    valid_mask = [
        m["valid"] for m in alignment_result["metadata"]
    ]
    valid_images = [
        img for img, v in zip(aligned_images, valid_mask) if v
    ]

    wcs_meta = get_wcs_metadata(ref_wcs, shape_out)
    pixel_scale_arcsec: float | None = wcs_meta.get("pixel_scale_arcsec")

    logger.info(
        "[pipeline] Aligned %d / %d frames (ref=%d, scale=%.2f\"/px)",
        n_aligned, len(fits_paths), ref_index,
        pixel_scale_arcsec or 0.0,
    )

    # ------------------------------------------------------------------
    # Stage 2: Median template and residuals
    # ------------------------------------------------------------------
    logger.info("[pipeline] Stage 2: median template + residual frames")
    template = median_template(valid_images)
    residuals = residual_frames(aligned_images, template)

    # Mean absolute residual map (over all frames, ignoring NaN)
    abs_stack = np.stack([np.abs(r) for r in residuals])
    with np.errstate(all="ignore"):
        mean_abs_resid = np.nanmean(abs_stack, axis=0)

    # Global stack statistics
    finite_residuals = np.concatenate(
        [r[np.isfinite(r)].ravel() for r in residuals]
    )
    if finite_residuals.size > 0:
        residual_median = float(np.median(finite_residuals))
        residual_mad = float(np.median(np.abs(finite_residuals - residual_median)))
        residual_std = residual_mad * 1.4826
    else:
        residual_median = 0.0
        residual_std = 1.0

    stack_stats = {
        "residual_std": round(residual_std, 6),
        "median_background": round(residual_median, 6),
        "n_valid_frames": len(valid_images),
    }

    logger.info(
        "[pipeline] Residual stats: std=%.4e, median_bg=%.4e, n_valid=%d",
        residual_std, residual_median, len(valid_images),
    )

    # ------------------------------------------------------------------
    # Stage 3: Detect candidate peaks
    # ------------------------------------------------------------------
    logger.info("[pipeline] Stage 3: peak detection (sigma=%.1f)", detection_sigma)
    peaks = _detect_peaks(mean_abs_resid, detection_sigma)
    logger.info("[pipeline] Found %d raw peaks above %.1f-sigma", len(peaks), detection_sigma)

    # ------------------------------------------------------------------
    # Stage 4-7: Evaluate each candidate
    # ------------------------------------------------------------------
    logger.info("[pipeline] Stage 4-7: candidate evaluation + motion tracking")

    # Load the reference FITS bundle for FLAGS / VARIANCE data
    ref_path = fits_paths[ref_index]
    try:
        ref_bundle = read_fits_bundle(ref_path)
        flags_map: np.ndarray | None = ref_bundle.get("flags")
        variance_map: np.ndarray | None = ref_bundle.get("variance")
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "[pipeline] Could not read reference bundle %s: %s", ref_path, exc,
        )
        flags_map = None
        variance_map = None

    candidates: list[dict[str, Any]] = []

    for peak_y, peak_x, peak_snr in peaks:
        try:
            # --- Sky coordinate conversion ---
            sky_coord = ref_wcs.pixel_to_world(float(peak_x), float(peak_y))
            ra_deg = float(sky_coord.ra.deg)
            dec_deg = float(sky_coord.dec.deg)

            # --- Stage 4–5: Flag & variance evaluation ---
            signal_flux = float(mean_abs_resid[peak_y, peak_x])

            flag_eval = evaluate_flags_at_pixel(
                flags_map, y=peak_y, x=peak_x, radius=stamp_radius,
            )
            variance_eval = evaluate_variance_noise(
                variance_map, y=peak_y, x=peak_x, signal_flux=signal_flux,
            )

            # --- PSF consistency ---
            stamp = _extract_stamp(mean_abs_resid, peak_y, peak_x, stamp_radius)
            psf_eval = evaluate_psf_consistency(stamp)

            # --- Epoch detection count ---
            n_epoch_detections = _count_epoch_detections(
                residuals, peak_y, peak_x, residual_std, detection_sigma,
            )

            # --- Stage 7: Motion tracking ---
            motion = _measure_motion(
                residuals, peak_y, peak_x,
                search_radius=stamp_radius + 1,
                pixel_scale_arcsec=pixel_scale_arcsec,
                residual_std=residual_std,
                detection_sigma=detection_sigma,
            )

            # --- Stage 6: Confidence score ---
            confidence = calculate_candidate_confidence(
                snr=peak_snr,
                psf_score=psf_eval["psf_score"],
                flag_penalty=flag_eval.get("penalty", 0.0),
                n_epoch_detections=n_epoch_detections,
                is_moving=motion["is_moving"],
            )

            # --- Classification ---
            classification = _classify_candidate(
                confidence, psf_eval, flag_eval, motion,
            )

            # Strip epoch_positions from the returned motion dict (internal)
            motion_out = {
                "is_moving": motion["is_moving"],
                "displacement_arcsec": motion["displacement_arcsec"],
                "direction_deg": motion["direction_deg"],
            }

            candidates.append({
                "ra_deg": round(ra_deg, 7),
                "dec_deg": round(dec_deg, 7),
                "pixel_y": peak_y,
                "pixel_x": peak_x,
                "snr": round(peak_snr, 2),
                "n_epoch_detections": n_epoch_detections,
                "confidence": confidence,
                "classification": classification,
                "psf_eval": psf_eval,
                "flag_eval": flag_eval,
                "variance_eval": variance_eval,
                "motion": motion_out,
            })

        except Exception:  # noqa: BLE001
            logger.exception(
                "[pipeline] Error evaluating candidate at pixel (%d, %d)",
                peak_y, peak_x,
            )
            continue

    # Sort final candidates: high_priority first, then by descending SNR
    classification_order = {"high_priority": 0, "low_priority": 1, "artifact": 2}
    candidates.sort(
        key=lambda c: (classification_order.get(c["classification"], 9), -c["snr"]),
    )

    logger.info(
        "[pipeline] Pipeline complete: %d candidates "
        "(high=%d, low=%d, artifact=%d)",
        len(candidates),
        sum(1 for c in candidates if c["classification"] == "high_priority"),
        sum(1 for c in candidates if c["classification"] == "low_priority"),
        sum(1 for c in candidates if c["classification"] == "artifact"),
    )

    return {
        "candidates": candidates,
        "n_candidates": len(candidates),
        "alignment": {
            "method": "wcs",
            "n_aligned": n_aligned,
            "ref_index": ref_index,
        },
        "stack_stats": stack_stats,
    }
