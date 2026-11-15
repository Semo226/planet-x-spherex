"""FastAPI backend for Planet X and SPHEREx.

Endpoints:
    GET  /api/health
    POST /api/query   {ra, dec, fov_size?, max_frames?}
    GET  /api/frame/{index}.png   (after a query populates the session cache)
    GET  /                    -> frontend static UI

Run from repo root:
    .venv\\Scripts\\python.exe -m uvicorn backend.main:app --port 8000
"""

from __future__ import annotations

import sys
import threading
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from backend.data_service import get_spherex_cutouts  # noqa: E402
from backend.catalog_service import crossmatch  # noqa: E402
from processing.diff_imaging import build_difference_stack  # noqa: E402
from processing.pipeline import run_candidate_pipeline  # noqa: E402

app = FastAPI(title="Planet X and SPHEREx API", version="0.3.0")

FRONTEND_DIR = ROOT / "frontend"
CACHE_DIR = Path(__file__).resolve().parent / "cache"

# In-memory "session": latest query result. Fine for single-user hackathon.
_lock = threading.Lock()
_session: dict = {"frames": [], "diff": None, "candidates": None}


class QueryRequest(BaseModel):
    ra: float = Field(..., ge=0, le=360)
    dec: float = Field(..., ge=-90, le=90)
    fov_size: float = Field(0.05, gt=0, le=1.0)
    max_frames: int | None = Field(60, ge=1, le=300)
    detector: str | None = Field(None, description="Optional SPHEREx detector filter (e.g. D1-D6)")
    bandpass: str | None = Field(None, description="Optional SPHEREx bandpass name")


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.post("/api/query")
def query(req: QueryRequest):
    """Query IRSA, download+render cutouts, return frame metadata."""
    global _session
    result = get_spherex_cutouts(
        req.ra, req.dec,
        fov_size=req.fov_size,
        max_frames=req.max_frames,
    )
    if result["error"] and not result["frames"]:
        raise HTTPException(status_code=502, detail=result["error"])
    # Strip filesystem paths from the API payload except the png filename.
    frames = []
    for f in result["frames"]:
        frames.append({
            "index": f["index"],
            "obs_id": f["obs_id"],
            "bandpass": f["bandpass"],
            "detector": f["detector"],
            "t_min": f["t_min"],
            "t_max": f["t_max"],
            "t_min_iso": f["t_min_iso"],
            "png": f"/api/frame/{f['index']}.png",
        })
    with _lock:
        _session = {"frames": result["frames"], "diff": None, "candidates": None,
                    "request": {"ra": req.ra, "dec": req.dec,
                                "fov_size": req.fov_size,
                                "detector": req.detector,
                                "bandpass": req.bandpass}}
    return JSONResponse({
        "ra": result["ra"], "dec": result["dec"],
        "fov_size": result["fov_size"],
        "n_results": result["n_results"],
        "n_epochs": result["n_epochs"],
        "detector": req.detector,
        "offline": result.get("offline", False),
        "frames": frames,
    })


@app.get("/api/frame/{index}.png")
def frame(index: int):
    with _lock:
        frames = _session["frames"]
    if index < 0 or index >= len(frames):
        raise HTTPException(status_code=404, detail="frame out of range")
    png = Path(frames[index]["png_path"])
    if not png.exists():
        raise HTTPException(status_code=404, detail="png missing")
    return FileResponse(png, media_type="image/png")


class ProcessRequest(BaseModel):
    align: bool = True
    align_method: str = "wcs"   # 'wcs' (Phase 5 celestial WCS reprojection) | 'phase' (Phase 3 FFT)
    method: str = "median"      # 'median' (only method implemented so far)
    detector: str | None = None # e.g. 'D3'; None = auto (most common)


@app.post("/api/process")
def process(req: ProcessRequest):
    """Build aligned difference-imaging stack for the current session.

    Frames are grouped by detector array — SPHEREx's six arrays have
    independent distortion solutions, so mixing them produces misaligned
    residuals. `detector=None` auto-selects the most common detector.
    """
    with _lock:
        frames = list(_session["frames"])
    if not frames:
        raise HTTPException(status_code=400,
                            detail="no frames loaded — run /api/query first")

    if req.detector is None:
        from collections import Counter
        detector = Counter(f["detector"] for f in frames).most_common(1)[0][0]
    else:
        detector = req.detector
    sel = [f for f in frames if f["detector"] == detector]
    if not sel:
        raise HTTPException(status_code=400,
                            detail=f"no frames for detector {detector}")

    fits_paths = [f["fits_path"] for f in sel]
    import hashlib
    tag = hashlib.sha256(
        (detector + req.align_method + fits_paths[0] + fits_paths[-1]).encode()).hexdigest()[:12]
    out_dir = CACHE_DIR / "diff" / tag
    result = build_difference_stack(
        fits_paths,
        out_dir,
        align=req.align,
        align_method=req.align_method,
    )

    diff_frames = []
    for i, f in enumerate(sel):
        diff_frames.append({
            "index": i,
            "global_index": f["index"],
            "obs_id": f["obs_id"],
            "bandpass": f["bandpass"],
            "detector": f["detector"],
            "t_min": f["t_min"],
            "t_min_iso": f["t_min_iso"],
            "png": f"/api/diff/{i}.png",
        })
    with _lock:
        _session["diff"] = {"dir": result["out_dir"],
                            "n_frames": result["n_frames"],
                            "limits": result["limits"],
                            "template_png": result["template_png"],
                            "detector": detector,
                            "align_method": result.get("alignment_method", req.align_method),
                            "frames": diff_frames}
    return {
        "detector": detector,
        "align_method": result.get("alignment_method", req.align_method),
        "n_frames": result["n_frames"],
        "limits": result["limits"],
        "template_png": "/api/diff/template.png",
        "frames": diff_frames,
    }


@app.post("/api/candidates")
def candidates():
    """Run the Phase 5 candidate detection pipeline on the current session frames.

    Uses WCS-aligned difference imaging, false positive filtering (FLAGS/VARIANCE/PSF),
    and motion detection to produce a scored candidate table.
    """
    with _lock:
        frames = list(_session["frames"])
        req_data = _session.get("request", {})
    if not frames:
        raise HTTPException(status_code=400,
                            detail="no frames loaded — run /api/query first")

    ra = req_data.get("ra", 0.0)
    dec = req_data.get("dec", 0.0)
    fits_paths = [f["fits_path"] for f in frames]

    try:
        result = run_candidate_pipeline(fits_paths, ra, dec)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Pipeline error: {exc}")

    # Attach frame timing metadata and cross-match candidates
    from backend.catalog_service import classify_candidate_target
    for c in result.get("candidates", []):
        c["target_ra"] = ra
        c["target_dec"] = dec
        c_ra = c.get("ra_deg", ra)
        c_dec = c.get("dec_deg", dec)
        try:
            cat = crossmatch(c_ra, c_dec, radius_arcsec=10.0)
            c["catalog_classification"] = cat.get("classification")
            c["classification"] = classify_candidate_target(c, cat)
            c["catalog_summary"] = cat.get("summary")
        except Exception:
            c["catalog_classification"] = "unqueried"

    with _lock:
        _session["candidates"] = result

    return result


@app.get("/api/diff/template.png")
def diff_template():
    with _lock:
        diff = _session.get("diff")
    if not diff:
        raise HTTPException(status_code=400, detail="run /api/process first")
    png = Path(diff["template_png"])
    if not png.exists():
        raise HTTPException(status_code=404, detail="template missing")
    return FileResponse(png, media_type="image/png")


@app.get("/api/diff/{index}.png")
def diff_frame(index: int):
    with _lock:
        diff = _session.get("diff")
    if not diff:
        raise HTTPException(status_code=400, detail="run /api/process first")
    png = Path(diff["dir"]) / f"diff_{index:04d}.png"
    if not png.exists():
        raise HTTPException(status_code=404, detail="diff frame missing")
    return FileResponse(png, media_type="image/png")


class CatalogRequest(BaseModel):
    ra: float = Field(..., ge=0, le=360)
    dec: float = Field(..., ge=-90, le=90)
    radius_arcsec: float = Field(10.0, gt=0, le=120)
    epoch: float = Field(2025.5)


@app.post("/api/catalog")
def catalog(req: CatalogRequest):
    """Cross-match position against Gaia DR3 + AllWISE."""
    return crossmatch(req.ra, req.dec, req.radius_arcsec, req.epoch)


@app.get("/api/presets")
def presets():
    return {"presets": [
        {"name": "WISE 0855", "ra": 133.795129, "dec": -7.245147,
         "note": "Y-dwarf, Phase 0 verified target"},
    ]}


# Static frontend (mounted last so /api wins)
if FRONTEND_DIR.exists():
    app.mount("/", StaticFiles(directory=str(FRONTEND_DIR), html=True),
              name="frontend")
