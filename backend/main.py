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

app = FastAPI(title="Planet X and SPHEREx API", version="0.2.0")

FRONTEND_DIR = ROOT / "frontend"
CACHE_DIR = Path(__file__).resolve().parent / "cache"

# In-memory "session": latest query result. Fine for single-user hackathon.
_lock = threading.Lock()
_session: dict = {"frames": [], "diff": None}


class QueryRequest(BaseModel):
    ra: float = Field(..., ge=0, le=360)
    dec: float = Field(..., ge=-90, le=90)
    fov_size: float = Field(0.05, gt=0, le=1.0)
    max_frames: int | None = Field(60, ge=1, le=300)


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
        _session = {"frames": result["frames"], "diff": None,
                    "request": {"ra": req.ra, "dec": req.dec,
                                "fov_size": req.fov_size}}
    return JSONResponse({
        "ra": result["ra"], "dec": result["dec"],
        "fov_size": result["fov_size"],
        "n_results": result["n_results"],
        "n_epochs": result["n_epochs"],
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
    method: str = "median"   # 'median' (only method implemented so far)
    detector: str | None = None  # e.g. 'D3'; None = auto (most common)


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
        (detector + fits_paths[0] + fits_paths[-1]).encode()).hexdigest()[:12]
    out_dir = CACHE_DIR / "diff" / tag
    result = build_difference_stack(fits_paths, out_dir, align=req.align)

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
                            "frames": diff_frames}
    return {
        "detector": detector,
        "n_frames": result["n_frames"],
        "limits": result["limits"],
        "template_png": "/api/diff/template.png",
        "frames": diff_frames,
    }


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
