# PROJECT_STATUS.md

## Pre-Phase-5 — Bug Fixes, UI Polish & Offline Fallback ✅

| Fix | Status | Notes |
|-----|--------|-------|
| A. Vertical orientation | ✅ | Verified via WCS: FITS array row 0 = northernmost declination → both `render.py` and `diff_imaging.py` write north-up PNGs identically; Original/Processed views share orientation, no flip possible |
| B. Network resilience | ✅ | SIA2 query: 3 retries w/ backoff; `_fetch` already retries 4× (covers WinError 10060/timeouts). **Offline fallback**: `session_manifest.json` persisted after each successful query + rebuilt from FITS headers (`scripts/rebuild_manifest.json` — cutout headers carry MJD-AVG/OBSID/DETECTOR). On IRSA failure → serves cached frames with `offline:true`. Verified live: IRSA 502 → 200 offline response, 133 frames, process+diff all work |
| C. VizieR fallback | ⚠️→✅ | Root cause = CDS-wide outage (tapvizier connect-timeout, SIMBAD 404, VizieR 404) — not a URL bug. Fallbacks retained; catalog response cache added (`catalog_cache/`, 7-day TTL, serves last good result when live queries fail) |
| D. Catalog card UI | ✅ | Spinner during lookup, color-coded badges, card-style match entries, warnings/disclaimer blocks |
| E. Detector explainer | ✅ | ⓘ note in metadata panel explaining D1–D6 array changes across epochs |

**Live verification (2026-11-15, IRSA degraded)**: `/api/query` → 200 `offline:true` (133 frames from cache); `/api/process` → 200 (D2, 27 residual frames); `/api/diff/*` + `/api/frame/*` → 200. App fully functional with zero external connectivity for the WISE 0855 demo path.

## Phase 4 — Catalogs, Cross-Matching & Metadata ✅ (degraded-service verified)

| Task | Status | Notes |
|------|--------|-------|
| `backend/catalog_service.py` | ✅ | Gaia DR3 (`astroquery.gaia` cone search) + AllWISE (`allwise_p3as_psd` via IRSA TAP) + VizieR fallbacks (I/355/gaiadr3, II/328/allwise); proper-motion propagation to epoch 2025.5; 3″ match tolerance; quality flags (RUWE, astrometric_excess_noise, cc_flags, ph_qual) |
| `POST /api/catalog` | ✅ | Returns classification (`known_object` / `unmatched_within_tolerance` / `discovery_candidate`), matches, flags, warnings, disclaimer |
| Resilience | ✅ | Queries run in parallel; 150 s timeout each; VizieR fallback; all failures → warnings, never 5xx |
| Frontend catalog panel | ✅ | Sidebar panel: status badge (Known Object / Candidate), match cards (sep, mags, pm, parallax, flags), warnings, disclaimer |
| Live verification | ⚠️ | **External catalog services were degraded during testing** (Gaia TAP 500/timeout, IRSA TAP 502/hang, VizieR+SIMBAD 404 — CDS-wide outage). Endpoint correctly returned `discovery_candidate` + warnings in ~4 min. Gaia path proven earlier: cone_search returned 1 source at WISE 0855 position. Re-test when services recover. |

## Phase 3 — Image Processing, Alignment & Difference Imaging ✅ (with caveats)

| Task | Status | Notes |
|------|--------|-------|
| `processing/diff_imaging.py` | ✅ | FITS stack loader, FFT phase-correlation alignment (subpixel + 2px sanity gate), median template, residual frames, diverging-colormap PNG render |
| `POST /api/process` | ✅ | Builds diff stack for session frames; auto-selects dominant detector (cross-detector mixing is invalid — independent distortion per array) |
| `GET /api/diff/{i}.png`, `/api/diff/template.png` | ✅ | Verified 200s |
| Frontend toggle | ✅ | "Original" / "Processed (Difference · D5)" segmented control; separate frame lists, slider/playback work in both modes; diff stack auto-builds after load |
| Bad-frame handling | ✅ | Corrupt/misshapen FITS → NaN placeholder, index alignment preserved |

**Verification (WISE 0855)**: 133 frames queried; D5 stack = 27 residual frames. Static suppression works mechanically (residual std ~77% of raw for D5). **Caveat (see docs/decisions/0001)**: phase-correlation shifts were mostly noise-driven (~8px median) and rejected; cutouts are already sky-registered, so effective alignment ≈ none. WISE 0855 motion is sub-pixel — it will not appear as a mover; residuals show noise/coverage/artifacts. True mover detection needs larger FOV / PSF matching / stacking (P1).

## MVP Readiness

P0 demo path works end-to-end: coordinates → IRSA query → cached cutouts → PNG animation → difference view. Ready for human review; difference-imaging depth is honestly limited (documented).

## Phase 1 — Backend Data Access ✅

| Task | Status | Notes |
|------|--------|-------|
| `backend/data_service.py` | ✅ | `get_spherex_cutouts(ra, dec, fov_size)` — SIA2 query over DR1+QR3 (wide+deep), threaded cutout downloads, local cache in `backend/cache/` |
| `processing/render.py` | ✅ | FITS IMAGE HDU → PNG (asinh stretch + colormap, 512px) |
| `backend/main.py` (FastAPI) | ✅ | `POST /api/query`, `GET /api/frame/{i}.png`, `GET /api/presets`, `GET /api/health`; serves frontend static |
| WISE 0855 backend test | ✅ | 270 SIA rows / 137 epochs; **132 frames downloaded+rendered** (~500 MB FITS, ~200 KB PNG each); 5 frames failed transiently (retried once — acceptable loss). Epochs span 2025-05-05 → 2025-11-14 |
| Cache | ✅ | `backend/cache/` (gitignored via `cache/`); re-queries reuse FITS/PNGs |

## Phase 2 — Frontend Basic Viewer ✅

| Task | Status | Notes |
|------|--------|-------|
| `frontend/` UI | ✅ | Vanilla HTML/CSS/JS (chose over React to keep single-server deploy; Node exists if we migrate later) |
| RA/Dec/FOV/max-frames inputs + WISE 0855 preset | ✅ | `index.html`, `app.js`, `style.css` |
| Viewer + metadata panel | ✅ | PNG frame, obs_id, detector, bandpass, ISO timestamp |
| Animation controls | ✅ | Play/Pause, frame slider, speed select (0.5–4×), frame counter |
| Live test | ✅ | Server running at http://127.0.0.1:8000; `/api/query` → 132 frames; `/api/frame/{i}.png` 200s; browser preview opened |

## Phase 0 — Preparation, Research & Scientific Feasibility

Status legend: ✅ done · 🔄 in progress · ⬜ pending · ⚠️ blocked/risk

| Task | Status | Notes |
|------|--------|-------|
| Directory structure | ✅ | frontend/, backend/, processing/, tests/, scripts/, sample-data/, docs/, phase-0/experiments/ |
| Root docs (AGENTS/README/requirements/.gitignore) | ✅ | Created |
| Python .venv + package install | ✅ | Python 3.12 venv; astropy 8.0.2, astroquery 0.4.11, pyvo 1.10.0, etc. |
| verify_env.py | ✅ | All 7 packages import OK |
| IRSA official-source research docs | ✅ | phase-0/research/{official_sources,spherex_data_structure,irsa_access_methods}.md |
| WISE 0855 query experiment | ✅ | **258 MEF rows / 129 epochs**, all spherex_dr1; see experiment_report.md |
| FITS sample download + HDU inspection | ✅ | wise0855_first_cutout.fits, 4.05 MB, 29×29 px, all 7 HDUs verified |
| Architecture docs (tech_stack, mvp_scope, team_roles) | ✅ | docs/architecture/, docs/mvp_scope.md |
| experiment_report.md | ✅ | phase-0/experiments/experiment_report.md |

## Phase 0 Exit Criteria
- [x] `.venv` exists and `verify_env.py` passes all imports
- [x] Research docs cite only official IRSA/SPHEREx documentation
- [x] WISE 0855 query returns real results — 129 epochs confirmed
- [x] FITS sample downloaded and headers/HDUs inspected
- [x] Architecture + MVP scope docs written
- [x] This file updated — **ready for human review**

## Key Findings (for reviewers)
- IRSA SIA2 collections live: `spherex_dr1`, `spherex_dr1_deep`,
  `spherex_dr1_cube`, `spherex_qr3`, `spherex_qr3_deep`.
- WISE 0855: 129 distinct DR1 epochs (t_min 60800.7→61174.3), all 6 detector
  arrays. `energy_bandpassname` = detector ID (D1–D6), NOT spectral band.
- Cutout service works: `?center=RA,Dec&size=deg` on any MEF `access_url`.
  Cutouts preserve all extensions (IMAGE/FLAGS/VARIANCE/ZODI/EPSF/WCS-WAVE).
- Transient HTTP 503 seen once on cutout download; use retries + UA header.

## Open Questions
- QR3 returned zero rows for this position — need to check whether QR3
  coverage is still being ingested or genuinely absent here.
- Consecutive frames often share `t_min` across different detector arrays
  (D1–D6) — decide whether the animation should dedupe to unique epochs or
  keep all detector frames.
- 5/137 frame downloads failed transiently in the big test — add a retry
  pass in `get_spherex_cutouts` for missing frames.
- t_min confirmed as MJD (converted to UTC ISO in API responses).
- Which spectral channel/band to default to in the viewer (WCS-WAVE is a
  per-column lookup; needs a small extraction helper).

## Known Risks
- Machine has no `git` on PATH — version control deferred.
- SPHEREx ingest is weekly; result counts will grow — cache accordingly.
- WISE 0855 is faint; pixel-level detection significance unverified
  (needs photometry — P1, not Phase 0).

## Log
- 2026-11-14 — Phase 0 COMPLETE. Structure + docs + venv done; IRSA research
  documented; WISE 0855 query returned 258 rows/129 epochs; 4 MB cutout FITS
  downloaded and inspected; architecture/MVP/team docs written.
- 2026-11-14 — Phase 1+2 COMPLETE. Backend `data_service.py` + FastAPI app
  verified on WISE 0855 (132 frames cached); vanilla-JS frontend viewer with
  animation controls live at :8000; full 137-frame load took ~5 min (downloads
  dominate; caching makes re-views instant).
- 2026-11-14 — Phase 3 COMPLETE. `processing/diff_imaging.py` + `/api/process`
  + `/api/diff/*` endpoints + frontend Original/Processed toggle. Detector-
  filtered median-subtraction stack (D5, 27 frames) live. Alignment caveat
  documented in docs/decisions/0001-difference-imaging.md.
- 2026-11-14 — Phase 4 COMPLETE (service-degraded). `catalog_service.py` +
  `/api/catalog` + frontend catalog panel. Gaia/AllWISE/VizieR all timed out
  or 404'd during verification — graceful-degradation path verified; real
  matches pending service recovery. CDS (VizieR/SIMBAD) returned site-wide
  404s tonight.
- 2026-11-15 — Pre-Phase-5 fixes COMPLETE. Orientation verified north-up;
  offline cached fallback implemented + verified during IRSA outage
  (manifest rebuilt from FITS headers: 133 frames); catalog cache + UI
  spinner/cards/detector note added.
