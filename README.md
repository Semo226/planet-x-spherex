# Planet X and SPHEREx

NASA Space Apps 2026 project. We build an interactive explorer for NASA's
SPHEREx all-sky infrared survey data, focusing on epoch-to-epoch variability
of faint, cold objects — starting with WISE 0855, a Y-dwarf ~7.2 pc away that
lacks strong SPHEREx coverage and is a stand-in for undiscovered
"Planet X"-type bodies.

## Team
- TBD — placeholder for the 6-member team roster (see `docs/architecture/team_roles.md`).

## Credits
- Data: NASA/IPAC Infrared Science Archive (IRSA) — SPHEREx survey.
- Challenge: NASA Space Apps 2026 "Planet X and SPHEREx".
- Built with Python (astropy, astroquery, FastAPI) and a React frontend.

## Repository Layout
- `frontend/` — React web UI (epoch animation viewer, coordinate input)
- `backend/` — FastAPI service wrapping IRSA queries
- `processing/` — FITS processing / cutout / light-curve pipeline
- `tests/` — automated tests
- `scripts/` — utility scripts
- `sample-data/` — small FITS cutouts used for demos
- `docs/` — architecture, decisions, research notes
- `phase-0/experiments/` — Phase 0 feasibility experiments and artifacts

## Quickstart (Phase 0)
```powershell
py -3.12 -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python phase-0/experiments/verify_env.py
python phase-0/experiments/wise0855_query.py
```
