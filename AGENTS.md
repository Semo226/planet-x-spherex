# AGENTS.md — Project Operating Rules

## Project
NASA Space Apps 2026 — "Planet X and SPHEREx".
Goal: Build a tool that lets users query SPHEREx IRSA archives for a target
coordinate, fetch/visualize observations over time (light-curve / epoch
animation), with WISE 0855 (a cold brown dwarf near the ecliptic) as the
flagship scientific test case.

## Core Working Rules
1. **Scientific honesty**: Never invent API endpoints, table names, or data
   values. Every query/API claim must be verified against official NASA/IPAC
   IRSA documentation or by running a real query. If a value is an estimate,
   label it as such.
2. **Official sources only for data access**: Use `astroquery.irsa`, IRSA
   SIA2 (`/sva/sia2api`), and IRSA TAP (`/tap`) services as documented in the
   official IRSA tutorials. Do not guess REST paths.
3. **Reproducibility**: All experiment scripts live under `phase-0/experiments/`;
   query outputs are saved as CSV/FITS artifacts next to the scripts.
4. **Large data**: Do not commit FITS files or caches. Sample data goes in
   `sample-data/` and is limited to small cutouts.
5. **Documentation discipline**: Record design decisions in
   `docs/decisions/`, architecture in `docs/architecture/`, and research
   findings in `phase-0/research/`.
6. **Phase gating**: Do not start implementation (frontend/backend features)
   until Phase 0 is complete and reviewed. Phase 0 = preparation, research,
   and scientific feasibility only.
7. **Secrets**: No API keys or credentials in the repo. IRSA is public; if any
   credential is ever needed, use environment variables only.

## Environment
- Python 3.12, project-local virtualenv at `.venv/`.
- See `requirements.txt` for pinned scientific packages.
- Verification: `python phase-0/experiments/verify_env.py`

## Status
See `PROJECT_STATUS.md` for the live phase tracker.
