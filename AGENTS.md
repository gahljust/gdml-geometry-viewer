# MOLLER Geometry Viewer

- Entry: `launch.py` → `geometry/run.py`. Browser: `geometry/index.html`, `view.js`, `deconvolution.js`.
- GDML: `geometry/export_geometry.cc`. Hits/cuts: `geometry/load_hits.C`.
- Fit: `root/Deconvolution.C`, `root/root_solver.hh`; bridge: `geometry/deconvolution_bridge.C`.
- Manual view export: `geometry/save-view.js`; no automatic analysis storage.
- PE resources: `app/responses.py`, `data/response-manifest.json`, `vendor/`.
- Visual guides: `docs/GUIDES.txt`, PDFs in `docs/`. Editable masters in `docs/visual/` stay local; preserve manual edits.
- Checks: `python3 tests/run.py --remoll /path/to/remoll`; scope in `tests/README.md`.
- Outputs follow `MOLLER_OUTPUT`; never store ROOT/builds in source.
- Keep text brief and the interface within one frame. No automatic publication.
