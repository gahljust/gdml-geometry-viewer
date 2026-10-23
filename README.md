# GDML Geometry Viewer

A small, local viewer for remoll’s six main detector rings and ShowerMax.
Watch Geant4 placements become meshes, inspect the executing code, and unwrap
all detector rings. Optionally color them by ROOT hit counts or stored rates.

## Run

Requires Python 3.10+, CMake 3.20+, a C++17 compiler, and Geant4 with GDML support.
ROOT is required for ROOT-file analysis. Use an environment where these tools
are available; no Docker or Python/JavaScript packages are needed.
Linux’s file picker needs Zenity or Python Tk.

```bash
git clone https://github.com/gahljust/gdml-geometry-viewer.git
cd gdml-geometry-viewer
python3 run.py
```

Select your **remoll checkout**, then **Load**. The first load builds a small
helper in `build/`. If Geant4 is not found, set `Geant4_DIR` to the directory
containing `Geant4Config.cmake`. On Windows, use `python` instead of `python3`.

**Pause / Step** follows the code; **Finish loading** skips the waits.
**Stop** cancels the current operation; **Reload** starts again.
**Quit** or **Ctrl+C** closes the app. Default address: <http://127.0.0.1:8785/>.
Use `--port 8786` for another port, or `--no-open` to skip opening the browser.

## ROOT values

Choose a file, then click **Analyze**. Files remain local.

- **Hits:** every stored `hit.det` record; no particle, energy or direction cuts.
- **Rate [Hz]:** sum the stored scalar `rate` once per hit, without renormalizing.
- Main detector IDs match GDML `DetNo`. ShowerMax uses **virtual planes only**;
  its physical quartz layers are display geometry, never a substitute for plane hits.
- Missing plane hits make ShowerMax unavailable. Missing/invalid rates disable Rate.
  Zero hits cannot distinguish a disabled detector from one that received no hits.

Uses tree `T`, split `hit.det` (or `vector<int>`), and optionally scalar double
`rate`. This is a hit/rate viewer, not PE response or deconvolution analysis.
Geometry comes from the chosen checkout, before macro-driven moves; the app
cannot verify that a ROOT file used that same geometry.

## Code

| File | Purpose |
| --- | --- |
| `run.py` | Local server, embedded Geant4 exporter and ROOT counting macro |
| `index.html` | Canvas viewer, placement trace, unwrapping and hit display |
| `test_*.py`, `test_display.js` | Lifecycle, synthetic ROOT inputs and display checks |

Geant4 reads `geometry/mollerMother.gdml` and its references. The exporter walks
the loaded detector hierarchy, combines parent transforms and checks vertices
against Geant4 navigation. The native reader waits at each displayed step;
this shows the loaded placements, not the parser’s internal operations.
Generated C++ and its displayed line numbers can be inspected in `build/`.
ShowerMax plane IDs come from `geometry/mollerParallel.gdml`.

## Checks

With the app running, set `REMOLL_DIR` to your remoll checkout, then run:

```bash
python3 test_lifecycle.py
python3 test_root.py
node test_display.js
```

Tests replace the current load; ROOT fixtures are small and removed afterward.
The ROOT test uses the checkout’s `develop` header. Node is needed only for the
display test. Tested on macOS and [Linux](https://github.com/gahljust/remoll/actions/runs/39910911862);
Windows has not been verified.
