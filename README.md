# MOLLER detector workspace

Load remoll geometry and ROOT files; inspect hits, PE signals and deconvolution.
Files stay on your computer.

In the terminal where ROOT and Geant4 are configured:

```sh
python3 launch.py --remoll "/path/to/remoll"
```

On Windows, use `py -3` instead of `python3`, from your configured compiler/ROOT
[terminal](https://root.cern/install/). Open http://127.0.0.1:8784/; **Quit** or Ctrl+C stops it.
No Python packages, Node, Docker, remoll rebuild or internet connection needed.
Requires **Python 3.9+, ROOT with ACLiC/C++17, Geant4 with GDML, CMake 3.20+ and a compatible C++17 compiler**.
Use the compiler/runtime environment matching your existing ROOT and Geant4 builds.

- `--check`: check installed tools; installs nothing.
- `--geometry-only`: build the helper, write geometry JSON, then exit.
- `--output /path/to/build`: choose a build/result directory without spaces (ROOT/ACLiC limitation; input filenames may contain spaces).
- `--geant4-dir /path/to/Geant4`: folder containing `Geant4Config.cmake`, if not found.
- `--no-open`: serve without opening a browser.

Geant4 discovery uses `Geant4_DIR`, remoll's `build*/CMakeCache.txt`, or CMake's
[normal CMake search](https://geant4.web.cern.ch/documentation/dev/ig_html/InstallationGuide/buildtools.html) (`CMAKE_PREFIX_PATH`). The first geometry load compiles one small
helper; subsequent loads reuse it. Native file dialogs are optional: paths can be
pasted directly. Linux dialogs use Zenity or Tk if already installed.

| Step | Code |
|---|---|
| Launch, build and serve | `launch.py` → `geometry/run.py` (Python standard library) |
| GDML → placed detector meshes | `geometry/export_geometry.cc` (Geant4, built by CMake) |
| ROOT files → cuts, histograms, PE | `geometry/load_hits.C` (ROOT/ACLiC) |
| ROOT files → deconvolution | `Deconvolution.C` → `root/Deconvolution.C` |
| Results → browser | `geometry/deconvolution_bridge.C`, JSON, `geometry/*.js` |
| Response models | `app/responses.py`, `data/`, `vendor/` |
| [Visual guides](docs/GUIDES.txt) (PDF) | [Loading](docs/run-flow.pdf) · [Analysis](docs/analysis-flow.pdf) · [Using the app](docs/workspace-guide.pdf) · [ROOT fields](docs/root-output-atlas.pdf) |

See [ROOT-only commands](docs/ANALYSIS.md) and
[response origins](THIRD_PARTY.md). Bundled or local PE models are selectable.
PE/current requires recorded detector entries; all-step files remain usable in Hits.

**Visual guides:** see [docs/](docs/GUIDES.txt) for diagrams of the code, analysis and app, with example cuts.

**Tests:** see [tests/README.md](tests/README.md) for what each test checks and how to run it.
Run all: `python3 tests/run.py --remoll "/path/to/remoll"`.

Standalone builds/results stay in `build/`; this research checkout uses
`rootfiles/local_work/experiments/detector_launcher/geometry-build`.
`MOLLER_OUTPUT` also overrides the location. **Save view** downloads one JSON;
recent selections use bounded memory, with no automatic saved-analysis archive.
ROOT/GDML and local response code run natively; use trusted files and code.

Uses standard CMake and ROOT build paths for macOS, Linux and Windows; the latest
end-to-end source-package check was on macOS, not a three-platform certification.
