# Tests

These tests check that the app calculates known answers, reads ROOT files correctly,
applies cuts and PE weights consistently, and rejects invalid inputs. They also
check loading, saved views and the downloadable package. Small synthetic ROOT
files provide controlled inputs, so you do not need a production dataset.

Use them after installation or code changes. Passing means the tested cases work;
it does not establish that your simulation settings or detector response models
are physically correct.

Run from the package directory, with ROOT configured:

```sh
python3 tests/run.py --remoll "/path/to/remoll"
```

Windows: use `py -3` instead of `python3`. Add `--verbose` for details.
Node.js is needed for `hits` and `saved-view` tests, not for the app.
The `pe` and `hits` tests need a remoll checkout (`--remoll` or `REMOLL_DIR`).

## Run one test

```sh
python3 tests/run.py --remoll "/path/to/remoll" --only solver
```

Replace `solver` with any name below. Each runs the listed file in `tests/`.

| `--only` | Test file | Checks |
|---|---|---|
| `solver` | `test_root_solver.C` | Known statistical σ = 1–5 ppb and MC σ = 0.1–0.5 ppb; invalid and singular matrices. |
| `root-inputs` | `test_root_inputs.C` | Synthetic ROOT metadata, incompatible runs, detector settings and malformed PE tables. |
| `responses` | `test_response_sources.py` | Bundled/local response resources and invalid inputs. |
| `cache` | `test_result_cache.py` | In-memory reuse, invalidation and no automatic disk cache. |
| `launch` | `test_launch.py` | Dependency discovery, paths and launcher options. |
| `package` | `test_package.py` | Exported files, hashes, launch command and exclusion of local artifacts. |
| `pe` | `test_weighted_maps.py` | Known 90 + 66 = 156 PE/s through ROOT; multi-file normalization and input rejection. |
| `hits` | `test_hit_loader.py` | Cuts, histogram bins, units, malformed files and recovery. |
| `saved-view` | `test_save_view.js` | Saved-view contents and display units. |

Full suite: **9 passed, 0 failed, 0 skipped** when prerequisites are available.
A failure exits nonzero; a skip means untested. Temporary files in
`build/checks-*` are removed afterward. Use `--work-dir /path/without/spaces`
to change that location. The export test uses a temporary sibling directory,
also removed afterward.

The solver's pure-component case uses H = I and C = diag(1,4,9,16,25) ppb².
These checks verify known cases and individual stages, not the physics accuracy
of arbitrary inputs. PE tests require geometry compatible with the bundled model.

## Geometry check

Start a separate test instance (do not use your active work session):

```sh
python3 launch.py --remoll "/path/to/remoll" --port 8786 --no-open
```

In another configured terminal:

```sh
python3 tests/test_geometry.py http://127.0.0.1:8786 "/path/to/remoll"
```

Checks 224 main tiles, 28 ShowerMax planes, valid geometry and repeat/error loads.
Stop the test server with Ctrl+C.
