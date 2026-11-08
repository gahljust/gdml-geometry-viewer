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
| `solver` | `test_root_solver.C` | Known statistical σ = 1–5 ppb and MC σ = 0.1–0.5 ppb; stable sparse moments; invalid and singular matrices. |
| `root-inputs` | `test_root_inputs.C` | Synthetic ROOT metadata, detector settings and PE tables; constant/near-constant scores, zero-hit events and inelastic cross-category covariance. |
| `responses` | `test_response_sources.py` | Bundled/local response resources and invalid inputs. |
| `cache` | `test_result_cache.py` | In-memory reuse, invalidation and no automatic disk cache. |
| `launch` | `test_launch.py` | Dependency discovery, paths and launcher options. |
| `package` | `test_package.py` | Exported files, hashes, launch command and exclusion of local artifacts. |
| `pe` | `test_weighted_maps.py` | Known 90 + 66 = 156 PE/s; actual energy tables through maps and deconvolution, bin boundaries, kinetic cuts, multi-file normalization and input rejection. |
| `hits` | `test_hit_loader.py` | Cuts, histogram bins, units, malformed files and recovery. |
| `saved-view` | `test_save_view.js` | Saved-view contents and display units. |
| `consistency` | `verify_consistency.py` | Fresh raw ROOT scans through the app and a separately validated reference; channel order, all rate/asymmetry sums, joint MC covariance, reference asymmetries, F, H, B, C, M, correlations and all fit uncertainties, in both weighting modes. |

The nine package checks can pass without establishing agreement with the current
campaign. The consistency check is explicitly **skipped** unless its reference
and raw-input receipts are supplied. A skip is not release verification.
A failure exits nonzero; a skip means untested. Temporary files in
`build/checks-*` are removed afterward. Use `--work-dir /path/without/spaces`
to change that location. The export test uses a temporary sibling directory,
also removed afterward.

The solver's pure-component case uses H = I and C = diag(1,4,9,16,25) ppb².
These checks verify known cases and individual stages, not the physics accuracy
of arbitrary inputs. PE tests require geometry compatible with the bundled model.

## Analysis consistency before release

```sh
python3 tests/run.py --only consistency \
  --reference-macro /path/to/validated/deconvolution_correlated.C \
  --reference-validation /path/to/validation.json \
  --input-receipt /path/to/identified-inputs.json \
  --reference-main-tables /path/to/accepted/main-tables \
  --reference-sm-tables /path/to/accepted/showermax-tables \
  --consistency-receipt /path/to/new-receipt.json
```

Requires NumPy and PyROOT. The validation receipt must bind the reference's SHA256
to a passed validation. The input receipt contains `inputs` with `path`, `sha256`
and `events` for each raw file, covering Møller, elastic and inelastic. Both
implementations run independently from those files, with the same intended
current, exposure and cuts. The reference uses independently supplied accepted PE
tables; the app uses its bundled tables, so table drift is also detectable.
Input identities and all analysis-source/table hashes
are recorded. Changed sources or inputs invalidate the run. Temporary builds and
ROOT outputs are removed; one new receipt remains. No reference is regenerated
from the app under test. Investigate discrepancies instead of updating expected
answers to make them pass.

This is general cross-implementation consistency, not a special Ring 5 test.
It covers the inputs and modes executed; extra generators, sampling schemes
unsupported by this app and unrepresented boundaries require separate coverage.
Two implementations can share a mistake, so the independent known-answer solver
and ROOT tests remain necessary. `verify_root.py` is now explicitly historical
regression only and requires an immutable baseline; it cannot use mutable
`latest.json.gz` or replace this release check.

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
