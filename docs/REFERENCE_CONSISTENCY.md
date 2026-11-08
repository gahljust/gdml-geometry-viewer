# Reference consistency investigation — 2026-11-08

## Correction

Fixed the **legacy deconvolution angular cut** to use physical Ring 5 Open
tiles for the **asymmetry reference in H**. Both main weighting modes retain
the rate-weighted, charged-crossing reference and the configured kinetic cut.
No detector signal cuts, PE tables, event weights or GLS equations changed.
The shared campaign website already uses physical membership and was untouched.

## Why the previous checks missed it

1. **Wrong comparison target.** `tests/verify_root.py` compared fresh app results
   to this launcher's saved `latest.json.gz`, not to the current campaign's
   independently produced data. This established regression agreement with the
   older implementation, not agreement with the current analysis definition.
2. **Solver agreement started too late.** Its Python comparison consumed the
   app scanner's captured moments. Both solvers therefore received the same
   reference selection. Correct linear algebra cannot detect that selection
   disagreement.
3. **The definition was known but not reconciled.** The October 25
   `archive/viewers/detector_launcher_legacy/docs/REGION_MEMBERSHIP_AUDIT.md`
   explicitly identified the angular Ring 5 reference, treated it as a historical
   convention and left it unchanged. The campaign's physical-region correction
   was not propagated to this active launcher. This was a consistency failure,
   not an undiscovered mathematical limitation of correlated fits.
4. **The recorded contract was incomplete.** The old verification contract listed
   particle cuts, weights and channel counts but omitted the reference-region
   selection. Tests did not require agreement on the entire raw-input analysis
   with the accepted campaign definition. The normal nine-test package suite
   did not run the real-data comparison at all.
5. **The baseline was mutable.** `latest.json.gz` is a working output; it has
   since been replaced by an unrelated one-file geometry-viewer result. The
   historical receipt preserves its original baseline hash, but that hash was
   not enforced on later invocations. A working output is not a test oracle.

The prior statement that the software “exactly matched” was too broad. The
recorded tests demonstrate agreement with the selected old baseline and solver
algebra, not consistency with the current shared viewer's complete analysis.

## Replacement check and evidence

`tests/verify_consistency.py` independently runs the app's public ROOT entry and
the standalone macro authenticated by its full-dataset validation receipt.
Both read identical, hash-checked raw ROOT inputs and intended settings. The
reference uses accepted campaign PE tables supplied independently of the app's
bundled tables, so table drift is not hidden by giving both sides the app's copy.
The reference macro was separately checked against the accepted 9M and 29.79M
campaign moments; the app under test does not generate its expected results.

The new check compares detector order, complete rate/asymmetry sums, full joint
MC input covariance, reference asymmetries, rates, F, its MC errors, H, B, C,
correlations, M, and experimental/MC/combined fit covariances in both modes.
It records source/table/input hashes and rejects changes during execution.
It leaves a compact receipt and removes temporary scans and compiled outputs.

Executed on the original validation sample: **375,000 histories, three files,
all 252 channels**, main crossings + ShowerMax PE and main PE + ShowerMax PE.

| Stage | Before correction | After correction |
|---|---|---|
| Detector order, rates, F, F MC errors, B, C, correlations | Agree | Agree |
| Rate/asymmetry sums and their MC covariance | Disagree | Agree |
| Reference asymmetries, H, M and all fit covariances | Disagree | Agree |

**30 comparisons pass after correction.** Before correction, 16 fail. This is
a general pipeline check that caught the defect, not a test of a hard-coded
Ring 5 example. Reference asymmetries agree exactly after correction; other
objects are compared numerically with recorded tolerances, not bitwise promises.

Local evidence under `rootfiles/local_work/`:

- `experiments/detector_launcher/reference_consistency/before.json`
- `experiments/detector_launcher/reference_consistency/after.json`
- `experiments/detector_launcher/reference_consistency/release-check.json`
  (final general check through the package runner, independent campaign tables).
- Historical: `experiments/detector_launcher/adversarial-validation.json`,
  `verification-contract.json`, `VERIFICATION.md` (preserved, not rewritten).
- Independent reference: `validation/correlated_deconvolution/20261108/validation.json`.

## Required verification scope

The package runner now exposes this comparison as `--only consistency`; an
ordinary run explicitly reports it untested without the required receipts.
Analysis releases require it in addition to the known-answer tests. Historical
regression now requires an explicit immutable baseline and a new output receipt;
it cannot silently use `latest.json.gz` or overwrite prior evidence.

This verifies consistency for the inputs and supported modes executed, not every
possible configuration. Additional generators, unsupported stratified samples
and unrepresented boundaries remain outside this run. Agreement between related
implementations does not replace independent mathematical tests. No simulation,
full-cohort reanalysis, repository cleanup or publication was performed.

All nine existing package checks also passed. The HTTP hit-loader check initially
could not bind localhost in the restricted environment; it passed when rerun with
local-server access. The final consistency check passed separately through the
same runner. Temporary test servers, builds and ROOT products were removed.
