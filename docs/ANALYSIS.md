# Analysis entry points

| Step | Code |
|---|---|
| GDML placements → detector faces | `geometry/export_geometry.cc` |
| ROOT branches, cuts, histograms | `geometry/load_hits.C` |
| PE resources and local overrides | `app/responses.py`, `root/pe_response.hh` |
| Input metadata, event sums, fit | `root/Deconvolution.C` |
| Covariance and GLS | `root/root_solver.hh` |
| ROOT results → browser JSON | `geometry/deconvolution_bridge.C` |

Detector selections use recorded IDs. ShowerMax signals use virtual planes
73001–73028, not internal quartz layers. Presets expose their ROOT cuts in the UI.
PE model support restrictions still apply. Cathode current is PE/s × elementary charge.

For one event, add the selected expected PE in each detector. Call the totals
in two detectors `a` and `b`, and the normalized event rate `w`. Accumulate:

- Detector rates: `R₁ = Σ w a`, `R₂ = Σ w b`.
- Second moments: `B₁₁ = Σ w a²`, `B₁₂ = B₂₁ = Σ w a b`, `B₂₂ = Σ w b²`.

Each sum runs over events, using that event’s `a`, `b` and `w`. The rate appears
once in each product. Repeat for all detector pairs to fill B.

Convert each entry to projected asymmetry covariance:
`Cᵣₛ = Bᵣₛ / (P² T Rᵣ Rₛ)`; multiply by `10¹⁸` for ppb².
Here r and s identify detectors, P is polarization and T is beam time in seconds.
Use the asymmetry-adjusted dilution matrix H to calculate `M = Hᵀ C⁻¹ H`.
The projected statistical uncertainty of each extracted asymmetry is the square
root of its diagonal entry in `M⁻¹`. Finite simulation uncertainty propagates
through H separately.

These are conditional forecasts. C is estimated from simulation and then treated
as fixed; its estimation uncertainty is not propagated. MC through H uses a
first-order approximation, which can fail when MC variations make the components
indistinguishable. Numerical invertibility alone does not establish reliable
uncertainty intervals. Expected PE products include variation between simulated
histories, but omit additional optical, photoelectron-counting and electronics
fluctuations. Combined adds only the displayed statistical and linear-MC terms;
it is not a complete experimental error budget. Production impact of these
approximations requires separate validation.

Deconvolution requires compatible run metadata and detector coverage.
Missing information is reported, not treated as zero. Compatibility compares saved
source, geometry, field hashes and ordered transport commands, including macro
calls. Included macro contents are not saved in ROOT; matching calls alone cannot
verify later edits to those files. Progress printing and known rate-corrected
beam-momentum proposal settings do not change transport identity. Neyman/fixed-quota
samples remain unsupported: their MC uncertainty requires stratified statistics.

The same steps can run separately, without the browser:

```sh
python3 launch.py --remoll "/path/to/remoll" --geometry-only --output build
root -n -l
```

At the ROOT prompt (all platforms):

```cpp
gSystem->SetBuildDir("build", true);
gSystem->Setenv("MOLLER_HIT_FACES", "build/surface-faces.tsv");
.L geometry/load_hits.C+
load_hits("input.root", "build/hits.json", "hit.pid == 11");
.L Deconvolution.C+
Deconvolution("inputs.txt", "build/deconvolution.root");
```

`load_hits` above counts selected hits. `Deconvolution` defaults to PE weighting
for main tiles and ShowerMax. `inputs.txt` lists one ROOT path per line.
Use forward slashes in Windows paths at the ROOT prompt.
The app calls these same analyzers and converts results to browser JSON; it does
not require the manual commands. See the macro for optional arguments and cuts.

## Corrections — 2026-11-02

Main-detector PE weighting previously passed kinetic energy (`hit.k`) to the
response table, whose documented input is total energy (`hit.e`). This could
choose the wrong bin or the unit fallback near threshold. Both PE maps and the
deconvolution now use `hit.e`; kinetic-energy cuts are unchanged. Tests use the
actual energy tables at low energies, a bin boundary and high energy. Previously
saved PE results need reanalysis; the effect on production totals is not yet measured.

Finite-MC covariance previously subtracted nearly equal raw moments. Constant
scores could produce a tiny negative variance and reject valid input. The ROOT
scanner now accumulates centered scalar/pair moments, including implicit zero
histories and cross-category covariance. This changes numerical evaluation, not
the estimator; negative or indefinite input covariance is still rejected.
