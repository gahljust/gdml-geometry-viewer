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

For each event, sum selected signals per detector, then accumulate
`R = Σ w q` and `B = Σ w q qᵀ` (one rate weight, not its square).
`Crs = Brs / (P² T Rr Rs)`, `M = Hᵀ C⁻¹ H`; statistical covariance is `M⁻¹`.
Finite simulation uncertainty propagates through H separately. Fit metadata checks
are stricter than hit viewing; missing information must not become zero.

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
