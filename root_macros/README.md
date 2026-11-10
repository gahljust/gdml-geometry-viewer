# ROOT deconvolution

Use raw remoll Møller, ep-elastic and ep-inelastic ROOT files with the required quartz-tile and ShowerMax-plane hits and saved `run_data` metadata.

## Select files

To select the files you want to analyze, edit `fnms` near the top of [deconvolution_correlated.C](deconvolution_correlated.C). List individual files:

```cpp
std::vector<std::string> fnms = {
    "/path/to/moller_01.root",
    "/path/to/moller_02.root",
    "/path/to/ep_elastic.root",
    "/path/to/ep_inelastic.root"
};
```

Or select multiple files with wildcards:

```cpp
std::vector<std::string> fnms = {
    "/path/to/moller_*.root",
    "/path/to/ep_elastic_*.root",
    "/path/to/ep_inelastic_*.root"
};
```

Relative paths start from the directory where you launch ROOT. Avoid overlapping patterns: duplicate files are rejected. The saved metadata identifies each generator; filenames do not determine it.

## Run

From this directory, choose one command. Replace `/path/to/tables` with your PE response folder.

Main crossing rates + ShowerMax PE (default weighting):

```sh
root -l -b -q 'deconvolution_correlated.C+(false,"/path/to/tables","mixed.root")'
```

PE weighting for both detectors:

```sh
root -l -b -q 'deconvolution_correlated.C+(true,"/path/to/tables","pe.root")'
```

Crossing rates for both, without PE tables:

```sh
root -l -b -q 'deconvolution_correlated.C+(false,"","rates.root",344,0.8,65,1,false)'
```

## Argument positions

The slots inside the parentheses are:

```text
                            1       2       3           4     5    6        7     8
deconvolution_correlated.C+(mainPE, tables, outputFile, days, pol, current, minK, showerMaxPE)
```

| Slot | What to enter | Default if omitted |
|---|---|---|
| **1 — Main PE** | `true` for main-detector PE weighting; `false` for selected crossing counts. | `false` |
| **2 — Tables** | PE response folder in double quotes, e.g. `"/path/to/tables"`. `""` uses automatic discovery when tables are needed. | `""` |
| **3 — Output** | Output ROOT filename in double quotes, e.g. `"result.root"`. It must not already exist. | `"deconvolution_correlated.root"` |
| **4 — Days** | Beam exposure in days. | `344` |
| **5 — Polarization** | Fraction, e.g. `0.8` means 80%, not `80`. | `0.8` |
| **6 — Current** | Beam current in µA, e.g. `65`. | `65` |
| **7 — Minimum energy** | Minimum crossing kinetic energy in MeV. Main PE mode requires this setting to remain `1`. | `1` |
| **8 — ShowerMax PE** | `true` for ShowerMax PE weighting; `false` for selected charged crossing counts. `false` requires slot 1 to be `false` too. | `true` |

You can omit trailing slots, but cannot skip a slot in the middle. For example, to change slot 8, supply slots 1–7 as well. Input ROOT files go in `fnms`, not in these slots.

Main PE mode requires ShowerMax PE enabled; the response model applies its own additional particle, direction and energy selections. Turning PE off does not remove detectors from the fit.

Compatible files are combined by event count within each interaction. Results are written to the named ROOT file.
