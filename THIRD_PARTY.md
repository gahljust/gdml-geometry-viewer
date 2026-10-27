# Response origins

- Main evaluator: remoll develop `15eb235236b998ccb64927ac0f1df485a5f2142c` (PRs 649/647).
- Main CSVs: MOLLER DocDB **Moller-doc-1662-v1**, `MainDet_ResponseFunction.zip`.
  Table hashes are in `data/response-manifest.json`.
- ShowerMax CSVs: remoll `analysis/showermax/response_lookup/data`;
  bundled adapter chooses the module by detector ID.
- `root/metadata_schema.hh`: read-only adaptation of remoll run metadata.

ROOT, Geant4 and third-party sources retain their upstream licensing.
Response models are not measured detector calibrations.
