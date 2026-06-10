# PR Summary

This PR adds `atomic_property` metadata to atom-resolved charge datasets written to HDF5 and stores fitted ESP/RESP charges in dedicated result groups instead of folding them into `scf` results.

## What changed

- `resultsio.py`
  - adds shared HDF5 metadata support for atom-resolved datasets
  - tags `nuclear_charges`, `esp/esp_charges`, and `resp/resp_charges` with `atomic_property`
- `EspChargesDriver`
  - writes fitted charges to a top-level `esp` group as `esp_charges`
  - keeps the implementation independent from `scf_results` rewrites
- `RespChargesDriver`
  - writes fitted charges to a top-level `resp` group as `resp_charges`
  - reuses the same simplified HDF5 write path as ESP
- tests
  - verify metadata tagging for charge datasets
  - verify ESP/RESP HDF5 persistence for both API-style runs and input-file runs where the driver performs the SCF step internally

## Resulting HDF5 layout

- `scf/...` remains reserved for SCF data
- `esp/esp_charges` stores ESP-fitted charges with `atomic_property = "ESP Charges"`
- `resp/resp_charges` stores RESP-fitted charges with `atomic_property = "RESP Charges"`

## Validation

- `python -m pytest tests/test_resultsio.py::test_atomic_property_metadata_for_charge_result_datasets tests/test_esp_charges.py tests/test_resp_charges.py -q`
  - `9 passed, 1 warning`
- direct CLI check with `python -m veloxchem tests/data/water_esp.inp`
  - confirmed top-level `esp` group
  - confirmed `esp/esp_charges`
  - confirmed `atomic_property: ESP Charges`