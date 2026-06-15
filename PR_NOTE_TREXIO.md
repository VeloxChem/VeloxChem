# TREXIO PR Note

## Summary

This branch adds a dedicated TREXIO reader/writer module for VeloxChem in `src/pymodule/trexio.py` and wires it into the Python package API.

The implementation now supports:

- writing molecules, Gaussian basis metadata, AO overlap, molecular orbitals, state energy, and ECP data to TREXIO,
- reading molecules, bases, orbitals, and state energy back from TREXIO,
- AO reordering between VeloxChem AO order and TREXIO canonical spherical order,
- basis metadata export aimed at external TREXIO consumers,
- focused regression coverage in `tests/test_trexio.py`, including restricted, unrestricted, and restricted open-shell orbital roundtrips,
- a development notebook for HDF5/TREXIO validation in `notebooks/trexio_example.ipynb`, including larger ROSCF and USCF cases.

## Main Changes

- Added `write_trexio`, `read_trexio`, `read_molecule`, `read_molecule_and_basis`, and `read_molecular_orbitals`.
- Added basis-shell and AO-permutation helpers to map VeloxChem AO ordering to TREXIO AO ordering.
- Exported primitive normalization factors separately via `basis_prim_factor` and raw contraction-like coefficients via `basis_coefficient`.
- Added ECP serialization/deserialization support.
- Added tests for molecule/basis roundtrip, restricted SCF roundtrip, unrestricted orbital roundtrip, restricted open-shell orbital roundtrip, and AO basis metadata layout.
- Added TREXIO ECP regression tests for basis roundtrip and raw ECP metadata roundtrip.
- Simplified the notebook workflow to keep only native HDF5 and TREXIO HDF5 outputs.
- Added heavier notebook regression cases for TEMPO ROSCF/6-31+G* and trityl USCF/cc-pVDZ.

## Validation

- Reviewed branch delta against `master` for `src/pymodule/trexio.py`, `tests/test_trexio.py`, and the TREXIO notebook.
- Ran `conda run -n vlxsrc python -m pytest -q tests/test_trexio.py`: 9 passed.
- Re-ran the notebook summary cell after removing text-backend output; it passed.
- Ran the appended notebook TREXIO roundtrip checks for TEMPO ROSCF/6-31+G* and trityl USCF/cc-pVDZ; both passed.

## Open Follow-Ups Before Merge

No blocking follow-ups identified from the current TREXIO test and notebook validation scope.

One remaining design boundary to keep in mind: TREXIO primitive normalization is implemented through `l=6`, consistent with the underlying VeloxChem basis-function code. For `l > 6`, TREXIO export now fails explicitly instead of writing incorrect metadata.

One separate caveat remains outside the current TREXIO basis/ECP validation scope: for ECP-bearing molecules, molecule-level charge reconstruction on plain `read_molecule()` needs separate clarification in VeloxChem itself. The TREXIO ECP tests added here validate basis/ECP serialization and raw ECP metadata, but do not assert molecular charge equality for the Au/H ECP fixture.

## Suggested Merge Message

Add TREXIO read/write support for molecules, bases, orbitals, and SCF metadata