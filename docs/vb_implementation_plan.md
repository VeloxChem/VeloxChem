# Valence Bond (VB) Driver Implementation Plan for VeloxChem

## Overview
This plan outlines a staged approach to implementing a general Valence Bond (VB) driver in VeloxChem, starting from minimal VB-SCF for H2 and building up to BOVB and user-friendly wavefunction specification. The approach is designed to fit VeloxChem's architecture and leverage existing NBO/NAO machinery where appropriate, but to keep VB as a true wavefunction method with non-orthogonal structure algebra.

---

## Phase 0: Driver Skeleton

### Tasks (Completed)
- [x] Create a new file: `src/pymodule/vbdriver.py`.
- [x] Define the following dataclasses:
  - `VbOrbital`: label, coefficients, center, kind.
  - `VbStructure`: label, occupation, spin, charge_pattern.
  - `VbComputeOptions`: mode, max_iter, conv_thresh, optimize_orbitals, include_ionic, include_bovb.
- [x] Implement a `VbDriver` class with a `compute` method matching the API:
  - Accepts: molecule, basis, structures, orbitals, reference_orbitals, options.
  - Returns a dictionary with keys: energy, structure_coefficients, overlap, Hamiltonian, weights, orbitals, diagnostics.
- [x] Add docstrings and placeholder logic for each class/method.

Phase 0 is complete. Proceeding to Phase 1.

### Example API
```python
from veloxchem.vbdriver import VbDriver, VbComputeOptions, VbOrbital, VbStructure

vb = VbDriver()
results = vb.compute(
    molecule, basis,
    structures=[...],
    orbitals=[...],
    reference_orbitals=scf.mol_orbs,
    options=VbComputeOptions(),
)
```

### Result Schema
The result dictionary should include:
- `energy`: total VB energy
- `structure_coefficients`: VB structure coefficients
- `overlap`: structure overlap matrix
- `Hamiltonian`: structure Hamiltonian matrix
- `weights`: structure weights (various definitions)
- `orbitals`: final orbital coefficients
- `diagnostics`: convergence and algebra diagnostics

### Notes
- The initial implementation should use placeholder logic and return mock values for all outputs.
- Add a docstring to each class and method describing its intended role.
---

## Phase 1: Minimal H2 VB-CI


### Tasks

- [x] Implement fixed-orbital VB-CI for H2 in minimal basis
  - [x] Add minimal working example for H2 (two 1s-like orbitals, two electrons)
  - [x] Build covalent and ionic structures explicitly (now three: covalent, ionic1, ionic2)
  - [x] Expand each structure into determinants (occupation patterns set)
  - [x] Compute all overlap and Hamiltonian elements from AO integrals
  - [x] Solve generalized eigenvalue problem `H c = E S c`
  - [x] Return all required results in the schema (energy, coefficients, overlap, Hamiltonian, weights, orbitals, diagnostics)
- [x] Validate against FCI and check structure weights, dissociation behavior, and matrix properties

#### Progress
- [x] Code scaffolding for VB-CI started in `vbdriver.py`
- [x] Minimal H2 test/example added (with three structures)
- [x] AO-integral-based overlap and Hamiltonian matrix construction implemented
- [x] FCI reference calculation integrated in test notebook (using AO integrals from VeloxChem)
- [x] Results compared: VBDriver placeholder returns zeros, FCI returns correct energy and coefficients (see notebook for output)

**Next:**
- Update VBDriver to return real AO-integral-based results (not placeholder)
- Ensure correct structure algebra and weights
- Analyze dissociation behavior and matrix properties

## Current Status (April 2026)

- The VBDriver skeleton and dataclasses are implemented and documented.
- Minimal H₂ VB-CI logic is complete: the driver computes overlap and Hamiltonian matrices from AO integrals, solves the generalized eigenproblem, and returns energies, coefficients, and diagnostics.
- The AO ERI C++/pybind11 binding is available and integrated, enabling general AO integral access from Python.
- The test notebook demonstrates and validates the minimal H₂ VB-CI implementation, with automated comparison to FCI (energies, coefficients, weights).
- The driver output matches FCI for the minimal case, confirming correct algebra and integration.
- Structure weights and dissociation behavior can now be analyzed and visualized.

## Current Status (April 2026)

- The VBDriver skeleton, dataclasses, and API are implemented and documented.
- Minimal H₂ VB-CI logic is complete: the driver computes overlap and Hamiltonian matrices from AO integrals, solves the generalized eigenproblem, and returns energies, coefficients, and diagnostics.
- AO ERI C++/pybind11 binding is available and integrated, enabling general AO integral access from Python.
- The test notebook demonstrates and validates the minimal H₂ VB-CI implementation, with automated comparison to FCI (energies, coefficients, weights).
- The driver output matches FCI for the minimal case, confirming correct algebra and integration.
- Robust AO-to-atom mapping is implemented for all geometries and basis sets.
- Dissociation scan logic and VB-SCF orbital optimization are implemented and run without error.
- **However, the dissociation scan does not yet show a physical minimum, and structure weights do not change as expected.**
    - This indicates a remaining bug in the structure algebra or orbital optimization logic (likely in the AO-to-atom mapping, orbital construction, or the optimization loop).
- Diagnostics and debugging are ongoing to resolve this issue.

## Next Steps (In Progress)

- [ ] Debug and fix structure algebra and orbital optimization so that dissociation scans yield a physical minimum and nontrivial structure weights.
- [ ] Generalize VBDriver to arbitrary diatomics and polyatomics
  - [ ] Support arbitrary numbers of atoms, electrons, and orbitals
  - [ ] Allow user-defined VB structures and occupations
  - [ ] Implement robust overlap and Hamiltonian construction for general cases
- [ ] Expand diagnostics and analysis
  - [ ] Compute and report structure weights, natural orbitals, charge/spin analysis
  - [ ] Provide tools for dissociation curves and structure mixing analysis
- [ ] Extend testing and validation
  - [ ] Add tests for more molecules and edge cases
  - [ ] Compare with literature and other codes for benchmark systems
- [ ] Improve documentation and examples
  - [ ] Document API and provide example notebooks
  - [ ] Add guidance for customizing VB structures and interpreting results

---

*Work is now moving to generalization and orbital optimization. This will enable the driver to handle arbitrary molecules and more flexible VB models, as well as provide deeper analysis and robust validation.*

---

## Phase 2: H2 VB-SCF (Complete)
- Add orbital optimization for H2.
- Use a simple parametrization (e.g., mixing angle between AO_H1 and AO_H2).
- Alternate or jointly optimize structure coefficients and orbital parameters.
- Use finite-difference or simple optimizer for orbital parameters.
- Validate energy, coefficients, and gradients.
- [x] Expose four-index AO electron repulsion integrals to Python via pybind11 binding (C++ layer)

---

## Phase 3: H2 Dissociation Curve (In Progress)
- Automate H2 dissociation scan with VB-SCF.
- Report total energy, covalent/ionic weights, structure overlap condition number, and compare to RHF/UHF/FCI.
- Ensure correct dissociation limit and stable structure algebra.

---

## Phase 4: Generalization (Planned)
- Generalize to arbitrary molecules, electrons, and user-defined structures.
- Robust overlap/Hamiltonian construction for general cases.
- Add tests/examples for diatomics and polyatomics.

---

## Phase 5+: (see original plan for details)
- General two-center bond VB, non-orthogonal engine, user-facing language, BOVB, NBO integration, etc.

---

*Phase 3 (H₂ dissociation curve) is now in progress. All logic remains in the driver; the notebook is for testing and validation only.*

---

## Phase 4: General Two-Center Bond VB
- Generalize to one active bond in larger molecules (e.g., LiH, HF, CH4, ethane).
- Support frozen inactive orbitals plus active VB orbitals.
- Allow user to specify active bonds and structure patterns.
- Validate with chemically meaningful test cases.

---

## Phase 5: General Non-Orthogonal VB Engine
- Implement general non-orthogonal Slater determinant formulas for overlap and Hamiltonian.
- Support arbitrary number of structures and orbitals.
- Add spin-adapted structures, linear dependence detection, and robust eigenvalue handling.
- Provide multiple structure weight definitions (Chirgwin-Coulson, Löwdin, etc.).

---

## Phase 6: User-Facing Wavefunction Language
- Design user input language for VB structures:
  - Level 1: canned constructors (e.g., `VbStructures.h2(include_ionic=True)`)
  - Level 2: chemical objects (bonds, lone pairs, radicals)
  - Level 3: explicit orbital occupation
- Internally normalize to orbital labels and occupation patterns.

---

## Phase 7: BOVB (Breathing-Orbital VB)
- Allow structure-specific orbital sets.
- Optimize orbitals for each structure.
- Compare VB-SCF and BOVB dissociation.
- Add constraints and diagnostics for breathing/polarization.

---

## Phase 8: Integration With NBO
- Use NBO for generating initial localized orbitals and structure suggestions.
- Allow user to convert NBO results to VB input.
- Provide analysis and reporting tools leveraging NBO/NAO machinery.

---

## General Principles
- Keep the algebra explicit and transparent.
- Prioritize correctness and stability of non-orthogonal structure machinery.
- Make all intermediate results (overlap, Hamiltonian, weights, diagnostics) available for inspection and testing.
- Grow the user interface only after the core engine is validated.

---

This plan is designed for incremental, testable progress, with each phase building a foundation for the next. The initial focus is on H2 for maximal transparency and validation, then generalization to more complex systems and user interfaces.