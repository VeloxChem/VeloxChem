# NboDriver implementation summary

This document summarizes the implemented `NboDriver` architecture in VeloxChem and the remaining implementation plan. The core NBO/NRA/NRT roadmap is now in place. The next work is focused on chemical interpretation, resonance-class grouping, acceptor reporting, and broader validation.

## Architecture

The driver is organized as a layered analysis pipeline:

```text
AO/SCF density and overlap
-> NAO/NPA transformation and diagnostics
-> MO composition in the NAO basis
-> NBO candidate generation
-> Lewis/resonance assignment
-> NRA/NRT density fitting
-> reports and structured diagnostics
```

NRA/NRT remains a post-processing density-fit layer on top of generated Lewis/resonance alternatives. It does not feed back into candidate generation, primary Lewis assignment, or score/ranking weights.

## Implemented capabilities

- Conservative NAO/NPA construction with explicit foundation invariants: `T^T S T`, `D_NAO = T^T S P S T`, electron count, density symmetry, and charge conservation.
- Deterministic handling of near-degenerate atom-block NAO rotations.
- MO-in-NAO composition analysis for restricted and unrestricted orbitals.
- NBO candidate generation for `CR`, `LP`, `BD(sigma)`, `BD(pi)`, `SOMO`, radical one-electron `LP`, radical one-electron `BD(pi)`, `BD*`, and `RY`.
- Lewis/NBO primary assignment with explicit candidate electron counts.
- Bond and pi-bond constraints through `NboConstraints` or dictionary input.
- Lewis/resonance alternatives with separate `sigma_nbo_list` and `pi_nbo_list` fields.
- Fixed/active partition fields for resonance analysis: `fixed_nbo_list`, `active_nbo_list`, `active_pi_nbo_list`, `active_lone_pair_nbo_list`, and `active_one_electron_nbo_list`.
- Lewis electron accounting, atom valence electron counts, formal-charge diagnostics, octet diagnostics, and transparent `score_terms`.
- Lone-pair donation in the resonance active space.
- Closed-shell NRA/NRT density fitting with `selected`, `pi`, `valence`, and `full` subspaces.
- `nra_report(level="summary")` and `nra_report(level="full")` report formatters.
- Molecule-independent pi-bond structure signatures for NRA/NRT reports and prior matching.
- Open-shell spin-resolved NAO densities and one-electron radical alternatives.
- Open-shell NRA/NRT fitting through a combined total-density and spin-density residual.
- `BD*` antibonding complements as same-subspace orthogonal partners to occupied `BD` candidates.
- One-center `RY` acceptor complements as candidate-only records.
- Donor-acceptor diagnostics from occupied donors to candidate-only `BD*` and `RY` acceptors using a density-coupling diagnostic.
- Optional prior-regularized NRA/NRT weights with explicit prior metadata.
- Regression coverage for foundation invariants, NBO counts, constraints, sigma/pi partitions, Lewis accounting, score terms, lone-pair donation, structure-pool invariants, open-shell SOMO/radical alternatives, antibonding/Rydberg complements, donor-acceptor diagnostics, NRA/NRT weights, NRA/NRT priors, open-shell spin-resolved NRA/NRT, labels, and reports.

## Current status: 2026-05-03

Recent implementation work added the general resonance machinery needed for aromatic, radical, anionic, and polar pi systems without molecule-specific rules.

- Pure pi resonance alternatives are generated through bounded exact pi matchings instead of all-subset enumeration. This removes the acene scaling cliff while preserving exact non-overlapping pi-bond patterns.
- Generic polar pi resonance units are handled through coupled pi/lone-pair assignments. Nitro-like fragments are represented as formal charge-separated alternatives with a positive center and a negative terminal atom, but the implementation is based on general polar X-Y-X connectivity and electronegativity patterns rather than named functional groups.
- Alternatives carry structured metadata for `pi_bonds`, `active_lone_pair_atoms`, `active_positive_atoms`, and `active_one_electron_atoms`.
- Alternatives are sorted by descending score/ranking weight before final ranks are assigned.
- `show_structures(...)` visualizes pi bonds, radical centers, lone-pair/negative centers, and positive centers.
- The compact smoke-test notebook contains allyl, benzene, acene, coronene, and 1,3,5-trinitrobenzene examples. The trinitrobenzene example is intended to run after geometry optimization and should continue to expect 16 simple formal resonance structures.
- `NboDriver` now obtains AO maps, NAO/NPA data, spin data, MO-in-NAO diagnostics, and candidate records from the shared `OrbitalAnalyzer` payload. The NBO driver owns the Lewis assignment, resonance alternatives, donor-acceptor diagnostics, NRA/NRT fitting, and reports.
- The current donor-acceptor layer reports density-coupling diagnostics to candidate-only `BD*` and `RY` acceptors. It is not yet a second-order perturbation-energy layer.
- NRA/NRT weights remain density-reconstruction weights and are intentionally separated from VB/wavefunction weights owned by `VbDriver`.

The main remaining gap is not the enumeration of individual structures. The next missing layer is grouping chemically equivalent structures into resonance classes so that degeneracies and class-level weights can be reported separately from raw localized-orbital weights.

## Weight terminology

The implementation keeps three weight concepts separate:

| Name | Meaning |
| --- | --- |
| Score/ranking weight | Softmax weight derived from the Lewis assignment score. |
| NRA/NRT density weight | Nonnegative density-fit weight from `results["nra"]`. |
| VB/wavefunction weight | State-mixing quantity owned by `VbDriver`; initially implemented for the H₂ two-orbital spin-adapted active-bond model, not part of NBO/NRA/NRT weighting. |

User priors are not a fourth physical weight. They are optional guidance used only in the regularized NRA/NRT density fit and are reported separately as prior metadata.

## Current scientific scope

The driver provides a transparent NBO/NRA/NRT analysis layer suitable for method development and interpretation. Several quantities are deliberately reported as diagnostics rather than final physical observables:

- Donor-acceptor diagnostics use the NAO density matrix as a coupling operator. They are not second-order perturbation energies because no NBO Fock matrix or energy denominator is used.
- NRA/NRT weights are density-reconstruction weights. They are not VB amplitudes or state-mixing coefficients.
- Pi-bond signatures are report annotations. The internal data model remains based on sigma/pi and fixed/active partitions, not named structure classes.
- NBO candidates are density-block natural orbitals in NAO space; a separate NHO directionality layer is outside the current implementation.

## Implementation plan

### 1. Resonance signatures and class metadata

Goal: give every Lewis/resonance alternative a stable, chemically meaningful identity that can be compared, grouped, tested, and reported.

Implementation steps:

1. Add a canonical resonance signature builder for each alternative using only general fields:
	- sorted `pi_bonds`
	- sorted `active_lone_pair_atoms`
	- sorted `active_positive_atoms`
	- sorted `active_one_electron_atoms`
	- optional sigma-bond changes if sigma resonance is later enabled
2. Store this signature on every alternative as structured data, not only as a display string.
3. Add a compact human-readable label derived from the signature for reports and debugging.
4. Keep the current individual `rank` and `weight` behavior unchanged so existing output remains stable.

Acceptance checks:

- Benzene alternatives have two distinct signatures.
- Allyl cation, radical, and anion alternatives differ by the correct pi/one-electron/lone-pair metadata.
- 1,3,5-trinitrobenzene has 16 distinct individual signatures.

### 2. Symmetry-equivalent resonance classes

Goal: group alternatives that are chemically equivalent under atom relabeling and report degeneracies separately from individual localized weights.

Implementation steps:

1. Build a molecular graph from atom labels and connectivity already available to the driver.
2. Generate graph automorphisms that preserve element labels and bond connectivity. Start with a conservative implementation suitable for small and medium organic resonance systems.
3. Apply each automorphism to the canonical resonance signature and choose the minimum transformed signature as the resonance-class key.
4. Add `resonance_class_id`, `class_signature`, and `class_degeneracy` to each alternative.
5. Add `results["resonance_classes"]` as a list of class records containing member ranks, degeneracy, representative rank, and class-level weights.
6. Keep this layer optional or lightweight if automorphism generation is too expensive for large systems.

Acceptance checks:

- Benzene reports one Kekule class with degeneracy 2 when symmetry grouping is enabled.
- 1,3,5-trinitrobenzene groups the 16 individual structures into symmetry classes consistent with two ring Kekule choices and equivalent nitro substitutions.
- Raw individual weights remain visible, while class summaries show the grouped interpretation.

### 3. Class-level weights and reporting

Goal: make output chemically readable when multiple equivalent structures have slightly different raw localized scores.

Implementation steps:

1. Compute class weight as the sum of member score/ranking weights.
2. Compute optional class mean and spread for diagnostics:
	- `class_weight_sum`
	- `class_weight_mean`
	- `class_weight_min`
	- `class_weight_max`
3. Add a compact resonance-class table to the NBO report above the individual structure table when classes are present.
4. Preserve the individual structure table sorted by descending individual weight.
5. Update `show_structures(...)` titles to include both individual rank and class id when class metadata exists.

Acceptance checks:

- Trinitrobenzene output explains why individual weights are not exactly equal while still showing grouped class degeneracies.
- Existing reports remain readable when no class grouping was requested or produced.

### 4. BD*/RY acceptor-candidate report

Goal: make the implemented acceptor space visible without mixing acceptor candidates into the selected occupied Lewis table.

Implementation steps:

1. Add a separate optional report section for candidate-only acceptors.
2. Include top `BD*` and `RY` candidates by occupation and donor-acceptor diagnostic strength.
3. Show donor, acceptor, acceptor type, atom/bond label, occupation, and density-coupling diagnostic.
4. Keep the section under `nbo_report(level="full")` first; consider a separate option later if output becomes too verbose.
5. Make clear in labels that these are acceptor candidates, not selected occupied Lewis NBOs.

Acceptance checks:

- Existing selected Lewis tables do not show `BD*` or `RY` as occupied assignments.
- Full reports expose the strongest `BD*` and `RY` acceptors in a compact table.
- Donor-acceptor diagnostics already covered by tests remain consistent.

### 5. Regression tests and notebook checks

Goal: lock down the general behavior added during this work before expanding NRT validation.

Implementation steps:

1. Add a focused test for pure pi matching:
	- benzene returns the two Kekule alternatives under allowed-pi constraints
	- acene-like systems avoid combinatorial explosion
2. Add a focused test for polar resonance:
	- a nitro-containing system returns the expected coupled positive/negative-center metadata
	- electron counts remain exact
3. Add a weight-order test:
	- alternatives are monotonically sorted by descending score/ranking weight after final rank assignment
4. Add a visualization metadata smoke test:
	- alternatives contain the fields needed by `show_structures(...)` for pi, lone-pair/negative, positive, and radical markers
5. Keep the notebook as a human-facing smoke test rather than the only validation path.

Acceptance checks:

- Tests fail if polar resonance falls back to one structure.
- Tests fail if final ranks no longer follow descending weight.
- Tests fail if active metadata needed for structure visualization disappears.

### 6. Later NRT/NRA validation

Goal: validate the density-fitting layer after the resonance-structure layer is stable.

Implementation steps:

1. Re-run closed-shell NRA/NRT tests on benzene, allyl ions, and nitro systems.
2. Re-run open-shell NRA/NRT tests on allyl radical and related radicals.
3. Compare selected, pi, valence, and full subspace fits for representative systems.
4. Check sensitivity to priors and regularization parameters.
5. Add publication-style examples only after the class grouping and acceptor report are stable.

Acceptance checks:

- NRA/NRT weights remain clearly reported as density-reconstruction weights, not VB amplitudes.
- Priors remain explicit metadata and do not silently change Lewis enumeration.
- Open-shell total-density and spin-density residuals remain controlled.

### 7. Longer-term method extensions

These remain useful, but they are lower priority than the interpretation and validation work above.

- Add a true NHO/hybrid-direction layer on top of the current NAO-space candidate vectors.
- Add an NBO Fock-matrix layer and second-order donor-acceptor perturbation energies.
- Expand optional user-facing structure annotations without adding molecule-specific logic to the core driver.
- Add broader validation notebooks and publication examples for additional radicals, zwitterionic systems, amides, nitro compounds, and transition-metal fragments.
- Investigate VB/BOND-style state coupling if wavefunction weights are desired.

### 8. Metal-ligand and coordination NBO analysis

Goal: consume `OrbitalAnalyzer` metal-ligand candidate records and make coordination chemistry interpretable without hard-coding organic Lewis assumptions into the analyzer.

Implementation steps:

1. Accept analyzer candidate metadata for metal-ligand sigma donation, pi donation, and metal-to-ligand back-donation. **Initial diagnostic exposure is implemented through `metal_ligand_diagnostics`.**
2. Add coordination-aware electron-accounting terms that can coexist with the current organic duet/octet diagnostics but do not force octet-style penalties on transition-metal centers.
3. Add scoring terms for ligand donor occupation, metal acceptor/back-donation occupation, and charge balance.
4. Keep coordination candidates separate from selected occupied organic NBO tables until their electron-count semantics are explicit.
5. Add report sections for metal center, ligand donor atoms, candidate donor/acceptor channels, d-manifold character, and strongest donation/back-donation diagnostics.
6. Validate first on small closed-shell donor complexes, then on pi-acceptor ligands and open-shell metal fragments.

Acceptance checks:

- Organic NBO output remains unchanged when no metal center is present.
- Simple metal-ligand systems expose coordination candidate records and donor/acceptor diagnostics without crashing Lewis assignment.
- Octet diagnostics remain meaningful for ligand atoms while metal-center diagnostics are reported separately.
- Donation and back-donation channels are visibly distinct in full reports.

### 9. Unified Orbital Analysis and Classification

Goal: Centralize all orbital analysis, NAO construction, and classification in a single `OrbitalAnalyzer` class/module. This ensures that both the NBO and VB drivers use the same, chemically meaningful set of orbitals and labels, eliminating code duplication and inconsistencies.

Current status: `NboDriver` and `VbDriver` now consume the same `OrbitalAnalyzer` payload (`nao_data`, `spin_data`, `mo_analysis`, and `orbital_candidates`). The communication layer has been confirmed: direct analyzer runs, NBO results, and VB diagnostics expose matching candidate records. The VB driver now uses those records to build the first generated H₂ one-active-bond model, while the NBO driver remains responsible for Lewis assignment, resonance alternatives, donor-acceptor diagnostics, NRA/NRT fitting, and reports. Candidate classification is now analyzer-owned implementation detail inside `orbitalanalyzerdriver.py`, not a second driver-facing classifier module. The analyzer also now emits metal-ligand `ML/sigma-acceptor` and `ML/pi-donor` diagnostic records, separating ligand-to-metal sigma donation from metal-to-ligand pi back-donation without forcing either channel into the primary Lewis assignment.

Metal-ligand update: `NboDriver.compute()` exposes analyzer `ML` records through `metal_ligand_diagnostics`, while the primary Lewis assignment remains organic/electron-counted and excludes `ML` records. `VbDriver.compute()` can explicitly activate these records for fixed-orbital sigma-only and sigma-plus-back-donation VB-CI diagnostics; this remains separate from NBO Lewis selection.

Implementation steps:

1. Refactor the current orbital classification logic into a standalone `OrbitalAnalyzer` class. **Completed for the shared payload and candidate records.**
2. Keep NBO-specific Lewis assignment and resonance machinery in `NboDriver`. **Completed architecturally.**
3. Use the analyzer payload as the input for VB active-space and default-structure construction. **Completed for the current fixed-orbital H₂, ethylene, and π-ladder checkpoints; future VB work should keep the same analyzer contract.**
2. `OrbitalAnalyzer` should:
   - Accept a molecule, basis, and (optionally) orbitals.
   - Run SCF if orbitals are not provided.
   - Compute NAOs (via NBO driver) if requested.
   - Classify all orbitals (core, valence, bonding, antibonding, lone pair, Rydberg, etc.).
   - Store canonical MOs, NAOs, classification results, and mapping between representations.
3. Update the NBO driver to accept orbitals and classification from `OrbitalAnalyzer`.
4. Ensure all NBO diagnostics and reports reference the unified classification and labels.
5. Document the new workflow and update all relevant examples and notebooks.

Immediate regression target:

- Run `OrbitalAnalyzer` once and verify that NBO candidate labels and atom/bond assignments are exactly the records used by the NBO primary assignment.
- Verify that VB diagnostics report the same candidate labels for the same molecule.
- Verify that no driver imports or calls a public `OrbitalClassifier` API.

Benefits:
- Guarantees that NBO analysis is always performed on a well-defined, classified set of orbitals.
- Enables seamless integration with the VB driver and other modules.
- Simplifies diagnostics and user interpretation.

## Tomorrow restart notes

Start from the current NBO scope and keep coordination chemistry diagnostic until its electron-count semantics are explicit:

1. Re-open `docs/metal_ligand_recognition.ipynb` and inspect the `NboDriver.compute()` payload for `metal_ligand_diagnostics` on Pd--NH3 and Pd--PH3.
2. Keep `ML/sigma-acceptor` and `ML/pi-donor` out of the primary selected Lewis table for now. They are candidate diagnostics, not occupied organic Lewis NBOs.
3. The next NBO implementation work should focus on resonance-class grouping, class-level reporting, and the separate `BD*`/`RY` acceptor-candidate report before production coordination-NBO logic.
4. For coordination chemistry, the next safe step is fuller reporting of metal center, ligand donor atoms, d-manifold character, and donation/back-donation diagnostics without changing the current Lewis assignment.
5. Longer-term coordination NBO/NRT work should add explicit metal-center electron-accounting rules and, later, true NBO Fock-matrix second-order donor-acceptor energies. The current donor-acceptor values are density-coupling diagnostics only.

Current safe stopping point: `NboDriver` consumes the shared analyzer payload, exposes metal-ligand diagnostics, preserves organic Lewis/NRA/NRT behavior, and leaves wavefunction weights to `VbDriver`.

---
