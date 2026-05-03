# OrbitalAnalyzer Documentation

## Overview
The `OrbitalAnalyzer` class provides a unified interface for orbital analysis, AO mapping, and chemically meaningful classification (core, sigma, pi, lone pair, etc.) in VeloxChem. It is the central entry point for all orbital diagnostics used by the NBO and VB drivers, ensuring consistent and robust orbital labeling and mapping for arbitrary molecules and basis sets.

## Features
- AO-to-atom and AO-to-angular-momentum mapping
- NAO/NPA payload construction when molecular orbitals are provided
- Shared orbital candidate records for core, sigma, pi, lone pair, radical, antibonding, and Rydberg orbitals
- Centralized logic for NBO/VB orbital diagnostics
- Used by both NBO and VB drivers for seamless integration

## Current status: 2026-05-03

`OrbitalAnalyzer` is now the shared orbital-recognition layer for the NBO and VB drivers.

Implemented and validated:

- AO-to-atom and AO-to-angular-momentum maps are provided through one public analyzer entry point.
- RHF/UHF molecular orbitals can be converted into a shared NAO/NPA payload with total-density and spin-density diagnostics.
- The analyzer returns canonical candidate records for `CR`, `LP`, `SOMO`, `BD(sigma)`, `BD(pi)`, `BD*`, and `RY` objects.
- `NboDriver` consumes this payload for candidate generation while retaining Lewis assignment, resonance enumeration, donor-acceptor diagnostics, NRA/NRT fitting, and reports.
- `VbDriver` consumes the same payload for active-space selection, active/inactive/frozen partitioning, and traceable fixed-orbital VB active spaces.
- Notebook and source-level checks confirm that direct analyzer results, NBO results, and VB diagnostics expose consistent candidate labels and atom/bond assignments.

Current boundaries:

- The analyzer identifies and labels orbital candidates; it does not choose a VB wavefunction model.
- The analyzer does not perform Lewis assignment, NRT fitting, BOVB, or orbital optimization.
- `orbitalclassifier.py` remains a private implementation helper. Driver-facing code should depend on `OrbitalAnalyzer`, not on a separate public classifier.

## API

### Initialization
```python
from veloxchem import OrbitalAnalyzer, OrbitalAnalyzerOptions
analyzer = OrbitalAnalyzer(molecule, basis, options=OrbitalAnalyzerOptions(threshold=1e-12))
```

### Running Analysis
```python
diagnostics = analyzer.run()
# diagnostics.ao_to_atom, diagnostics.ao_to_l
```

With molecular orbitals, `run()` returns the full shared payload:

```python
analyzer = OrbitalAnalyzer(molecule, basis, mol_orbs=scf.mol_orbs)
analysis = analyzer.run()
# analysis.nao_data, analysis.mo_analysis, analysis.orbital_candidates
```

### Options
- `threshold`: numerical threshold for AO analysis.
- `include_mo_analysis`: include MO composition in the NAO basis.
- `include_nbo_candidates`: include shared `CR`, `LP`, `BD`, `SOMO`, `BD*`, and `RY` candidate records.
- `mo_analysis_top` and `mo_analysis_threshold`: control the printed/recorded MO-in-NAO composition.
- `lone_pair_min_occupation`, `rydberg_max_occupation`, `bond_min_occupation`, `bond_min_atom_weight`, `pi_min_occupation`, and `conjugated_pi_max_path`: candidate-generation thresholds shared by NBO and VB.

### Diagnostics
- `ao_to_atom`: NumPy array mapping AO index to atom index
- `ao_to_l`: NumPy array mapping AO index to angular momentum quantum number
- `nao_data`: shared NAO/NPA data when molecular orbitals are available
- `spin_data`: spin-resolved NAO density diagnostics
- `mo_analysis`: molecular-orbital composition in the NAO basis
- `orbital_candidates`: shared NBO-like orbital candidates used by NBO and VB

## Usage Example
```python
from veloxchem import OrbitalAnalyzer
analyzer = OrbitalAnalyzer(molecule, basis)
diagnostics = analyzer.run()
print('AO to atom:', diagnostics.ao_to_atom)
print('AO to l:', diagnostics.ao_to_l)
```

## Extension Points
- Add higher-level query helpers for selecting active bonds, lone pairs, radicals, and acceptor partners from the existing candidate records.
- Keep NAO construction, population analysis, and detailed classification (core, valence, bonding, antibonding, lone pair, Rydberg, etc.) centralized here for use by all drivers.
- Do not move VB active-space or wavefunction algebra into the analyzer. The analyzer should identify and label orbital candidates; `VbDriver` decides how those candidates become active spaces, structures, spin-adapted matrix elements, and VB weights.

## Future improvements

### Metal-ligand recognition

The next major analyzer-level chemistry extension is recognition of metal-ligand bonding patterns. This should remain an orbital-recognition layer, not a Lewis/VB decision layer.

Planned capabilities:

1. **Metal-center and ligand-donor detection**
	- Detect transition-metal and main-group coordination centers from element labels, connectivity, and coordination number.
	- Identify ligand donor atoms and donor orbital candidates from `LP`, `SOMO`, pi-donor, and low-occupation acceptor spaces.
	- Preserve hapticity and bridging metadata where the connectivity graph supports it.

2. **Coordination-bond candidate records**
	- Add neutral candidate metadata for metal-ligand sigma donation, pi donation, and pi back-donation.
	- Keep the existing `type`/`subtype`/`atoms`/`occupation`/`coefficients` schema, extending it with optional fields such as `metal_atom`, `ligand_atoms`, `coordination_mode`, and `donor_acceptor_role`.
	- Avoid adding selected Lewis assignments in the analyzer; expose candidates only.

3. **d-orbital and local-frame diagnostics**
	- Report metal-centered d-manifold character in the NAO basis.
	- Add local-frame information suitable for distinguishing sigma, pi, delta, and back-donation channels when the geometry is sufficiently well defined.
	- Keep orientation diagnostics robust under near-degenerate d-orbital rotations.

4. **Downstream contracts**
	- Let `NboDriver` decide whether metal-ligand candidates enter Lewis/NRT alternatives.
	- Let `VbDriver` decide whether metal-ligand candidates become active orbitals in a VB active space.
	- Add tests that verify ligand-field candidate metadata without requiring a final NBO or VB interpretation.

Acceptance checks:

- Simple coordination complexes expose metal atom, donor atom, and candidate orbital metadata deterministically.
- Closed-shell donor-only complexes produce ligand-to-metal donation candidates without forcing them into occupied Lewis assignments.
- Systems with plausible back-donation expose metal d to ligand pi-star diagnostic channels when the NAO density supports them.
- Existing organic `BD(sigma)`, `BD(pi)`, `LP`, `SOMO`, `BD*`, and `RY` candidate labels remain unchanged.

## Alignment
- The `OrbitalAnalyzer` is the single source of truth for orbital analysis in VeloxChem, used by both the NBO and VB drivers.
- Keeps all orbital logic, diagnostics, and mapping consistent and robust across the codebase.
- `orbitalclassifier.py` is not a separate public classifier. It is only a private implementation helper for candidate construction; drivers should consume the `OrbitalAnalyzer` payload.

## Shared NBO/VB Analysis Contract

`OrbitalAnalyzer` returns one structured analysis payload that both `NboDriver` and `VbDriver` consume directly. This avoids separate NBO/VB orbital-recognition paths and makes the analyzer the single place where NAOs, orbital classes, labels, and atom/bond metadata are defined.

### Payload

`OrbitalAnalyzer.run(...)` returns:

- `ao_to_atom` and `ao_to_l` maps
- `nao_data` with NAO transform, density, populations, atom map, and angular-momentum map
- `mo_analysis` in the NAO basis when molecular orbitals are provided
- `orbital_candidates`: canonical records for `CR`, `LP`, `BD(sigma)`, `BD(pi)`, `SOMO`, `BD*`, and `RY`
- stable orbital type/subtype/index records and atom/bond metadata suitable for both NBO reporting and VB active-space construction

### Driver Integration

1. `NboDriver` calls `OrbitalAnalyzer` for AO maps, NAO construction, MO-in-NAO analysis, and NBO candidate generation, while keeping Lewis assignment, resonance enumeration, NRA/NRT fitting, and NBO reports in the NBO driver.
2. `VbDriver` calls the same `OrbitalAnalyzer` payload instead of running its own classifier or indirectly using NBO-specific internals. It translates selected analyzer candidates into `VbOrbital` objects and default VB structures.
3. Shared records should use neutral field names (`type`, `subtype`, `atoms`, `occupation`, `coefficients`, `source`, `label`) so they remain valid for both NBO and VB workflows.

### Existing NAO Payloads

For workflows that already have NAO data, `OrbitalAnalyzer.classify_nao_data(...)` classifies that payload through the same shared candidate builder. This keeps the public VB and NBO interfaces on `OrbitalAnalyzer` even when the candidate implementation lives in a private helper module.

### Communication Check

The current NBO/VB communication layer has been validated by running the analyzer, NBO driver, and VB driver on the same compiled source and checking that the same candidate labels and atom/bond assignments are visible through each path.

Follow-up source-level regression tests should verify that:

- `NboDriver` and `VbDriver` see the same orbital labels and atom/bond assignments
- no VB code imports a separate orbital classifier
- the NBO selected Lewis table and the VB default active-orbital list are derived from the same analyzer candidates

The current VB active-space builder is the first consumer of this contract. For H₂, it selects the analyzer H-H `BD(sigma)` candidate, builds one active bond, generates covalent and ionic structures, and passes those structures to the VB two-electron algebra without reclassifying orbitals. The next analyzer-facing VB extensions should keep the same rule: use analyzer candidates for chemical selection, but keep VB Hamiltonian/overlap algebra inside `VbDriver`.

## VB Sigma/Pi Separation Checkpoint

`VbDriver` now exposes a driver-level view of the analyzer candidates without duplicating analyzer logic. The diagnostics include:

- `candidate_partitions`: groups candidates into core, sigma bonds, pi bonds, lone pairs, radicals, antibonds, Rydberg, and other records.
- `bond_candidate_summary`: a compact table of two-center `BD(sigma)` and `BD(pi)` records.
- `suggested_active_candidates`: the highest-occupation pi candidate followed by the highest-occupation sigma candidate, when present.

This is intended for the ethylene checkpoint: confirm that the C=C pi bond is separated from the C=C and C-H sigma framework before turning that pi candidate into a VB active space. The analyzer remains responsible for identifying the candidates; `VbDriver` remains responsible for deciding which candidate is active, which candidates are frozen/inactive, and how structures and matrix elements are generated.

The first polyatomic consumer of this partition is the ethylene one-active-π prototype. With `VbComputeOptions(active_candidate_label="BD_pi_4", active_candidate_subtype="pi", freeze_inactive_orbitals=True)`, `VbDriver`:

- selects the analyzer C=C `BD(pi)` candidate as the active two-electron/two-orbital subsystem,
- projects that candidate onto the two carbon NAO blocks to form atom-centered π active orbitals,
- leaves the C=C sigma and C-H sigma candidates inactive in diagnostics, and
- embeds the active π Hamiltonian in the frozen RHF reference density, computed as the total HF density minus the selected active pair density.

This keeps sigma recognition in the analyzer and VB wavefunction algebra in `VbDriver`, while providing a chemically traceable route from shared candidates to a frozen-sigma active π model.

## Multi-Center Pi Active-Space Roadmap

The ethylene prototype deliberately activates only one analyzer `BD(pi)` candidate, so it remains a two-electron/two-orbital problem. The next VB layer should use the same analyzer payload to build larger fixed-orbital π active spaces while keeping sigma/core orbitals frozen at the HF reference level.

Target systems and active spaces:

- allyl cation: 2 electrons in 3 π orbitals, singlet;
- allyl radical: 3 electrons in 3 π orbitals, doublet;
- allyl anion: 4 electrons in 3 π orbitals, singlet;
- butadiene: 4 electrons in 4 π orbitals, singlet;
- benzene: 6 electrons in 6 π orbitals, singlet.

For these systems, `OrbitalAnalyzer` should continue to provide recognition only:

- pi-capable atom identification through NAO atom/angular-momentum maps,
- `BD(pi)` candidates and their atom labels,
- `SOMO`/spin diagnostics for open-shell cases,
- `LP` and `CR` candidates for frozen/inactive partitions.

`VbDriver` should own the multi-center interpretation:

- selecting the active π atom set,
- constructing one localized π active orbital per active center,
- assigning active electron count and spin sector,
- generating resonance/ionic VB structures,
- building the active-space Hamiltonian and overlap matrices,
- embedding the active model in the frozen HF sigma/core density.

This separation should remain strict. BOVB and orbital relaxation should not be added to the analyzer; they belong in `VbDriver` as future work on top of the validated fixed-orbital multi-center π active-space engine.
