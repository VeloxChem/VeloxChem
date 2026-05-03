# Valence Bond Driver Implementation Plan

## Purpose

This document tracks the current implementation strategy for the VeloxChem valence bond (VB) driver. The driver is being built incrementally from small, inspectable fixed-orbital models toward general non-orthogonal VB, BOVB, VBB, and eventually compact user-facing VB specifications.

The guiding principle is:

> Validate chemically transparent fixed-orbital active spaces before adding orbital optimization or breathing-orbital models.

BOVB is still an important target, but it should be introduced only after the active/inactive orbital partition, fixed-orbital Hamiltonian, spin algebra, and diagnostics are reliable for more than a two-orbital test case.

## Current architecture

### Shared orbital analysis

`OrbitalAnalyzer` is the shared owner of orbital recognition. It provides:

- AO-to-atom and AO-to-angular-momentum maps.
- NAO/NPA data from an RHF/UHF reference.
- Shared candidate records for `CR`, `BD(sigma)`, `BD(pi)`, `LP`, `SOMO`, `BD*`, and `RY` candidates.
- Metal-ligand `ML` diagnostic records with separate `sigma-acceptor` and `pi-donor` channels.
- Candidate labels and atom assignments consumed consistently by NBO and VB.

Candidate classification is implemented inside `orbitalanalyzerdriver.py` as analyzer-private helper logic. Public workflows should not depend on a separate classifier module.

### VB driver responsibilities

`VbDriver` consumes the analyzer payload and owns all VB-specific logic:

- active-space selection,
- active/inactive/frozen candidate partitioning,
- VB active orbital construction,
- VB structure generation,
- overlap and Hamiltonian construction,
- generalized eigenvalue solution,
- structure weights and diagnostics.

The analyzer identifies orbitals; the VB driver decides how those orbitals become a wavefunction.

For metal-ligand systems, the analyzer records are intentionally diagnostic: `ML/sigma-acceptor` marks ligand-to-metal sigma donation, while `ML/pi-donor` marks metal-to-ligand pi back-donation into ligand pi-type nonbonding/acceptor space. `VbDriver` can now explicitly select these records as fixed-orbital active-space seeds. The analyzer does not activate them automatically, and NBO still keeps them out of the primary Lewis table.

## Current status: 2026-05-03

`VbDriver` is currently a fixed-orbital VB development driver with shared analyzer input, validated two-electron spin-adapted algebra, frozen-HF embedding, and a determinant-CI fallback for larger π active spaces.

Implemented and validated:

- Shared `OrbitalAnalyzer` consumption for active-orbital selection and active/inactive/frozen candidate diagnostics.
- H₂ one-active-bond VB-CI with covalent/ionic structures, spin-adapted two-electron algebra, generalized overlap metric pruning, and Chirgwin-Coulson/Löwdin weights.
- H₂ one-angle VB-SCF diagnostic path retained as a prototype, while fixed-orbital VB-CI remains the validated baseline.
- Ethylene one-active-π VB-CI with frozen sigma/core HF embedding.
- Fixed-orbital multi-center π active spaces through allyl cation, allyl radical, allyl anion, butadiene, and benzene.
- Orthonormal determinant-CI fallback for larger and open-shell active spaces.
- Chemically readable determinant labels and grouped spatial-occupation resonance diagnostics.
- Graph-generated spin-adapted chemical resonance/CSF template projection diagnostics for allyl radical, allyl anion, butadiene, and benzene.
- Graph-automorphism averaged displayed resonance weights with raw unsymmetrized diagnostics retained.
- Singlet determinant-CI root selection by alpha/beta exchange symmetry before CSF projection.
- CSF template phase convention matched to the determinant-CI alpha/beta occupation-string representation.
- Explicit metal-ligand fixed-orbital active spaces through `active_metal_ligand_channels`, including sigma-only and combined sigma-plus-pi back-donation determinant-CI scans.
- Notebook validation and source-level regression tests for the fixed-orbital π resonance diagnostics, including butadiene captured-subspace weight and singlet exchange parity.

Current boundaries:

- The multi-center π engine is still fixed-orbital; it is not BOVB.
- Larger π systems currently solve a determinant-CI reference and then project chemically compact CSF templates onto that root; this is not yet a true compact CSF Hamiltonian.
- The current BOVB/VBB work is planned but not implemented.
- `VbDriver` owns wavefunction algebra; `OrbitalAnalyzer` only supplies recognition and candidate metadata.

## Implemented capabilities

### H₂ one-active-bond VB-CI

Status: implemented and validated as the two-electron/two-orbital baseline.

Features:

- Automatic generated active space for H₂ when no explicit structures/orbitals are provided.
- Covalent plus two ionic structures.
- Spin-adapted two-electron singlet algebra.
- Metric-pruned generalized eigenproblem.
- Chirgwin-Coulson and Löwdin weights.
- AO ERIs through `FockDriver.compute_eri(..., eri_thresh=1.0e-12)`.
- H₂ dissociation notebook comparison to RHF and broken-symmetry UHF.

Important interpretation:

- Chirgwin-Coulson weights are useful but can become negative or larger than one in a strongly nonorthogonal structure basis.
- Löwdin weights are better for bounded visual diagnostics.

### Ethylene one-active-π VB-CI

Status: implemented and validated as the first polyatomic fixed-orbital checkpoint.

The current ethylene model is still a **2-electron / 2-orbital** active problem:

- active candidate: analyzer C=C `BD(pi)`, e.g. `BD_pi_4`,
- inactive sigma framework: C=C sigma plus four C-H sigma candidates,
- active orbitals: carbon-projected π components of the analyzer candidate,
- structures: π covalent plus two π ionic structures,
- embedding: inactive/core/sigma environment kept at the RHF reference level.

The embedding model is:

```text
inactive density = total RHF density - selected active pair density
```

The inactive density contributes an effective one-electron field and a constant frozen-HF energy. This keeps the sigma framework at the HF level while the selected π pair is treated by the VB-CI structure algebra.

User-verified ethylene checkpoint:

- active-space model: `one-active-pi-bond`
- active candidate: `BD_pi_4` / `pi`
- inactive sigma candidates: `BD_sigma_3`, `BD_sigma_7`, `BD_sigma_9`, `BD_sigma_11`, `BD_sigma_13`
- frozen embedding: `True`
- frozen electrons: about `14`
- active reference electrons: about `2`
- structure weights are finite and symmetric for the two equivalent ionic π structures.

## Why orbital relaxation is not the immediate next step

BOVB introduces structure-specific orbital relaxation. It should now follow the fixed-orbital π validation work, but only after the determinant/CSF boundary and compact resonance diagnostics remain covered by source-level tests.

Before production BOVB, the driver must be able to:

1. select a multi-center π system,
2. construct more than two active π orbitals,
3. support more than two active electrons,
4. generate chemically meaningful π VB structures,
5. handle closed-shell and open-shell spin cases,
6. keep sigma/core/inactive orbitals frozen at the HF reference level,
7. report stable weights and diagnostics.

Those fixed-orbital capabilities are now in place through benzene. The next bridge is a compact spin-adapted CSF Hamiltonian layer, followed by carefully constrained BOVB and VBB method work.

## Fixed-orbital multi-center π engine

### Scope

The implemented ladder generalizes from the ethylene 2e/2π model to fixed-orbital π active spaces of the form:

```text
n_active_electrons / n_active_pi_orbitals
```

The current implementation remains fixed-orbital and embedded in the frozen HF sigma/core density. Orbital optimization and BOVB come later.

### Required engine extensions

1. **Active π orbital selection**
   - Use analyzer `BD(pi)`, pi-capable one-center candidates, and NAO atom/angular maps to identify the π-active atoms.
   - Build one localized π active orbital per active atom or per selected π center.
   - Preserve traceability from each active orbital to atom index, NAO indices, and analyzer candidate labels.

2. **Flexible active electron count**
   - Ethylene: 2e/2π.
   - Allyl cation: 2e/3π.
   - Allyl radical: 3e/3π.
   - Allyl anion: 4e/3π.
   - Butadiene: 4e/4π.
   - Benzene: 6e/6π.

3. **General active determinant or CSF basis**
   - Start with determinant-based active-space CI for robustness.
   - Add spin-adapted CSFs once determinant algebra is stable.
   - Support singlets first, then doublets for radicals.

4. **General overlap and Hamiltonian construction**
   - Remove the current restriction to `(1,1)`, `(2,0)`, and `(0,2)` two-orbital structures.
   - Support arbitrary active occupation strings over non-orthogonal active orbitals.
   - Keep metric pruning and Löwdin weights.

5. **Frozen HF sigma/core embedding**
   - Continue using the total HF density minus active reference density as the inactive environment.
   - Report frozen electron count, active electron count, frozen constant energy, and active candidate labels.

6. **VB structure generation**
   - Generate chemically labeled π structures from graph/connectivity and active electron count.
   - Keep the first generator conservative; add richer ionic and long-bond structures later.

## Validation ladder

### 1. Allyl cation: 2e/3π singlet

This is the recommended next implementation target.

Why first:

- closed-shell singlet,
- only two active π electrons,
- three π centers,
- tests multi-center active orbital construction without open-shell complications,
- tests sigma/core frozen embedding beyond ethylene.

Expected model:

- active atoms: the three allyl carbon centers,
- active orbitals: one π orbital on each carbon,
- inactive framework: C-C sigma and C-H sigma candidates plus core candidates,
- active electron count: `2`,
- spin: singlet,
- minimal structures: left and right π-bond resonance forms, with optional ionic/polar extensions later.

Acceptance checks:

- three active π orbitals are produced,
- frozen electron count equals total electrons minus two,
- sigma candidates remain inactive/frozen,
- resonance structures are symmetry-equivalent for symmetric allyl cation,
- weights are normalized and finite,
- overlap metric rank is stable.

### 2. Allyl radical: 3e/3π doublet

Why second:

- adds open-shell spin handling,
- requires one unpaired π electron,
- tests `SOMO`/spin-density information from `OrbitalAnalyzer`.

Expected model:

- active electron count: `3`,
- spin: doublet,
- resonance forms with the radical center delocalized over terminal carbons,
- frozen sigma/core environment.

Acceptance checks:

- doublet spin sector is selected correctly,
- radical populations are chemically sensible,
- left/right radical resonance weights are symmetric in symmetric geometry.

### 3. Allyl anion: 4e/3π singlet

Why third:

- tests a filled 3-center π system with excess electron density,
- requires lone-pair-like/ionic π structures,
- remains closed-shell but has more active electrons than ethylene or allyl cation.

Expected model:

- active electron count: `4`,
- spin: singlet,
- resonance forms with negative charge delocalized over terminal carbons,
- optional ionic/lone-pair-like π structures.

Acceptance checks:

- active electron count is four,
- charge-delocalized structures are generated and labeled,
- terminal symmetry is respected.

### 4. Butadiene: 4e/4π singlet

Why before benzene:

- bridge from 3-center allyl to cyclic/aromatic systems,
- tests a linear conjugated system with two π bonds,
- simpler than benzene while requiring more than one π pair.

Expected model:

- active electron count: `4`,
- active orbitals: four carbon π orbitals,
- structures: two covalent π bonds plus charge-separated/resonance alternatives.

Acceptance checks:

- four active π orbitals are found,
- sigma framework remains frozen,
- dominant structure is chemically sensible.

### 5. Benzene: 6e/6π singlet

Why after butadiene:

- tests cyclic conjugation and aromatic resonance,
- requires Kekulé and Dewar-style structure generation,
- larger overlap/Hamiltonian matrix.

Expected model:

- active electron count: `6`,
- active orbitals: six carbon π orbitals,
- structures: at least two Kekulé structures, then Dewar structures,
- sigma framework frozen.

Acceptance checks:

- two Kekulé structures are degenerate by symmetry,
- Dewar structures are present and have lower weight in first fixed-orbital model,
- structure weights are normalized and stable,
- no BOVB or orbital relaxation is needed to pass the fixed-orbital checkpoint.

### User-validated fixed-orbital π ladder output

The notebook has now run the fixed-orbital π ladder through benzene with frozen-HF sigma/core embedding:

| System | Active space | Spin | Active orbitals | Structures/determinants | Active reference electrons | Frozen electrons | Status |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| allyl cation | 2e/3π | singlet | 3 | 6 spin-adapted structures | 2.0 | 20.0 | validated |
| allyl radical | 3e/3π | doublet | 3 | 9 determinants | 3.0 | 20.0 | validated |
| allyl anion | 4e/3π | singlet | 3 | 9 determinants | 4.0 | 20.0 | validated |
| butadiene | 4e/4π | singlet | 4 | 36 determinants | 4.0 | 26.0 | validated |
| benzene | 6e/6π | singlet | 6 | 400 determinants | 6.0 | 36.0 | validated |

Interpretation:

- The active/inactive electron partition is stable for all current checkpoints.
- The 2e/3π allyl cation remains the chemically transparent spin-adapted VB benchmark.
- Larger and open-shell systems currently use an orthonormal determinant-CI fallback. This is a robust fixed-orbital active-space checkpoint, but the determinants are not yet compact chemically labeled VB resonance structures.
- Determinant labels now include chemically readable spatial occupations, such as radical center, lone-pair/doubly occupied center, empty center, and open-shell singlet center lists.
- Determinant-CI diagnostics also report grouped resonance-structure weights by spatial occupation. This is the first interpretability layer before true spin-adapted CSF Hamiltonians.
- A second interpretability layer now generates graph-based spin-adapted chemical resonance templates and projects the determinant-CI wavefunction onto them. The current template set covers allyl radical centers, allyl anion lone-pair/π-bond forms, butadiene paired π-bond forms, and benzene Kekulé/Dewar-style pairings.
- For singlet determinant-CI active spaces, the driver now selects the lowest alpha/beta exchange-symmetric root before building resonance diagnostics. This avoids projecting spin-adapted CSF templates onto the exchange-antisymmetric M_S = 0 root that appeared first for butadiene in the determinant-only spectrum.
- CSF template phases are now generated directly in the determinant-CI alpha/beta occupation-string representation instead of applying second-quantized creation-order signs a second time. This fixes the butadiene paired π-bond template capture and keeps the displayed resonance weights compact enough for validation.
- Chemical-resonance diagnostics now separate the compact CSF-template subspace capture from the normalized interpretation within that subspace: `chemical_resonance_subspace_weight` reports how much of the determinant-CI root is captured by the generated templates, while `chemical_resonance_weights` are normalized Löwdin-style weights over the retained template metric. These reported weights are graph-automorphism averaged for the current allyl, butadiene, and benzene template graphs so symmetry-equivalent chemical resonance structures carry equal displayed weights; the raw template-metric values remain available as `chemical_resonance_unsymmetrized_weights` to diagnose orbital/reference symmetry breaking.
- The next scientific improvement should promote these diagnostic CSF projections into an optional compact spin-adapted CSF Hamiltonian where chemically useful, while retaining the determinant-CI fallback as a validation reference.

## Future improvements: compact CSFs, BOVB, and VBB

The fixed-orbital π engine is now stable enough to define the next VB roadmap. The next work should proceed in layers so that each new orbital-relaxation or user-facing feature can be checked against the existing determinant-CI and CSF-projection diagnostics.

### 1. Compact spin-adapted CSF Hamiltonians

Goal: promote the current diagnostic CSF projections into an optional working Hamiltonian basis where chemically useful.

Implementation steps:

1. Reuse the graph-generated CSF templates for allyl, butadiene, and benzene.
2. Build the Hamiltonian and overlap matrices directly in the retained CSF-template metric.
3. Compare CSF-Hamiltonian roots against determinant-CI roots for the same fixed active orbitals.
4. Keep determinant-CI as the reference fallback and as a regression oracle.
5. Report captured-subspace weight whenever a compact CSF basis is used so users can tell whether the compact basis is faithful.

Acceptance checks:

- Butadiene paired π-bond CSFs retain nonzero captured-subspace weight and exchange-symmetric singlet parity.
- Benzene Kekulé/Dewar template counts and symmetry-orbit weights remain stable.
- Compact CSF energies are explicitly labeled as compact-basis approximations when they do not span the full determinant-CI space.

### 2. BOVB: breathing-orbital valence bond

BOVB should be introduced after the fixed-orbital multi-center π engine and the determinant/CSF regression tests remain stable for allyl, butadiene, and benzene checkpoints.

Recommended BOVB order:

1. Revisit ethylene with two π structures and structure-specific breathing orbitals.
2. Apply BOVB to allyl cation after fixed-orbital 2e/3π works.
3. Extend to butadiene and benzene only after the active-space generator and spin handling are stable.

BOVB requirements:

- structure-specific orbital sets,
- constraints to preserve sigma/core frozen orbitals,
- orbital normalization and nonorthogonality controls,
- stable orbital optimization with line search/trust region,
- diagnostics comparing fixed-orbital VB-CI, VB-SCF, and BOVB.

Implementation steps:

1. Start with ethylene and H₂-like two-structure/two-orbital cases, where the fixed-orbital baseline is already transparent.
2. Represent each VB structure with its own allowed active-orbital set while sharing frozen core/sigma orbitals where requested.
3. Add constraints that preserve active-center identity and prevent collapse into unrestricted HF-like orbitals.
4. Optimize orbital parameters and structure coefficients with a robust alternating or joint optimizer.
5. Add diagnostics for structure-specific orbital overlaps, breathing magnitude, orbital gradients, line-search progress, and energy lowering relative to fixed-orbital VB-CI.
6. Only then extend to allyl and butadiene π systems.

Acceptance checks:

- BOVB energy is no higher than the corresponding fixed-orbital VB-CI energy for the same structure set.
- Structure-specific orbitals remain normalized and chemically attached to their intended centers.
- Frozen sigma/core electron count remains stable.
- BOVB output clearly distinguishes fixed-orbital weights, optimized-orbital weights, and any compact CSF projection diagnostics.

### 3. VBB: Valence Bond BOND method

VBB refers to **Valence Bond BOND**, where the second `BOND` stands for **Breathing Orbital Naturally Delocalized**. This is the Linares et al. VBB method and should be documented as a VB electronic-structure method, not as a generic builder/block layer.

The name should stay separate from BOVB:

- BOVB relaxes structure-specific orbitals in breathing-orbital valence bond models.
- VBB uses breathing orbitals that are naturally delocalized according to the Linares et al. method definition.
- A future user-facing constructor/helper layer may still be useful, but it should not be called VBB.

Planned VBB work:

1. **Method definition and notation**
   - Add a concise implementation note for the Linares et al. VBB formalism.
   - Define the relation between VBB structures, breathing orbitals, and naturally delocalized orbital constraints.
   - Keep VBB terminology distinct from generic VB input builders.

2. **VBB orbital model**
   - Introduce naturally delocalized breathing-orbital parameterizations for selected VB structures.
   - Preserve active-center identity, spin coupling, and frozen sigma/core constraints.
   - Compare the allowed VBB orbital space to the corresponding BOVB and fixed-orbital VB spaces.

3. **VBB Hamiltonian and optimization**
   - Reuse the compact CSF or determinant reference machinery for fixed-orbital checks.
   - Add VBB-specific orbital-optimization variables and constraints only after BOVB/fixed-orbital baselines are stable.
   - Report orbital-delocalization diagnostics, breathing amplitudes, structure weights, and energy lowering relative to fixed-orbital VB-CI and BOVB where applicable.

4. **Validation and reporting**
   - Start from two-orbital H₂/ethylene-style benchmarks before applying VBB to allyl, butadiene, and benzene.
   - Verify that VBB energies and diagnostics are labeled separately from BOVB and compact CSF results.
   - Refuse ambiguous VBB specifications unless the Linares et al. constraints can be applied deterministically.

Acceptance checks:

- H₂ and ethylene VBB prototypes reproduce the expected Linares et al. limiting behavior for the chosen benchmark definitions.
- Allyl, butadiene, and benzene VBB extensions remain traceable to the current fixed-orbital π ladder and compact CSF templates.
- Output clearly labels fixed-orbital VB-CI, compact CSF, BOVB, and VBB quantities as distinct models.

### 4. Metal-ligand VB active spaces

`OrbitalAnalyzer` now exposes metal-ligand candidate records and `VbDriver` supports optional fixed-orbital VB active spaces involving ligand donation and back-donation.

Current status:

- `active_metal_ligand_channels=("sigma-acceptor",)` builds a sigma-only donor/acceptor active space.
- `active_metal_ligand_channels=("sigma-acceptor", "pi-donor")` builds a combined sigma-plus-back-donation active space.
- The combined model is currently a determinant-CI active space with the sigma donor/acceptor pair and the pi donor/acceptor pair active.
- The notebook `docs/metal_ligand_recognition.ipynb` scans Pd--NH3 and Pd--PH3 distances and compares `E_sigma_only` with `E_sigma_plus_pi`.

Planned steps:

- Treat metal-ligand candidates as analyzer-provided active-orbital suggestions, not as automatic selected VB structures.
- Add a pedagogical dissociation sequence for metal-ligand bonds analogous to H2: RHF, broken-symmetry UHF, fixed-orbital VB-CI, and later BOVB/VBB correlation-orbital relaxation.
- Add metal-ligand VBB benchmarks only after the organic VBB method layer is defined and validated.
- Keep determinant-CI fallback available because compact spin coupling for metal centers may be system dependent.
- Add diagnostics separating ligand-field interpretation from organic π resonance interpretation.

## Immediate implementation checklist

Current fixed-orbital π implementation status:

- [x] Document current H₂ and ethylene status.
- [x] Define the fixed-orbital π active-space roadmap.
- [x] Add an `active_pi_atoms` selection option to `VbComputeOptions`.
- [x] Build a multi-center `VbActiveSpace` representation with active electron count, spin, and arbitrary active π orbitals.
- [x] Implement active π orbital construction for more than two atoms.
- [x] Keep the 2e/3π allyl cation as a spin-adapted two-electron VB checkpoint.
- [x] Add an orthonormal determinant-CI fallback for larger/open-shell fixed-orbital π spaces.
- [x] Add chemically readable determinant labels and grouped spatial-occupation resonance weights.
- [x] Add graph-based spin-adapted resonance/CSF projection diagnostics for allyl radical, allyl anion, butadiene, and benzene.
- [x] Extend frozen-HF embedding from one selected pair density to a general active reference density.
- [x] Add notebook cells for allyl cation and the validation ladder through benzene.
- [x] Add source-level smoke tests for allyl cation, allyl radical, allyl anion, butadiene, and benzene.
- [x] Add source-level regression tests for resonance/CSF labels, counts, projection weights, butadiene captured-subspace weight, and singlet exchange parity.
- [ ] Replace selected determinant-CI diagnostics by true spin-adapted CSF Hamiltonians where chemically useful.
- [ ] Add BOVB prototype for a two-orbital benchmark.
- [ ] Add VBB method notes and prototypes for H₂/ethylene-style benchmarks, then extend to allyl, butadiene, and benzene.
- [x] Add fixed-orbital metal-ligand active-space support after analyzer metal-ligand recognition exists.
- [ ] Add H2 and metal-ligand dissociation notebooks comparing RHF, UHF, fixed-orbital VB-CI, and future BOVB correlation/orbital-relaxation effects.

## Tomorrow restart notes

Start from the validated fixed-orbital baseline, not from BOVB/VBB:

1. Re-open `docs/metal_ligand_recognition.ipynb` for the current Pd--NH3/Pd--PH3 real-SCF workflow and the sigma-only versus sigma-plus-pi VB-CI scan.
2. The metal-ligand entry points to remember are `active_metal_ligand_channels=("sigma-acceptor",)` for the sigma-only model and `active_metal_ligand_channels=("sigma-acceptor", "pi-donor")` for the combined sigma/back-donation model.
3. The best next documentation/example task is a clean dissociation comparison: RHF, broken-symmetry UHF, fixed-orbital VB-CI, and clearly marked future BOVB/VBB placeholders.
4. The best next implementation task is still compact spin-adapted CSF Hamiltonians for selected organic π examples, using determinant-CI as the regression oracle.
5. Do not claim BOVB, VBB, or production metal-ligand VB is implemented. The current metal-ligand model is a fixed-orbital determinant-CI prototype seeded by analyzer diagnostics.

Current safe stopping point: fixed-orbital H2, ethylene, π-ladder, and metal-ligand sigma/sigma-plus-pi active-space diagnostics are in place; orbital relaxation and production coordination VB remain roadmap items.

## Documentation and notebook policy

- The notebook is for validation and interpretation, not for core logic.
- All active-space construction and Hamiltonian logic should live in `VbDriver`.
- Notebook cells should print concise diagnostics:
  - active model,
  - active electron/orbital counts,
  - active atoms and labels,
  - inactive sigma/core labels,
  - frozen electron count,
  - energy,
  - overlap rank/eigenvalues,
  - weights by structure.

## Current short-term decision

The initial fixed-orbital multi-center π engine is now in place through benzene. The immediate next VB step is to keep the determinant-CI fallback as a reference while promoting selected CSF diagnostics into compact spin-adapted Hamiltonians. BOVB and VBB should follow only with explicit regression tests against the fixed-orbital baseline.
