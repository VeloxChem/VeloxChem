# NboDriver Technical Notes

These notes describe the implemented algorithms in VeloxChem's `NboDriver`. The emphasis is on reproducibility, explicit electron accounting, and clear separation between candidate generation, Lewis assignment, and NRA/NRT density fitting.

## Pipeline

The computation proceeds through the following layers:

```text
AO overlap and SCF density
-> NAO orthonormalization and NPA
-> MO analysis in the NAO basis
-> NBO candidate generation
-> primary Lewis assignment
-> active-space Lewis/resonance alternatives
-> optional NRA/NRT density fitting
-> reports and diagnostics
```

The AO overlap, density, NAO/NPA, spin data, MO-in-NAO analysis, and candidate-construction layers are obtained through the shared `OrbitalAnalyzer` payload. The NBO driver owns the Lewis assignment, resonance alternatives, donor-acceptor diagnostics, NRA/NRT fitting, and reports. The layers are intentionally one-way. NRA/NRT weights are computed after alternatives are available and do not alter the candidate pool or primary assignment.

## NAO/NPA construction

The input AO overlap matrix is `S`, and the SCF AO density is `P`. The NAO transformation `T` is assembled from atom-local AO blocks and canonicalized so equivalent atom-block rotations are deterministic. The defining invariant is

```text
T^T S T = I.
```

The total NAO density is

```text
D_NAO = T^T S P S T.
```

For unrestricted references, alpha and beta densities are transformed separately:

```text
D_alpha_NAO = T^T S P_alpha S T
D_beta_NAO  = T^T S P_beta  S T
D_spin_NAO  = D_alpha_NAO - D_beta_NAO.
```

The implementation records diagnostics for orthonormality, density symmetry, electron count, trace conservation, and charge conservation. NPA populations are atom-block traces of `D_NAO`; NPA charges are nuclear charge minus the corresponding atom population.

## MO-in-NAO analysis

Molecular orbitals are projected into the orthonormal NAO basis. For restricted references the occupied and virtual orbital compositions are reported from the closed-shell coefficient set. For unrestricted references, alpha and beta channels are handled independently and the spin density is retained for radical diagnostics.

MO composition is diagnostic. It does not determine the Lewis assignment directly, but it provides interpretable orbital character in the same basis used for NBO candidate vectors.

## Candidate construction

Candidates are normalized vectors in NAO space with explicit metadata:

```text
n(c) = c^T D_NAO c
m(c) = c^T D_spin_NAO c.
```

The shared analyzer candidate layer generates these candidate families for both NBO and VB consumers:

| Family | Role |
| --- | --- |
| `CR` | Core orbitals from atom-local occupied NAO space. |
| `LP` | Lone-pair candidates from one-center valence density blocks. |
| `SOMO` | Singly occupied open-shell candidates. |
| `BD(sigma)` | Two-center sigma bonding candidates. |
| `BD(pi)` | Two-center pi bonding candidates. |
| one-electron `LP`/`BD(pi)` | Radical active-space candidates with `electron_count = 1`. |
| `BD*` | Antibonding complements to bonding candidates. |
| `RY` | One-center Rydberg/acceptor complements with low occupation. |
| `ML` | Metal-ligand diagnostic records emitted by `OrbitalAnalyzer`; not selected as occupied organic Lewis NBOs. |

Occupied Lewis assignment uses candidates with explicit electron counts. Candidate-only acceptors (`BD*` and `RY`) are retained for diagnostics and reports but are not inserted into occupied Lewis structures. VB workflows consume the same candidate records as labels and initial active-orbital suggestions, without depending on NBO Lewis-assignment internals.

Two-center candidates are built from atom-pair density blocks. Sigma and pi candidates are kept separate in metadata and downstream result fields. Pi candidates enter the resonance active space; sigma bonds and non-participating lone pairs normally form the fixed framework unless constraints or electron accounting require otherwise.

## Primary Lewis assignment

Primary assignment selects a chemically conservative occupied NBO list with exact electron counting. The selected structure includes:

- `nbo_list`: all occupied objects in report order.
- `sigma_nbo_list`: core, sigma, and fixed non-pi framework objects.
- `pi_nbo_list`: occupied pi objects.
- `fixed_nbo_list`: occupied objects held fixed during resonance enumeration.
- `active_nbo_list`: occupied objects considered in the resonance active space.
- `active_pi_nbo_list`: active pi objects.
- `active_lone_pair_nbo_list`: active lone-pair donation objects.
- `active_one_electron_nbo_list`: radical one-electron active objects.

Lewis accounting is evaluated per atom. Each selected candidate contributes its `electron_count` to its participating atoms according to its type. The driver records valence electron counts, formal-charge diagnostics, duet/octet deviations, missing/extra electron diagnostics, and score terms.

## Constraints

`NboConstraints` can require or forbid sigma bonds and can require or restrict pi bonds. Constraints are normalized to atom-index tuples and are applied during candidate selection and alternative enumeration. Violations are surfaced through score terms rather than hidden state.

## Lewis/resonance alternatives

Alternatives are generated from the active space by exact electron-count enumeration. The fixed framework contributes a known electron count, and active choices are accepted only if the target active electron count is recovered. This allows closed-shell pi resonance, lone-pair donation, and radical one-electron alternatives to use the same mechanism.

Alternatives are ranked by an explicit Lewis score. The score combines candidate occupations, bond/constraint satisfaction, formal-charge diagnostics, octet/duet diagnostics, and consistency penalties. The ranking weight is

```text
w_score,k = exp(-beta s_k) / sum_j exp(-beta s_j),
```

where `s_k` is the Lewis score and `beta` is `NboComputeOptions.lewis_weight_beta`.

These weights are useful for ordering alternatives. They are not NRA/NRT density weights and are not VB amplitudes.

## Resonance signatures and classes

Each alternative receives a molecule-independent structured resonance signature built from the one-based report fields:

```text
pi_bonds
active_lone_pair_atoms
active_positive_atoms
active_one_electron_atoms
```

The same information is stored as a compact display label such as `pi:1-2,3-4,5-6` or `pi:1-2; lp:3; pos:2`. These labels are report annotations and prior-matching handles; the actual data model remains the explicit fixed/active and sigma/pi partitions.

Symmetry-equivalent alternatives are grouped with conservative graph automorphisms that preserve element labels and connectivity. Each alternative stores:

- `resonance_signature`
- `resonance_label`
- `class_signature`
- `class_label`
- `resonance_class_id`
- `class_degeneracy`

The top-level `results["resonance_classes"]` list reports one record per class with member ranks, member labels, degeneracy, and class-level score/ranking weight summaries:

- `class_weight_sum`
- `class_weight_mean`
- `class_weight_min`
- `class_weight_max`

For benzene with only ring pi bonds allowed, the two localized Kekule alternatives have distinct resonance signatures but share one class with degeneracy 2. This separates localized alternative ranking from chemically equivalent class interpretation.

## NRA/NRT density fitting

When requested, NRA/NRT constructs a model density for each alternative and fits the SCF density as a convex combination of those model densities. For a selected NAO subspace, the closed-shell problem is

```text
min_w || d - A w ||^2
subject to w_k >= 0, sum_k w_k = 1.
```

`d` is the vectorized target density in the chosen subspace, and column `A_k` is the vectorized model density for alternative `k`. The implementation reports the fitted weights, residual norm, rank diagnostics, structure metadata, and the subspace used.

Supported subspaces are:

| Subspace | Description |
| --- | --- |
| `selected` | NAOs touched by the selected alternatives. |
| `pi` | Pi-active NAO support. |
| `valence` | Non-core valence NAO support. |
| `full` | Full NAO density. |

The fit uses nonnegative normalized weights. The reported `nra_weight` values are density-reconstruction weights, not wavefunction weights.

## Prior-regularized NRA/NRT

Optional priors add a regularization target without replacing the density objective:

```text
min_w || d - A w ||^2 + lambda || w - w0 ||^2
subject to w_k >= 0, sum_k w_k = 1.
```

Priors can be matched by alternative index or recognized label. The normalized prior vector, matched structures, unmatched entries, mode, and regularization strength are stored in `results["nra"]["prior"]`. Each structure may also carry `prior_weight`.

The regularization is deliberately explicit because prior-guided weights are still density-fit weights. A nonzero prior can influence the fitted distribution, but the residual and model densities remain visible for interpretation.

## Open-shell treatment

For unrestricted references, the driver keeps total and spin NAO densities. Radical candidates carry one-electron `electron_count` metadata, and doublet active spaces are constrained so alternatives preserve the intended radical electron count.

Open-shell NRA/NRT with `nra_spin_fit="total_spin"` concatenates total-density and spin-density residuals:

```text
min_w || [d_total; d_spin] - [A_total; A_spin] w ||^2
subject to w_k >= 0, sum_k w_k = 1.
```

This gives alternatives weights that reconstruct both charge distribution and spin distribution in the selected subspace.

## Donor-acceptor diagnostics

`BD*` candidates are constructed as orthogonal antibonding complements in the same two-center subspace as their parent `BD` candidates. `RY` candidates are one-center low-occupation complements. Both remain candidate-only acceptors.

For each occupied donor and candidate-only acceptor, the driver reports a density-coupling diagnostic in the NAO basis. This ranking identifies plausible donation channels and keeps zero or small channels visible when useful for regression testing.

The diagnostic is not the mature NBO second-order perturbation expression. That expression requires an NBO Fock matrix and energy denominator, which are outside the current implementation.

In full NBO reports, `BD*` and `RY` acceptors are printed in a separate candidate-only section. This keeps them visible for interpretation while preserving the invariant that selected occupied Lewis tables contain only occupied Lewis NBOs.

## Metal-ligand diagnostics

`OrbitalAnalyzer` can emit metal-ligand `ML` records for coordination diagnostics, including ligand-to-metal sigma donation and metal-to-ligand pi/back-donation channels. `NboDriver.compute()` exposes these records through `metal_ligand_diagnostics`, and `nbo_report(level="full")` prints them in a diagnostic-only metal-ligand section. They are not used in occupied organic Lewis assignment.

The current boundary is deliberate:

- organic Lewis assignment remains governed by explicit electron counts, valence diagnostics, and duet/octet checks;
- `ML` records are diagnostic candidate records until their electron-count semantics are defined for NBO selection;
- current full reports present metal center, ligand atom, coordination mode, channel, occupation, donation strength, back-donation strength, and role separately from the selected organic Lewis table;
- future reports should add richer d-manifold character and metal-aware electron-accounting interpretation.

## Reporting and labels

Reports are generated from structured result dictionaries. The primary report presents foundation diagnostics, selected NBOs, candidate counts, Lewis accounting, resonance classes, Lewis/resonance alternatives, and donor-acceptor diagnostics. The NRA report presents structure weights, residuals, spin-fit metadata, and prior metadata.

Structure signatures are generated directly from the pi-bond and active-center metadata, for example `pi:1-2,3-4` or `pi:1-2; lp:3; pos:2`. They are useful for display, prior matching, and resonance-class grouping, but the algorithmic representation remains sigma/pi and fixed/active partition fields. The implementation does not contain molecule-specific structure-name classifiers.

## Validation coverage

The regression suite covers:

- NAO orthonormality, charge conservation, and deterministic canonicalization.
- Restricted and unrestricted MO-in-NAO analysis.
- Candidate counts and candidate metadata.
- Required/forbidden bond and pi-bond constraints.
- Sigma/pi partition fields and fixed/active fields.
- Lewis electron accounting, formal-charge diagnostics, octet diagnostics, and score terms.
- Lone-pair donation alternatives.
- Resonance signatures, symmetry-equivalent classes, class-level weights, and final rank ordering.
- Polar-resonance metadata needed by `show_structures(...)`.
- Closed-shell NRA/NRT structure pools and reports.
- NRA/NRT prior regularization and metadata.
- Open-shell SOMO candidates, one-electron alternatives, and total-plus-spin NRA/NRT.
- `BD*`, `RY`, donor-acceptor diagnostics, candidate-only acceptor report sections, and diagnostic-only `ML` report sections.

## Current status: 2026-05-04

`NboDriver` is currently a structured NBO/NRA/NRT analysis driver built on the shared `OrbitalAnalyzer` payload.

Implemented and validated:

- Shared analyzer integration for AO maps, NAO/NPA data, MO-in-NAO analysis, spin-density diagnostics, and candidate records.
- Conservative occupied NBO selection with explicit electron-count accounting.
- Sigma/pi/fixed/active partitions for Lewis and resonance analysis.
- General pi-resonance alternatives for organic conjugated systems, including closed-shell, radical, anionic, and polar pi cases.
- Candidate-only `BD*` and `RY` acceptors plus donor-acceptor density-coupling diagnostics.
- Resonance signatures, symmetry-equivalent resonance classes, class degeneracies, and class-level score/ranking weight summaries.
- Full-report candidate-only acceptor sections for `BD*` and `RY` diagnostics.
- Full-report diagnostic-only metal-ligand sections for analyzer-owned `ML` records.
- Closed-shell NRA/NRT density fitting with selected, pi, valence, and full subspaces.
- Prior-regularized NRA/NRT with explicit prior metadata.
- Open-shell total-density plus spin-density NRA/NRT fitting.
- Structured reports and regression coverage for the current pipeline.

Metal-ligand checkpoint for the current notebook work:

- `OrbitalAnalyzer` emits `ML/sigma-acceptor` and `ML/pi-donor` diagnostic records through the shared payload.
- `NboDriver.compute()` exposes these records through `metal_ligand_diagnostics`, and `nbo_report(level="full")` prints them in a diagnostic-only section.
- The real Pd--NH3/Pd--PH3 notebook now treats B3LYP constrained-scan and HF single-point total energies as the validated reference dissociation curves, plotted as `E(R) - E(5.0 Angstrom)`.
- The validated grid dissociation energies are taken only from those reference total-energy curves. Metal-ligand VB-SCF/BOVB sigma values are kept outside NBO interpretation and are documented as VB diagnostics, not NBO or NRT observables.

Important boundaries:

- NBO score/ranking weights are Lewis-ranking diagnostics, not density-fit weights.
- NRA/NRT weights are density-reconstruction weights, not VB wavefunction amplitudes.
- Donor-acceptor diagnostics are not yet second-order perturbation energies because an NBO Fock matrix and orbital-energy denominator layer are not implemented.
- Metal-ligand chemistry is recognized as analyzer-owned diagnostic `ML` records, exposed through `metal_ligand_diagnostics`, and printed in full NBO reports. Metal-aware electron-accounting interpretation is not yet implemented.

## Scientific scope

The implemented driver is a coherent, documented NBO/NRA/NRT analysis layer. It intentionally keeps several advanced NBO concepts out of scope until their required mathematical layers are present: NHO directionality, NBO Fock-matrix perturbation energies, broad automatic NRT labeling, and VB/wavefunction state weights.

## Future improvements

The remaining NBO work should focus on interpretation, reporting, validation breadth, and mathematically distinct NBO layers that are not present yet.

1. **Metal-ligand NBO reporting and interpretation**
	- Consume analyzer `ML` records for sigma donation, pi donation, and back-donation diagnostics.
	- Extend the initial full-report coordination section with d-manifold character and richer ligand-channel interpretation.
	- Keep `ML` diagnostics separate from occupied organic Lewis tables until their electron-count semantics are explicit.

2. **NHO and directional hybrids**
	- Add a natural hybrid orbital layer on top of NAO-space candidates.
	- Provide stable hybrid directionality and percent s/p/d character.
	- Preserve deterministic behavior for nearly degenerate local rotations.

3. **NBO Fock matrix and second-order donor-acceptor energies**
	- Transform a Fock-like operator into the NBO basis.
	- Add donor-acceptor perturbation energies only after orbital energies and denominators are well defined.
	- Keep the current density-coupling diagnostic available as a separate lightweight ranking.

4. **Broader validation and examples**
	- Extend regression tests for benzene, allyl ions/radicals, nitro systems, zwitterions, amides, and transition-metal fragments.
	- Keep notebooks as human-facing demonstrations, not as the only source of validation.
