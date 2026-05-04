# NboDriver Theory and Implementation Connection

This note describes the theoretical model behind VeloxChem's `NboDriver` and maps each concept to the implemented data structures. It is meant to sit between the user-facing API guide and the developer implementation notes: the focus is scientific interpretation, with enough implementation detail to make the equations traceable in code and result dictionaries.

## Conceptual scope

Natural Bond Orbital analysis rewrites an SCF density in a chemically localized orbital language. The central idea is not to change the underlying wavefunction, but to express its one-particle density in orbitals that resemble core orbitals, lone pairs, bonds, antibonds, Rydberg acceptors, and resonance alternatives.

`NboDriver` implements this interpretation as a density-based pipeline:

```text
SCF density -> NAO/NPA -> NBO candidates -> Lewis alternatives -> NRA/NRT density weights
```

The implementation is deliberately explicit about what each layer means:

- NAO/NPA is a population and orbital-basis transformation layer.
- NBO candidates are localized density-block natural orbitals in NAO space.
- Lewis alternatives are electron-counted selections from the occupied candidate pool.
- NRA/NRT weights are convex density-reconstruction weights over alternatives.
- Donor-acceptor diagnostics are density-coupling diagnostics to unoccupied candidate-only acceptors.

This separation prevents three common ambiguities: a Lewis ranking weight is not an NRA/NRT density weight, an NRA/NRT density weight is not a VB wavefunction coefficient, and a candidate-only acceptor is not part of an occupied Lewis structure.

## AO density to NAO density

The SCF calculation provides an AO overlap matrix `S` and an AO density matrix `P`. Natural Atomic Orbitals are represented by a transformation matrix `T` that orthonormalizes atom-local orbital blocks:

```text
T^T S T = I.
```

The density in the NAO basis is

```text
D_NAO = T^T S P S T.
```

For unrestricted references, the alpha and beta densities are transformed separately:

```text
D_alpha_NAO = T^T S P_alpha S T
D_beta_NAO  = T^T S P_beta  S T
D_spin_NAO  = D_alpha_NAO - D_beta_NAO.
```

The total density describes charge distribution. The spin density describes where the unpaired-electron character resides.

### Implementation connection

The NAO layer appears in `results["nao"]`. The NPA layer appears in `results["npa"]`. The implementation records invariants such as orthonormality, density symmetry, electron count, trace conservation, and charge conservation, so the theoretical transformation can be checked numerically for every run.

The same NAO density is the source for candidate occupations, Lewis accounting, donor-acceptor diagnostics, and NRA/NRT target vectors.

## Natural populations and charges

Natural Population Analysis partitions the NAO density into atom-local populations. For atom `A`, the population is the trace of the NAO density over the NAOs belonging to that atom:

```text
N_A = Tr_A(D_NAO).
```

The corresponding NPA charge is

```text
q_A = Z_A - N_A,
```

where `Z_A` is the nuclear charge. This is a density partition, not a formal oxidation-state assignment. It reflects how the SCF density is distributed over the NAO basis.

### Implementation connection

NPA charges and populations are stored in `results["npa"]`. Lewis formal-charge diagnostics later in the pipeline are separate bookkeeping quantities based on selected NBO electron counts; they should not be confused with NPA charges.

## NBO candidates as density-block natural orbitals

An NBO candidate is a normalized vector `c` in NAO space. Its total and spin occupations are evaluated as

```text
n(c) = c^T D_NAO c
m(c) = c^T D_spin_NAO c.
```

The implementation constructs candidates by diagonalizing chemically meaningful density blocks:

- one-center blocks for core and lone-pair candidates;
- two-center blocks for sigma and pi bond candidates;
- open-shell blocks for singly occupied and one-electron radical candidates;
- complementary subspaces for antibonding `BD*` and Rydberg-like `RY` acceptors.

The theory is density-local: candidates are chosen because the SCF one-particle density has large or interpretable natural occupations in those localized subspaces.

### Implementation connection

Candidates are stored in `results["nbo_candidates"]`. Each candidate record carries its type, atoms, electron count, occupation, vector coefficients, and metadata such as bond kind or parentage.

Occupied Lewis candidates include `CR`, `LP`, `SOMO`, `BD(sigma)`, `BD(pi)`, and radical one-electron candidates. Candidate-only acceptors include `BD*` and `RY`; these are retained for donor-acceptor diagnostics and reports but excluded from occupied Lewis alternatives.

## Sigma and pi separation

The driver keeps sigma and pi content as explicit partitions rather than assigning broad named structure classes internally. A sigma framework contains core orbitals, sigma bonds, and fixed lone-pair content. The pi active space contains pi bonds and participating lone pairs or one-electron radical objects.

The theoretical reason is simple: resonance is usually a redistribution of active pi and lone-pair electron density over a largely fixed sigma skeleton. Keeping the partition explicit preserves the chemical model without imposing labels too early.

### Implementation connection

Primary assignments and alternatives expose:

- `sigma_nbo_list`
- `pi_nbo_list`
- `fixed_nbo_list`
- `active_nbo_list`
- `active_pi_nbo_list`
- `active_lone_pair_nbo_list`
- `active_one_electron_nbo_list`

Reports may show generic pi-bond signatures such as `pi:1-2,3-4`, but the algorithmic representation remains these sigma/pi and fixed/active fields.

## Lewis electron accounting

A Lewis alternative is an electron-counted selection of occupied candidates. Each selected candidate contributes an explicit number of electrons:

```text
CR, LP, BD        -> typically 2 electrons
SOMO, radical LP -> 1 electron
radical BD(pi)   -> 1 electron
```

The selected candidates determine atom-local Lewis valence counts, formal-charge diagnostics, and duet/octet deviations. These diagnostics are bookkeeping checks on the selected Lewis model, not direct observables.

A valid alternative must recover the required active electron count once the fixed framework is included. This is the key constraint that allows closed-shell resonance, lone-pair donation, and radical alternatives to share the same machinery.

### Implementation connection

The selected structure is stored in `results["primary_assignment"]`. Resonance alternatives are stored in `results["alternatives"]`. Each alternative carries lists of selected candidates, atom accounting, score terms, and a Lewis score.

Constraints supplied through `NboConstraints` guide this selection by requiring or forbidding bonds and by restricting pi-bond choices.

## Lewis scores and ranking weights

The Lewis score is a transparent ranking objective. It combines terms such as candidate occupation, required/forbidden bond satisfaction, formal-charge balance, duet/octet diagnostics, and active-space consistency. If `s_k` is the score for alternative `k`, the score-derived ranking weight is

```text
w_score,k = exp(-beta s_k) / sum_j exp(-beta s_j),
```

where `beta` is `NboComputeOptions.lewis_weight_beta`.

This is a ranking model. It says which Lewis alternatives are preferred by the implemented scoring objective. It does not say how much of the SCF density is reconstructed by an alternative, and it is not a wavefunction expansion coefficient.

### Implementation connection

The score and score terms are carried by the primary assignment and by alternatives. The score-derived alternative weight is the `weight` field on alternatives. NRA/NRT weights are stored separately under `results["nra"]` and as `structures[*]["nra_weight"]`.

## NRA/NRT as density reconstruction

Natural Resonance Analysis / Natural Resonance Theory is implemented as a density-fitting problem over Lewis alternatives. Each alternative produces a model density vector `A_k` in a chosen NAO subspace. The SCF target density is `d`. The closed-shell fit is

```text
min_w || d - A w ||^2
subject to w_k >= 0, sum_k w_k = 1.
```

The constraints make the fit a convex mixture of alternative densities. The resulting weights answer a density question: how should the available alternatives be combined to reconstruct the SCF density in the chosen subspace?

The subspace can be the selected active support, the pi support, the valence support, or the full NAO density. A smaller subspace emphasizes the resonance region; a larger subspace makes the reconstruction more global.

### Implementation connection

NRA/NRT is enabled with `NboComputeOptions(include_nra=True)`. The subspace is selected through `nra_subspace`, with values such as `"selected"`, `"pi"`, `"valence"`, and `"full"`.

The fit results are stored in `results["nra"]`. Structure-specific weights appear as `structures[*]["nra_weight"]`. The residual norm, subspace metadata, structure labels, and rank diagnostics are reported so the user can judge whether the fitted weights are chemically meaningful for the chosen model space.

## Prior-regularized NRA/NRT

A prior can be used when a user has external information about expected resonance weights. The prior does not replace the density target. It modifies the objective:

```text
min_w || d - A w ||^2 + lambda || w - w0 ||^2
subject to w_k >= 0, sum_k w_k = 1.
```

Here `w0` is the normalized prior vector and `lambda` is the prior strength. As `lambda` increases, the fit is pulled toward the prior; as `lambda` approaches zero, the fit returns to the density-only solution.

### Implementation connection

Priors are supplied through `NboComputeOptions.nra_prior_weights` and `NboComputeOptions.nra_prior_strength`. The prior metadata is stored in `results["nra"]["prior"]`, and each matched structure may carry `prior_weight`.

A prior-guided `nra_weight` remains a density-fit weight. It should be interpreted together with the reported residual and prior strength.

## Open-shell total-plus-spin fitting

For radicals and other unrestricted references, total density alone may not distinguish where the unpaired electron resides. The implemented open-shell NRA/NRT fit therefore includes both total-density and spin-density targets:

```text
min_w || [d_total; d_spin] - [A_total; A_spin] w ||^2
subject to w_k >= 0, sum_k w_k = 1.
```

This connects radical Lewis alternatives to both charge reconstruction and spin reconstruction. The spin term is essential when alternatives differ primarily by the location of the unpaired electron.

### Implementation connection

Open-shell fitting is controlled by `NboComputeOptions.nra_spin_fit`, with `"total_spin"` as the implemented mode. Spin-resolved NAO densities are stored in the NAO layer, and one-electron active candidates appear in `active_one_electron_nbo_list`.

## Donor-acceptor diagnostics

Classical NBO analysis often interprets delocalization through donor-to-acceptor interactions. In the current driver, occupied NBOs act as donors and candidate-only `BD*` or `RY` orbitals act as acceptors. The implemented coupling is a density-based diagnostic in NAO space.

This diagnostic identifies plausible delocalization channels, but it is not a second-order perturbation energy. A true perturbative donor-acceptor energy would require an NBO Fock matrix and orbital-energy denominator.

### Implementation connection

Diagnostics are stored in `results["donor_acceptor_diagnostics"]`. They connect occupied donor candidates to `BD*` and `RY` acceptors. The report layer can display these channels, but the acceptors remain outside the occupied Lewis lists.

## Interpreting the three weight types

The driver deliberately separates three weight concepts:

| Weight | Location | Interpretation |
| --- | --- | --- |
| Lewis ranking weight | `alternatives[*]["weight"]` | Softmax ranking from the Lewis score. |
| NRA/NRT density weight | `results["nra"]["weights"]`, `structures[*]["nra_weight"]` | Convex density-reconstruction coefficient. |
| VB/wavefunction weight | `VbDriver` results | State-coupling quantity owned by the VB driver; not computed or reinterpreted by NBO/NRA/NRT. |

Keeping these weights separate is part of the scientific contract of the implementation. It avoids presenting a density fit or a Lewis ranking heuristic as a wavefunction population.

## Reading results as theory objects

A useful way to read an `NboDriver` result is by layer:

1. Check `results["nao"]` and `results["npa"]` to verify the density transformation and atom populations.
2. Inspect `results["nbo_candidates"]` to see the localized orbital possibilities implied by the density.
3. Read `results["primary_assignment"]` as the selected Lewis model.
4. Read `results["alternatives"]` as electron-counted resonance possibilities.
5. Read `results["nra"]` as the density-reconstruction fit over those possibilities.
6. Read `results["donor_acceptor_diagnostics"]` as a diagnostic map of donor-to-acceptor delocalization channels.

This layer-by-layer interpretation is the most faithful connection between the theory and the implementation.

## Current boundaries

The implemented theory is a density-based NBO/NRA/NRT model. The following concepts are intentionally outside the current mathematical layer:

- NHO hybrid directionality as a separate hybrid-orbital construction.
- NBO Fock-matrix second-order donor-acceptor energies.
- General automatic classification of every resonance pattern into named chemical classes.
- VB or wavefunction state-mixing weights.
- Metal-ligand `ML` records are diagnostic orbital-recognition objects. They can be used to follow ligand-to-metal sigma donation and metal-to-ligand pi back-donation in notebooks, but they are not yet coordination-aware Lewis/NRT structures and they are not VB dissociation energies.

These boundaries do not weaken the implemented model; they make clear which quantities are available from the present density analysis and which require additional theoretical machinery.

## Current metal-ligand interpretation checkpoint

The current Pd--NH3/Pd--PH3 notebook separates three layers that should not be conflated:

1. B3LYP constrained-scan and HF single-point total energies provide the validated reference potential-energy curves, plotted as `E(R) - E(5.0 Angstrom)`.
2. Analyzer/NBO `ML/sigma-acceptor` and `ML/pi-donor` records provide density-based channel diagnostics along the same geometries.
3. VB-SCF/BOVB sigma metal-ligand traces are downstream wavefunction diagnostics and are not currently accepted as validated dissociation-energy curves.

This separation is the current scientific state: NBO owns density interpretation, VB owns wavefunction diagnostics, and only total-energy reference methods are presently used for quantitative Pd--ligand dissociation energies.
