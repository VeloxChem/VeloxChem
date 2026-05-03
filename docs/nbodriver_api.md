# NboDriver API Guide

`NboDriver` provides a structured Natural Bond Orbital analysis pipeline for VeloxChem wavefunctions. The driver consumes the shared `OrbitalAnalyzer` payload for Natural Atomic Orbitals (NAOs), Natural Population Analysis (NPA) data, MO-in-NAO diagnostics, and NBO-like candidates, then performs Lewis/resonance assignment, donor-acceptor diagnostics, and optional NRA/NRT density-fit weights.

The API is intentionally transparent: every assignment carries explicit electron counts, atom participation, score terms, and density-fit metadata. This makes the results suitable for regression testing, method development, and chemical interpretation without hiding the distinction between Lewis ranking weights, density-fit weights, and future wavefunction weights.

## Basic use

```python
import veloxchem as vlx
from veloxchem.nbodriver import NboDriver, NboComputeOptions

molecule = vlx.Molecule.read_molecule_string("""
O  0.000000  0.000000  0.000000
H  0.758602  0.000000  0.504284
H -0.758602  0.000000  0.504284
""", units="angstrom")

basis = vlx.MolecularBasis.read(molecule, "sto-3g")
scf = vlx.ScfRestrictedDriver()
scf.compute(molecule, basis)

nbo = NboDriver()
results = nbo.compute(
    molecule,
    basis,
    scf.mol_orbs,
    options=NboComputeOptions(include_nra=True),
)

print(nbo.primary_report(results))
print(nbo.nra_report(results, level="summary"))
```

The `compute` method accepts `mode="npa"` and `mode="primary"`. Both names run the same current pipeline; `primary` is retained as a readable alias for callers that expect a Lewis-assignment result.

## Options

`NboComputeOptions` controls candidate generation, Lewis scoring, and NRA/NRT post-processing.

| Option | Default | Meaning |
| --- | ---: | --- |
| `max_bond_order` | `3` | Maximum generated bond multiplicity. |
| `bond_cutoff` | implementation default | Population/coupling threshold for bond candidates. |
| `pi_coupling_cutoff` | implementation default | Minimum pi-coupling signal for pi candidates. |
| `lone_pair_cutoff` | implementation default | Minimum population for lone-pair candidates. |
| `rydberg_max_occupation` | `0.50` | Maximum occupation for one-center `RY` acceptor candidates. |
| `max_alternatives` | implementation default | Maximum Lewis/resonance alternatives retained. |
| `lewis_weight_beta` | `4.0` | Softmax sharpness for score/ranking weights. |
| `include_nra` | `False` | Build NRA/NRT density-fit weights for alternatives. |
| `nra_subspace` | `"selected"` | NAO subspace used in the density fit: `selected`, `pi`, `valence`, or `full`. |
| `nra_fit_metric` | `"frobenius"` | Residual norm reported for the density fit. |
| `nra_prior_weights` | `None` | Optional prior weights by structure index or label. |
| `nra_prior_strength` | `0.0` | Regularization strength for prior-guided NRA/NRT. |
| `nra_prior_mode` | `"regularized"` | Prior handling mode. |
| `nra_spin_fit` | `"total_spin"` | Open-shell NRA/NRT target using total and spin density residuals. |

## Constraints

`NboConstraints` provides optional structure guidance without changing the data model.

```python
from veloxchem.nbodriver import NboConstraints

constraints = NboConstraints(
    required_bonds=((0, 1),),
    forbidden_bonds=((2, 3),),
    required_pi_bonds=((0, 1),),
    allowed_pi_bonds=((0, 1), (2, 3)),
)

results = nbo.compute(molecule, basis, scf.mol_orbs, constraints=constraints)
```

A plain dictionary with the same keys is also accepted. Atom indices are zero-based in the input and one-based in most human-readable reports.

## Result schema

`compute` returns a dictionary. The most important top-level entries are:

| Key | Meaning |
| --- | --- |
| `nao` | NAO transformation, density matrices, atom-block metadata, and diagnostics. |
| `npa` | Natural populations, charges, and per-atom summaries. |
| `mo_nao` | Molecular-orbital composition in the NAO basis. |
| `nbo_candidates` | Generated `CR`, `LP`, `BD`, `SOMO`, `BD*`, and `RY` candidates. |
| `primary_assignment` | Selected Lewis/NBO assignment and score diagnostics. |
| `alternatives` | Lewis/resonance alternatives used for reporting and NRA/NRT fitting. |
| `donor_acceptor_diagnostics` | Occupied-donor to `BD*`/`RY` acceptor density-coupling diagnostics. |
| `nra` | Optional NRA/NRT density-fit weights and residuals. |
| `orbital_analysis` | The shared `OrbitalAnalysisResult` used by NBO and VB; includes AO maps, NAO data, spin data, MO analysis, and candidate records. |

Candidate records include the candidate type, subtype, atom indices, electron count, occupation, vector coefficients in the NAO basis, and chemically relevant metadata such as sigma/pi character or acceptor parentage. These are the same records exposed to `VbDriver` for analyzer-driven orbital recognition.

Primary assignments and alternatives expose the partition used by the resonance analysis:

| Field | Meaning |
| --- | --- |
| `nbo_list` | Complete occupied Lewis/NBO list for the structure. |
| `sigma_nbo_list` | Occupied sigma/core/lone-pair framework. |
| `pi_nbo_list` | Occupied pi candidates. |
| `fixed_nbo_list` | Occupied candidates kept fixed during resonance enumeration. |
| `active_nbo_list` | Occupied candidates in the resonance active space. |
| `active_pi_nbo_list` | Active pi candidates. |
| `active_lone_pair_nbo_list` | Active lone-pair donation candidates. |
| `active_one_electron_nbo_list` | Active one-electron radical candidates. |

The internal representation deliberately uses these mechanical fields rather than named chemical structure classes. Reports use molecule-independent pi-bond signatures such as `pi:1-2,3-4`; these signatures are annotations for display and prior matching, not separate structure types.

## Equations and invariants

The NAO transformation matrix `T` is built to satisfy

```text
T^T S T = I
```

where `S` is the AO overlap matrix. The total NAO density is

```text
D_NAO = T^T S P S T
```

for AO density matrix `P`. For unrestricted references the spin density is

```text
D_spin_NAO = D_alpha_NAO - D_beta_NAO
```

Candidate occupations are evaluated as

```text
n(c) = c^T D_NAO c
m(c) = c^T D_spin_NAO c
```

where `c` is a normalized candidate vector in NAO space. The implementation records diagnostics for orthonormality, density symmetry, trace conservation, electron count, and charge conservation.

## Lewis alternatives and ranking weights

Lewis/resonance alternatives are generated by exact electron-count accounting in the active space. Candidate-only acceptors such as `BD*` and `RY` are not inserted into the occupied Lewis list.

Each alternative has a transparent Lewis score assembled from terms such as bond coverage, formal-charge balance, octet/duet diagnostics, constraint satisfaction, and active-space consistency. The printed alternative `weight` is a score/ranking weight:

```text
w_score,k = exp(-beta s_k) / sum_j exp(-beta s_j)
```

where `s_k` is the Lewis score and `beta = options.lewis_weight_beta`. This weight ranks Lewis alternatives; it is not an NRA/NRT density weight and not a VB wavefunction weight.

## NRA/NRT density-fit weights

When `include_nra=True`, the driver fits the SCF density as a convex combination of alternative-structure densities. For closed-shell systems the fit solves

```text
min_w || d - A w ||^2
subject to w_k >= 0, sum_k w_k = 1
```

where `d` is the target NAO density vector and each column of `A` is an alternative density vector in the requested subspace. For open-shell systems with `nra_spin_fit="total_spin"`, the target concatenates total-density and spin-density residuals.

The fitted weights are reported as `results["nra"]["weights"]` and as `structures[k]["nra_weight"]`. These are density-reconstruction weights. They should be interpreted as how strongly each alternative contributes to reconstructing the SCF density within the chosen subspace, not as wavefunction amplitudes.

## Prior-guided NRA/NRT

Priors can guide the density fit without replacing the density target:

```python
options = NboComputeOptions(
    include_nra=True,
    nra_prior_weights={"pi:1-2,3-4": 0.8, "pi:1-4,2-3": 0.2},
    nra_prior_strength=0.05,
)
```

With regularization, the fitted objective becomes

```text
min_w || d - A w ||^2 + lambda || w - w0 ||^2
subject to w_k >= 0, sum_k w_k = 1
```

The prior vector `w0`, regularization strength, and matched structures are reported in `results["nra"]["prior"]`. Each structure may also carry `prior_weight`. Priors are user guidance, not a separate physical population analysis.

## Donor-acceptor diagnostics

The driver generates `BD*` antibonding complements and one-center `RY` acceptor complements as candidate-only objects. Occupied donors are coupled to these acceptors through a density-coupling diagnostic in the NAO basis. The diagnostic is useful for identifying plausible donation channels, but it is not an NBO second-order perturbation energy because the current driver does not use an NBO Fock matrix or energy denominator.

## Reports

`primary_report(results)` summarizes the selected Lewis/NBO assignment, candidate counts, atom accounting, and donor-acceptor diagnostics.

`nra_report(results, level="summary")` prints compact NRA/NRT weights, residuals, spin-fit metadata, and prior metadata. `level="full"` includes structure details and active-space content.

## Scope notes

The current implementation is a complete structured NBO/NRA/NRT analysis layer for the documented API. Scientific quantities with mature NBO-specific definitions remain clearly separated from diagnostics: donor-acceptor couplings are density diagnostics and NRA/NRT weights are density-fit weights. VB/wavefunction weights are owned by `VbDriver`; the NBO API does not compute or reinterpret them.
