# NboDriver API guide

This document describes the current clean-restart `NboDriver` API in VeloxChem.  The implementation currently provides:

1. natural atomic orbital / natural population analysis (NAO/NPA),
2. molecular-orbital composition in the NAO basis,
3. first-pass NBO candidate generation,
4. a deterministic primary Lewis-like assignment,
5. first-pass Lewis/resonance alternatives over compatible pi-bond matchings,
6. optional first-pass NRA/NRT density-fit weights for closed-shell alternatives.

The present code is intentionally transparent and conservative.  It is not yet a complete replacement for all mature NBO/NRT/VB functionality.

## Public import

The driver is exported through the top-level VeloxChem Python module:

```python
import veloxchem as vlx

nbo = vlx.NboDriver()
```

The public symbols are:

```python
vlx.NboDriver
vlx.NboComputeOptions
vlx.NboConstraints
```

## Main workflow

The normal workflow is:

1. build a `Molecule`,
2. build a `MolecularBasis`,
3. run SCF,
4. call `NboDriver.compute(molecule, basis, scf.mol_orbs)`,
5. print or inspect NPA, MO, and NBO reports,
6. optionally inspect `results["nra"]` when `include_nra=True`.

Mathematically, the driver starts from the AO overlap matrix $S$ and an AO density matrix $P$.  For restricted closed-shell SCF,

$$
P = 2 C_{\alpha} f_{\alpha} C_{\alpha}^{T},
$$

where $C_{\alpha}$ is the alpha MO coefficient matrix and $f_{\alpha}$ is the diagonal occupation matrix stored by VeloxChem.  For unrestricted SCF,

$$
P = C_{\alpha} f_{\alpha} C_{\alpha}^{T}
  + C_{\beta} f_{\beta} C_{\beta}^{T}.
$$

The resulting NAO transformation $T$ is used to obtain an orthonormal NAO density matrix

$$
D^{\mathrm{NAO}} = T^T S P S T,
\qquad
T^T S T = I.
$$

NAO populations are the diagonal elements

$$
n_i = D^{\mathrm{NAO}}_{ii}.
$$

Natural atomic charges are then reported as

$$
q_A = Z_A - \sum_{i \in A} n_i,
$$

where $Z_A$ is the nuclear charge and the sum runs over NAOs assigned to atom $A$.

## `NboDriver.compute()`

Signature:

```python
results = nbo.compute(
    molecule,
    basis,
    mol_orbs,
    mode="npa",
    options=None,
    constraints=None,
)
```

### Arguments

| Argument | Type | Meaning |
| --- | --- | --- |
| `molecule` | `vlx.Molecule` | Molecular geometry, atom labels, charge, multiplicity, connectivity, and electron count. |
| `basis` | `vlx.MolecularBasis` | AO basis used in the SCF calculation. |
| `mol_orbs` | VeloxChem molecular orbitals | SCF molecular orbital object, normally `scf_drv.mol_orbs`. |
| `mode` | `str` | Currently accepts `"npa"` and `"primary"`.  Both run the same current pipeline. |
| `options` | `None`, `dict`, or `NboComputeOptions` | Numerical and reporting options. |
| `constraints` | `None`, `dict`, or `NboConstraints` | User constraints for bond and pi-bond selection. |

### Return value

`compute()` returns a dictionary.  Important keys are:

| Key | Meaning |
| --- | --- |
| `nao_transform` | AO-to-NAO coefficient matrix $T$. |
| `nao_density_matrix` | $D^{\mathrm{NAO}} = T^T S P S T$. |
| `nao_overlap_matrix` | $T^T S T$, expected to be close to the identity. |
| `nao_populations` | Diagonal populations $n_i = D^{\mathrm{NAO}}_{ii}$. |
| `nao_atom_map` | Zero-based atom index assigned to each NAO. |
| `nao_l_map` | Angular momentum label for each NAO, with $s=0$, $p=1$, $d=2$, ... |
| `natural_charges` | Natural charge array $q_A$. |
| `mo_analysis` | MO-in-NAO decomposition data. |
| `nbo_candidates` | Full generated candidate pool. |
| `nbo_list` | Current selected primary NBO list. |
| `primary` | Primary Lewis-like assignment metadata. |
| `alternatives` | Lewis/resonance alternatives and current score/ranking weights. |
| `nra` | Optional Natural Resonance Analysis density-fit weights, present only when `include_nra=True`. |
| `diagnostics` | Electron count, orthonormality, equivalence groups, constraint summary, and other checks. |
| `provenance` | API version, mode, options, and constraints used. |

## Options

Options may be supplied either as a dictionary or as an `NboComputeOptions` instance.

```python
options = {
    "include_diagnostics": True,
    "include_mo_analysis": True,
    "include_nbo_candidates": True,
    "include_lewis_assignment": True,
    "mo_analysis_top": 6,
    "mo_analysis_threshold": 1.0e-2,
    "lone_pair_min_occupation": 1.50,
    "bond_min_occupation": 1.20,
    "bond_min_atom_weight": 0.10,
    "pi_min_occupation": 0.20,
    "conjugated_pi_max_path": 2,
    "max_alternatives": 12,
    "lewis_weight_beta": 4.0,
    "include_nra": False,
    "nra_subspace": "selected",
    "nra_fit_metric": "frobenius",
}
```

The most important numerical definitions are:

$$
n(c) = c^T D^{\mathrm{NAO}} c,
\qquad
\|c\|_2 = 1,
$$

where $n(c)$ is the occupation of a normalized one- or two-center candidate vector $c$ in the NAO basis.

The score-weight sharpness parameter `lewis_weight_beta` enters the current softmax ranking model:

$$
w_k = \frac{\exp\{\beta(s_k - s_{\max})\}}
           {\sum_l \exp\{\beta(s_l - s_{\max})\}},
\qquad
s_{\max} = \max_l s_l.
$$

Larger $\beta$ gives more score/ranking weight to the highest-scoring alternatives.

The optional NRA layer is a post-processing fit of ideal closed-shell Lewis
density matrices from `alternatives` to the actual NAO density. It does not
replace the current score/ranking weights. The first implemented fit metric is
`"frobenius"`. The available subspace choices are:

| `nra_subspace` | Meaning |
| --- | --- |
| `"selected"` | Fit only NAOs used by the selected alternatives. |
| `"pi"` | Fit p-type NAOs on atoms participating in alternative pi bonds. |
| `"valence"` | Fit selected non-core NBO support. |
| `"full"` | Fit the full NAO density matrix. |

When `include_nra=True`, the result contains:

```python
results["nra"] = {
    "subspace": "selected",
    "fit_metric": "frobenius",
    "weights": [...],
    "residual_norm": ...,
    "relative_residual": ...,
    "structures": [...],
    "warnings": [...],
}
```

Each NRA structure reports its `nra_weight`, the corresponding current
`score_weight`, the underlying `score`, its `pi_bonds`, and a single-structure
residual norm. The NRA implementation currently targets closed-shell singlet
Lewis alternatives.

## Constraints

Constraints may be supplied either as a dictionary or as an `NboConstraints` instance.

```python
constraints = {
    "required_bonds": [],
    "forbidden_bonds": [],
    "required_pi_bonds": [],
    "allowed_pi_bonds": [],
    "forbidden_pi_bonds": [],
    "fixed_lone_pairs": [],
    "fixed_core": [],
    "fragment_locks": [],
    "formal_charges": {},
}
```

Atom pairs may be given as one-based or zero-based pairs.  In user-facing examples, one-based indexing is preferred:

```python
constraints = {
    "required_pi_bonds": [(1, 2), (3, 4), (5, 6)],
}
```

The currently active constraints are:

| Constraint | Current effect |
| --- | --- |
| `required_bonds` | Forces available two-center bond candidates for these atom pairs into the primary assignment when slots permit. |
| `forbidden_bonds` | Excludes two-center bond candidates for these atom pairs. |
| `required_pi_bonds` | Requires pi candidates for these atom pairs in resonance alternatives, if candidates exist. |
| `allowed_pi_bonds` | Allows non-sigma pi candidates for these pairs, including through-conjugated or transannular pairs. |
| `forbidden_pi_bonds` | Excludes pi candidates for these atom pairs. |

The remaining constraint fields are reserved hooks for the next scoring/search layer.

## Reports

The driver stores the most recent `molecule` and `results`, so reports can be called directly after `compute()`.

### NPA report

```python
nbo.npa_report(level="summary")
nbo.npa_report(level="standard")
nbo.npa_report(level="full")
```

Report levels:

| Level | Meaning |
| --- | --- |
| `none` | Return no text. |
| `summary` | Natural charges and core/valence/Rydberg populations. |
| `standard` | Summary plus effective natural electron configurations. |
| `full` | NAO occupancies plus all standard output. |

### MO-in-NAO report

The MO transformation is

$$
C^{\mathrm{NAO}} = T^T S C^{\mathrm{AO}}.
$$

For a normalized canonical MO $\psi_m$,

$$
\psi_m = \sum_i C^{\mathrm{NAO}}_{im}\,\chi^{\mathrm{NAO}}_i,
\qquad
\sum_i |C^{\mathrm{NAO}}_{im}|^2 \approx 1.
$$

The atom weight for atom $A$ is

$$
W_{A m} = \sum_{i \in A} |C^{\mathrm{NAO}}_{im}|^2.
$$

Use:

```python
nbo.mo_report(level="occupied")
nbo.mo_report(level="all")
```

### NBO report

```python
nbo.nbo_report(level="summary")
nbo.nbo_report(level="full")
```

The report order is deterministic:

$$
\mathrm{CR} \rightarrow \mathrm{BD}(\sigma) \rightarrow \mathrm{BD}(\pi)
\rightarrow \mathrm{LP} \rightarrow \mathrm{RY} \rightarrow \mathrm{BD^*}/\mathrm{other}.
$$

### NRA/NRT density-fit data

The first NRA/NRT implementation is exposed as structured data in `results["nra"]`, not as a formatted public report yet. Use `include_nra=True` and inspect the returned dictionary directly:

```python
nra = results["nra"]
print(nra["weights"])
print(nra["residual_norm"], nra["relative_residual"])
for structure in nra["structures"]:
    print(structure["nra_weight"], structure["score_weight"], structure["pi_bonds"])
```

The next reporting step is to add `nbo.nra_report(level="summary")` and `nbo.nra_report(level="full")`.

## Example notebook cell: H2C=O

This is a complete notebook-style Python cell for a restricted HF/STO-3G calculation followed by NBO analysis.

```python
import veloxchem as vlx

h2co_xyz = """
C   0.0000000000   0.0000000000   0.0000000000
O   0.0000000000   0.0000000000   1.2167228600
H   1.0639020000   0.0000000000  -0.2816590000
H  -1.0639020000   0.0000000000  -0.2816590000
"""

molecule = vlx.Molecule.read_str(h2co_xyz)
basis = vlx.MolecularBasis.read(molecule, "sto-3g")

scf = vlx.ScfRestrictedDriver()
scf.xcfun_label = "hf"
scf_results = scf.compute(molecule, basis)

nbo = vlx.NboDriver()
nbo.verbose = False
results = nbo.compute(
    molecule,
    basis,
    scf.mol_orbs,
    options={
        "include_diagnostics": True,
        "include_mo_analysis": True,
        "include_nbo_candidates": True,
        "include_lewis_assignment": True,
    },
)

print(nbo.npa_report(level="full", return_text=True))
print(nbo.mo_report(level="occupied", return_text=True))
print(nbo.nbo_report(level="full", return_text=True))

print("Electron count:", results["diagnostics"]["electron_count"])
print("NAO orthonormality error:", results["diagnostics"]["orthonormality_error"])
```

## Example notebook cell: constrained benzene pi structures

This example allows classical Kekulé and Dewar pi patterns.  The non-ring para pairs are included through `allowed_pi_bonds`; without this, non-sigma transannular pi candidates are not selected automatically.

```python
import veloxchem as vlx

benzene_xyz = """
C   1.396792   0.000000   0.000000
C   0.698396   1.209951   0.000000
C  -0.698396   1.209951   0.000000
C  -1.396792   0.000000   0.000000
C  -0.698396  -1.209951   0.000000
C   0.698396  -1.209951   0.000000
H   2.490290   0.000000   0.000000
H   1.245145   2.156660   0.000000
H  -1.245145   2.156660   0.000000
H  -2.490290   0.000000   0.000000
H  -1.245145  -2.156660   0.000000
H   1.245145  -2.156660   0.000000
"""

mol = vlx.Molecule.read_str(benzene_xyz)
bas = vlx.MolecularBasis.read(mol, "sto-3g")

scf = vlx.ScfRestrictedDriver()
scf.xcfun_label = "hf"
_ = scf.compute(mol, bas)

ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
dewar_para_pi_bonds = [(1, 4), (2, 5), (3, 6)]

nbo = vlx.NboDriver()
nbo.verbose = False
res = nbo.compute(
    mol,
    bas,
    scf.mol_orbs,
    constraints={"allowed_pi_bonds": ring_pi_bonds + dewar_para_pi_bonds},
    options={
        "max_alternatives": 24,
        "pi_min_occupation": 0.05,
        "conjugated_pi_max_path": 2,
    },
)

for alt in res["alternatives"]:
    pi_text = ", ".join(f"{i}-{j}" for i, j in alt["pi_bonds"])
    print(
        f"rank={alt['rank']:2d} "
        f"score_weight={100.0 * alt['weight']:7.2f}% "
        f"score={alt['score']:10.5f} "
        f"pi={pi_text}"
    )
```

## Example notebook cell: first NRA/NRT density-fit check

This example uses the same benzene pi-structure pool but requests NRA/NRT density-fit weights in the pi subspace. The resulting `nra_weight` values are density-fit weights; `score_weight` remains the existing Lewis-ranking weight.

```python
import numpy as np
import veloxchem as vlx

# Reuse the benzene geometry, basis, SCF result, and pi-bond lists from the previous example.
nbo = vlx.NboDriver()
nbo.verbose = False
res = nbo.compute(
    mol,
    bas,
    scf.mol_orbs,
    constraints={"allowed_pi_bonds": ring_pi_bonds + dewar_para_pi_bonds},
    options={
        "include_nra": True,
        "nra_subspace": "pi",
        "nra_fit_metric": "frobenius",
        "max_alternatives": 24,
        "pi_min_occupation": 0.05,
        "conjugated_pi_max_path": 2,
    },
)

nra = res["nra"]
weights = np.array(nra["weights"])

assert np.all(weights >= -1.0e-12)
assert abs(float(np.sum(weights)) - 1.0) < 1.0e-10
assert np.isfinite(nra["residual_norm"])
assert np.isfinite(nra["relative_residual"])

print(
    f"subspace={nra['subspace']} metric={nra['fit_metric']} "
    f"relative_residual={nra['relative_residual']:.6e}"
)

for structure in nra["structures"]:
    pi_text = ", ".join(f"{i}-{j}" for i, j in structure["pi_bonds"])
    print(
        f"rank={structure['rank']:2d} "
        f"nra_weight={100.0 * structure['nra_weight']:7.3f}% "
        f"score_weight={100.0 * structure['score_weight']:7.3f}% "
        f"pi={pi_text}"
    )
```

## Interpreting current resonance score weights

The current alternatives are first-pass Lewis-like alternatives, not Natural Resonance Analysis weights and not full VB state mixing.  If an alternative has score $s_k$, its printed score/ranking weight is the softmax model

$$
w_k = \frac{e^{\beta(s_k-s_{\max})}}{\sum_l e^{\beta(s_l-s_{\max})}}.
$$

This is useful for development and ranking, but it should be described as a heuristic score/ranking weight, not as a physical NRA or VB weight.

## Common diagnostics

Useful checks are:

```python
diag = results["diagnostics"]
print(diag["electron_count"])
print(diag["orthonormality_error"])
print(diag["mo_nao_max_normalization_error"])
print(diag["nbo_candidate_counts"])
print(diag["primary_nbo_counts"])
print(diag["constraint_summary"])

if "nra" in results:
    print(results["nra"]["residual_norm"])
    print(results["nra"]["relative_residual"])
```

Expected identities are:

$$
\mathrm{Tr}(D^{\mathrm{NAO}}) \approx N_e,
\qquad
\|T^T S T - I\|_F \ll 1,
\qquad
\sum_A q_A = Q_{\mathrm{molecule}}.
$$
