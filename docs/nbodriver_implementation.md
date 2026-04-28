# NboDriver implementation notes

This document records the current implementation details of the clean-restart `NboDriver`.  It is intended for developers who need to understand, modify, or extend the code.

The central design choice is separation of layers:

$$
\mathrm{AO/SCF}\;\longrightarrow\;\mathrm{NAO/NPA}\;\longrightarrow\;
\mathrm{MO\ analysis}\;\longrightarrow\;\mathrm{NBO\ candidates}\;\longrightarrow\;
\mathrm{Lewis/resonance\ assignment}\;\longrightarrow\;\mathrm{NRA/NRT\ density\ fit}.
$$

Each layer stores structured data in `results`, rather than mixing analysis, assignment, and printing.

## 1. Inputs and AO density

The driver receives a `Molecule`, a `MolecularBasis`, and a VeloxChem molecular-orbital object.  The AO overlap matrix is computed as

$$
S_{\mu\nu} = \langle \chi_\mu | \chi_\nu \rangle.
$$

The total AO density matrix $P$ is built from the SCF orbitals.  For restricted orbitals,

$$
P_{\mu\nu} = 2 \sum_i f_i C_{\mu i} C_{\nu i}.
$$

For unrestricted orbitals,

$$
P_{\mu\nu} = \sum_i f_i^\alpha C^\alpha_{\mu i} C^\alpha_{\nu i}
           + \sum_i f_i^\beta  C^\beta_{\mu i}  C^\beta_{\nu i}.
$$

The electron-count invariant in a nonorthogonal AO basis is

$$
N_e = \mathrm{Tr}(SP).
$$

This invariant is used to guard conservative density symmetrization.

## 2. AO metadata

The function `_ao_shell_map()` builds two integer arrays:

$$
a(\mu) = \text{atom index of AO } \mu,
\qquad
\ell(\mu) = \text{angular momentum of AO } \mu.
$$

The angular momentum encoding is

$$
s=0,\quad p=1,\quad d=2,\quad f=3,\quad g=4,\quad h=5.
$$

These maps are deliberately kept simple.  They are used to maintain atom labels during NAO construction and to identify p-type spaces for pi candidates.

## 3. Equivalent atoms and local frames

Equivalent atom groups are detected from element labels and local connectivity/distance signatures.  For atom $A$, a simplified signature is

$$
\sigma_A = \left(Z_A, \mathrm{sort}\{(Z_B, \lfloor r_{AB}/\tau \rceil): B \in \mathcal{N}(A)\}\right),
$$

where $\tau$ is the distance tolerance and $\mathcal{N}(A)$ is the connectivity-neighbor set.

Local frames are also built for each atom.  If atom $A$ has neighbor directions $\hat{b}_1, \hat{b}_2, \ldots$, the first axis is

$$
\hat{x}_A = \hat{b}_1.
$$

If a second bond direction exists, the second axis is the component of $\hat{b}_2$ orthogonal to $\hat{x}_A$:

$$
\hat{y}_A = \frac{\hat{b}_2 - (\hat{b}_2 \cdot \hat{x}_A)\hat{x}_A}
                 {\|\hat{b}_2 - (\hat{b}_2 \cdot \hat{x}_A)\hat{x}_A\|}.
$$

The third axis is

$$
\hat{z}_A = \frac{\hat{x}_A \times \hat{y}_A}{\|\hat{x}_A \times \hat{y}_A\|},
\qquad
\hat{y}_A \leftarrow \frac{\hat{z}_A \times \hat{x}_A}{\|\hat{z}_A \times \hat{x}_A\|}.
$$

At present these frames are diagnostic and a future hook for oriented NHO construction.

## 4. Conservative equivalent-atom density symmetrization

The current symmetrizer only acts on equivalent groups with identical all-s AO layouts.  This is intentionally conservative: it stabilizes equivalent hydrogens without rotating p/d functions between atoms.

For an equivalent group with $g$ atoms, AO permutation matrices $\Pi_k$ are built from cyclic atom permutations.  The symmetrized density is

$$
P_{\mathrm{sym}} = \frac{1}{g}\sum_{k=1}^{g} \Pi_k^T P \Pi_k.
$$

The result is made symmetric:

$$
P_{\mathrm{sym}} \leftarrow \frac{1}{2}(P_{\mathrm{sym}} + P_{\mathrm{sym}}^T).
$$

The implementation checks electron conservation:

$$
\left|\mathrm{Tr}(S P_{\mathrm{sym}}) - \mathrm{Tr}(SP)\right| < \epsilon.
$$

If the check fails, the original density is retained.  This avoids the earlier electron-loss bug caused by row/column averaging in a nonorthogonal basis.

## 5. Orthonormal NAO construction

The current NAO construction is intentionally simple and robust.

### 5.1 Löwdin orthogonalized AO basis

Diagonalize the AO overlap matrix:

$$
S = U s U^T.
$$

The square Löwdin inverse square root is

$$
X = S^{-1/2} = U s^{-1/2} U^T.
$$

It satisfies

$$
X^T S X = I.
$$

The density in this orthonormal AO-like basis is

$$
D^{(X)} = X^T S P S X.
$$

The symmetrized numerical form is

$$
D^{(X)} \leftarrow \frac{1}{2}\left(D^{(X)} + D^{(X)T}\right).
$$

### 5.2 Atom-block density diagonalization

For each atom $A$, collect the Löwdin AO indices

$$
\mathcal{I}_A = \{\mu : a(\mu)=A\}.
$$

The atom block is

$$
D_A = D^{(X)}[\mathcal{I}_A,\mathcal{I}_A].
$$

Then solve

$$
D_A R_A = R_A n_A.
$$

The eigenvectors are sorted by descending occupation.  Near-degenerate subspaces are canonicalized deterministically by dominant angular momentum and AO index.  The block rotations are assembled into a block-diagonal matrix

$$
R = \bigoplus_A R_A.
$$

The final AO-to-NAO transformation is

$$
T = X R.
$$

This guarantees

$$
T^T S T = R^T X^T S X R = R^T R = I.
$$

The final NAO density is

$$
D^{\mathrm{NAO}} = T^T S P S T.
$$

The population of NAO $i$ is

$$
n_i = D^{\mathrm{NAO}}_{ii}.
$$

This atom-block approach prevents the label stealing observed when the full basis is diagonalized globally.

## 6. NPA classification and charges

Core-orbital counts are currently simple closed-shell counts from nuclear charge:

$$
n_{\mathrm{core}}(Z) =
\begin{cases}
0, & Z \le 2,\\
1, & 2 < Z \le 10,\\
5, & 10 < Z \le 18,\\
9, & 18 < Z \le 36,\\
18, & Z > 36.
\end{cases}
$$

For first-row atoms, the core orbital is chosen preferentially from the highest-populated s-type NAO to avoid unphysical labels such as `Cor(1p)`.

The atom population is

$$
N_A = \sum_{i \in A} n_i.
$$

The natural charge is

$$
q_A = Z_A - N_A.
$$

The molecular charge sum should obey

$$
\sum_A q_A = \sum_A Z_A - \sum_A N_A
          = Z_{\mathrm{tot}} - N_e
          = Q_{\mathrm{mol}}.
$$

Valence electron configurations are reported by angular momentum:

$$
N_{A\ell}^{\mathrm{val}} = \sum_{i \in A,\; \ell_i=\ell,\; i\notin\mathrm{core}} n_i.
$$

## 7. MO-in-NAO analysis

Canonical MOs are transformed to the NAO basis as

$$
C^{\mathrm{NAO}} = T^T S C^{\mathrm{AO}}.
$$

For MO $m$, the NAO coefficient vector is $c_m^{\mathrm{NAO}}$.  The normalization diagnostic is

$$
\delta_m = \left|\sum_i |c^{\mathrm{NAO}}_{im}|^2 - 1\right|.
$$

The reported maximum normalization error is

$$
\delta_{\max} = \max_m \delta_m.
$$

The atom contribution to MO $m$ is

$$
W_{Am} = \sum_{i\in A}|c^{\mathrm{NAO}}_{im}|^2.
$$

The report keeps the largest NAO and atom weights according to `mo_analysis_top` and `mo_analysis_threshold`.

## 8. NBO candidate generation

Candidate generation is intentionally separate from final assignment.  Every candidate is stored as a normalized vector $v$ in the NAO basis:

$$
\|v\|_2 = 1.
$$

Its occupation is

$$
n(v) = v^T D^{\mathrm{NAO}} v.
$$

### 8.1 Core candidates

Core candidates are one-center NAOs:

$$
v_i = e_i,
\qquad
n(v_i)=D^{\mathrm{NAO}}_{ii}.
$$

They are stored as type `CR` and subtype `core`.

### 8.2 Lone-pair candidates

For non-hydrogen atoms, non-core valence NAOs with

$$
n_i \ge n_{\mathrm{LP,min}}
$$

are stored as one-center lone-pair candidates.  The current default is

$$
n_{\mathrm{LP,min}} = 1.50.
$$

These are stored as type `LP` and subtype `lone-pair`.

### 8.3 Connected atom-pair candidates

For connected atom pairs $A-B$, the candidate builder removes selected core and obvious one-center lone-pair NAOs from the bonding search space.  The remaining atom-local index sets are

$$
\mathcal{L}_A,\qquad \mathcal{L}_B.
$$

The pair index set is

$$
\mathcal{P}_{AB} = \mathcal{L}_A \cup \mathcal{L}_B.
$$

The pair density block is

$$
D_{AB} = D^{\mathrm{NAO}}[\mathcal{P}_{AB},\mathcal{P}_{AB}].
$$

The algorithm solves

$$
D_{AB} u_k = n_k u_k.
$$

A two-center candidate is accepted if

$$
n_k \ge n_{\mathrm{BD,min}},
$$

and both atomic side weights satisfy

$$
w_A = \sum_{i\in\mathcal{L}_A}|u_{ki}|^2 \ge w_{\mathrm{atom,min}},
\qquad
w_B = \sum_{i\in\mathcal{L}_B}|u_{ki}|^2 \ge w_{\mathrm{atom,min}}.
$$

Defaults are

$$
n_{\mathrm{BD,min}} = 1.20,
\qquad
w_{\mathrm{atom,min}} = 0.10.
$$

For each connected atom pair, the first accepted candidate is labeled

$$
\mathrm{BD}(\sigma),
$$

and later accepted candidates for the same pair are labeled

$$
\mathrm{BD}(\pi).
$$

### 8.4 Non-sigma pi candidates

Non-sigma pi candidates are generated from p-type NAOs.  An atom is pi-capable if

$$
\exists i\in A \quad \ell_i = 1.
$$

For a non-connected pair $A-B$, a p-space candidate may be generated if either:

1. the graph distance $d(A,B)$ is within `conjugated_pi_max_path`, or
2. the pair is explicitly requested by `allowed_pi_bonds` or `required_pi_bonds`.

The graph-distance condition is

$$
d(A,B) \le d_{\max}.
$$

The p-space block is diagonalized exactly like a connected pair block, but only the top accepted p-p candidate is retained.  It is stored as

$$
\mathrm{BD}(\pi),\quad \mathrm{non\_sigma}=\mathrm{True}.
$$

The current primary assignment refuses non-sigma pi candidates unless the pair is explicitly allowed or required by user constraints.  This protects normal Lewis assignments while permitting allyl terminal and Dewar-style tests.

## 9. Primary Lewis-like assignment

The electron-pair target is

$$
N_{\mathrm{pair}} = \mathrm{round}\left(\frac{N_e}{2}\right).
$$

The selection order is deterministic:

1. all core candidates,
2. required bond candidates,
3. required pi-bond candidates,
4. remaining lone-pair and bond candidates by score/order until the target is reached.

The candidate sort key is effectively

$$
\mathrm{CR} \prec \mathrm{LP} \prec \mathrm{BD}(\sigma) \prec \mathrm{BD}(\pi),
$$

with higher occupation before lower occupation inside comparable classes.

The primary occupation sum is

$$
O = \sum_{v\in\mathcal{S}} n(v),
$$

where $\mathcal{S}$ is the selected candidate set.  The current primary score is

$$
s_{\mathrm{primary}} = O - 2\left| |\mathcal{S}| - N_{\mathrm{pair}} \right|.
$$

This is a deliberately simple first-pass objective.

## 10. Lewis/resonance alternatives

The current alternative enumerator keeps the non-pi part of the primary assignment fixed and varies compatible pi-bond matchings.

Let the fixed set be

$$
\mathcal{F} = \{v\in\mathcal{S}_{\mathrm{primary}} : v \text{ is not } \mathrm{BD}(\pi)\}.
$$

The number of pi slots is

$$
N_\pi = N_{\mathrm{pair}} - |\mathcal{F}|.
$$

The pi pool $\mathcal{P}$ contains `BD(pi)` candidates passing forbidden/allowed constraints.  The enumerator chooses combinations

$$
\mathcal{C}_k \subset \mathcal{P},
\qquad
|\mathcal{C}_k| = N_\pi.
$$

A combination is compatible if it is a simple matching:

$$
\forall (A,B),(C,D)\in\mathcal{C}_k,
\quad
\{A,B\}\cap\{C,D\}=\varnothing
\quad\text{unless the two pairs are identical, which is disallowed.}
$$

The candidate selected set for alternative $k$ is

$$
\mathcal{S}_k = \mathcal{F}\cup\mathcal{C}_k.
$$

Its occupation sum is

$$
O_k = \sum_{v\in\mathcal{S}_k} n(v).
$$

The current alternative score is

$$
s_k = O_k + 0.05\,N^{\mathrm{allowed}}_k
      - 2\left||\mathcal{S}_k|-N_{\mathrm{pair}}\right|,
$$

where $N^{\mathrm{allowed}}_k$ is the number of selected pi pairs also present in `allowed_pi_bonds`.  The small allowed-pair bonus is a practical tie-breaker for user-enabled alternatives.

The alternatives are sorted by decreasing $s_k$, truncated to `max_alternatives`, and converted to weights:

$$
w_k = \frac{\exp\{\beta(s_k-s_{\max})\}}
           {\sum_l\exp\{\beta(s_l-s_{\max})\}},
\qquad
s_{\max}=\max_l s_l.
$$

The weights satisfy

$$
\sum_k w_k = 1.
$$

These are development weights for ranking Lewis alternatives, not final VB amplitudes.

## 11. NRA/NRT density-fit layer

The current NRA/NRT layer is optional and runs only when `include_nra=True`. It is implemented as post-processing on the already enumerated `alternatives` list, so it does not affect candidate generation, primary assignment, or the existing score/ranking weights.

For each alternative $k$, the selected NBO coefficient lists are reconstructed into normalized NAO-space vectors $v_{ki}$. The first implementation assumes closed-shell pair occupation for all selected candidates:

$$
D_k^{\mathrm{Lewis}} = 2\sum_i v_{ki}v_{ki}^T.
$$

The target is the actual NAO density $D^{\mathrm{NAO}}$. The fit chooses nonnegative weights that sum to one:

$$
\min_w \left\|D^{\mathrm{NAO}} - \sum_k w_kD_k^{\mathrm{Lewis}}\right\|_F^2,
\qquad
w_k \ge 0,
\qquad
\sum_k w_k = 1.
$$

The dependency-free active-set solver first solves the equality-constrained least-squares problem, removes structures with negative provisional weights, and refits on the remaining active set until all weights are nonnegative. A one-structure pool returns weight 1.0 by construction.

The implemented subspaces are:

| `nra_subspace` | Current behavior |
| --- | --- |
| `"selected"` | Uses the union of NAOs appearing in the selected alternatives. |
| `"pi"` | Uses p-type NAOs on atoms involved in alternative pi bonds; falls back to `"selected"` if empty. |
| `"valence"` | Uses selected non-core candidate support. |
| `"full"` | Uses the full NAO density. |

For the Frobenius metric, symmetric matrix blocks are vectorized over the upper triangle with off-diagonal entries scaled by $\sqrt{2}$, so the vector norm matches the full symmetric-matrix Frobenius norm.

The `results["nra"]` dictionary contains:

```python
{
      "subspace": "pi",
      "fit_metric": "frobenius",
      "weights": [...],
      "residual_norm": ...,
      "relative_residual": ...,
      "structures": [
            {
                  "rank": ...,
                  "pi_bonds": [...],
                  "nra_weight": ...,
                  "score_weight": ...,
                  "score": ...,
                  "residual_norm": ...,
            },
      ],
      "warnings": [...],
}
```

The layer currently targets closed-shell singlet alternatives. Open-shell/radical NRA is intentionally deferred until one-electron candidates and spin-resolved assignment are available.

## 12. Report formatting

The NBO report separates candidate counts from primary counts.  It sorts selected NBOs by chemical reporting priority:

$$
\mathrm{CR},\quad \mathrm{BD}(\sigma),\quad \mathrm{BD}(\pi),\quad
\mathrm{LP},\quad \mathrm{RY},\quad \mathrm{BD^*}/\mathrm{other}.
$$

For a bond candidate between atoms $A$ and $B$, atom labels are reported heavier-to-lighter.  This keeps reports visually stable, for example reporting C-H rather than H-C.

There is not yet a public `nra_report()` formatter. NRA/NRT results are inspected directly through `results["nra"]`.

## 13. Important invariants

The following identities are expected for a healthy run:

### NAO orthonormality

$$
\|T^T S T - I\|_F \ll 1.
$$

### Electron conservation

$$
\mathrm{Tr}(D^{\mathrm{NAO}})
= \mathrm{Tr}(T^T S P S T)
= \mathrm{Tr}(P S T T^T S)
\approx N_e.
$$

Because $T$ spans the full AO space in the current square construction,

$$
T T^T S \approx S^{-1}S = I,
$$

so the NAO trace should match the AO electron count.

### MO normalization

For each canonical MO,

$$
\sum_i |C^{\mathrm{NAO}}_{im}|^2 \approx 1.
$$

### Charge conservation

$$
\sum_A q_A = Q_{\mathrm{mol}}.
$$

### NRA/NRT simplex weights

When NRA is requested for a closed-shell structure pool,

$$
w_k^{\mathrm{NRA}} \ge 0,
\qquad
\sum_k w_k^{\mathrm{NRA}} = 1,
$$

and the residual norms should be finite.

## 14. Current limitations

The implementation is a clean development scaffold, not the final scientific endpoint.  Known limitations include:

1. The primary Lewis score is occupation-driven and simple.
2. Formal charges are stored in `NboConstraints` but are not yet deeply coupled into scoring.
3. `fixed_lone_pairs`, `fixed_core`, and `fragment_locks` are reserved hooks.
4. Radical/open-shell alternatives are still mostly pair-based diagnostics.
5. `alternatives[*]["weight"]` values are softmax ranking weights, not NRA or VB weights.
6. NRA/NRT density weights are implemented only as a first closed-shell density fit; there is no spin-resolved radical NRA yet.
7. There is no public `nra_report()` formatter yet.
8. Antibonding NBOs and Rydberg complements are not yet constructed as a full NBO basis.
9. NHO/hybrid directionality is not yet a separate layer; current candidates are density-block natural orbitals in NAO space.

## 15. Suggested next implementation steps

The natural next steps are:

1. Add an explicit valence/electron-accounting model per atom.
2. Add octet and hypervalence penalties:

$$
E_{\mathrm{octet}} = \sum_A \lambda_A \max(0, N_A^{\mathrm{val}} - N_A^{\mathrm{allowed}})^2.
$$

3. Activate formal-charge scoring:

$$
E_{\mathrm{charge}} = \sum_A \gamma_A(q_A^{\mathrm{Lewis}} - q_A^{\mathrm{target}})^2.
$$

4. Add radical-specific one-electron alternatives, where a singly occupied candidate has occupancy target

$$
n(v) \approx 1
$$

instead of pair occupancy

$$
n(v) \approx 2.
$$

5. Add `nra_report(level="summary")` and `nra_report(level="full")` formatters.
6. Expand NRA validation to carboxylate, nitrate, ozone, nitrobenzene, and amides.
7. Separate sigma framework selection from pi-resonance enumeration.
8. Add optional prior-regularized NRA weights.
9. Replace or supplement the current softmax score with a documented VB/BOND-inspired model when available.

## 16. Developer checklist

When changing the implementation, verify at least:

$$
\mathrm{Tr}(D^{\mathrm{NAO}}) \approx N_e,
\qquad
\|T^T S T - I\|_F \ll 1,
\qquad
\max_m\left|\sum_i |C^{\mathrm{NAO}}_{im}|^2 - 1\right| \ll 1.
$$

Also inspect:

```python
results["diagnostics"]["nbo_candidate_counts"]
results["diagnostics"]["primary_nbo_counts"]
results["alternatives"]
```

When `include_nra=True`, also verify:

```python
weights = np.array(results["nra"]["weights"])
assert np.all(weights >= -1.0e-12)
assert abs(float(np.sum(weights)) - 1.0) < 1.0e-10
assert np.isfinite(results["nra"]["residual_norm"])
assert np.isfinite(results["nra"]["relative_residual"])
```

For benchmark-style development examples, keep checking H2C=O, water, methane, ethylene, benzene, allyl cation, allyl radical, allyl anion, and constrained benzene Kekulé/Dewar alternatives.
