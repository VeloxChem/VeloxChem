# NboDriver improvement plan and Natural Resonance Analysis roadmap

This document collects proposed improvements for the clean-restart `NboDriver`, with special focus on a Natural Resonance Analysis (NRA) layer.

The NRA API elements described here are partly implemented. The current `NboDriver` API supports `include_nra`, `nra_subspace`, `nra_fit_metric`, and `results["nra"]` as a first closed-shell density-fitting layer. Higher-level interfaces such as `nra_report()` and chemically richer NRA/NRT options remain proposed next-step interfaces, not implemented API.

Current implementation status:

| Area | Status |
| --- | --- |
| NAO/NPA invariants | Implemented and covered by automated tests. |
| NBO count sanity checks | Implemented for water, methane, ethylene, and benzene. |
| Bond and pi-bond constraints | Implemented for required/forbidden bonds and required/allowed/forbidden pi bonds. |
| Lewis/resonance alternatives | Implemented as compatible pi-bond matchings with score/ranking weights. |
| NRA/NRT density weights | Implemented as optional closed-shell nonnegative density fitting in selected/pi/valence/full subspaces. |
| NRA reporting | Not implemented; inspect `results["nra"]` directly. |
| Radical/open-shell NRA | Not implemented; current radical behavior is diagnostic only. |

The guiding principle is to preserve the current clean architecture:

$$
\mathrm{AO/SCF}
\rightarrow
\mathrm{NAO/NPA}
\rightarrow
\mathrm{MO\ analysis}
\rightarrow
\mathrm{NBO\ candidates}
\rightarrow
\mathrm{Lewis/resonance\ assignment}
\rightarrow
\mathrm{NRA}.
$$

The next stages should remain modular.  NRA should not be folded into candidate generation or report formatting; it should be a separate post-processing layer built on top of the resonance-structure alternatives.

## 1. Stabilize the current foundation

The current implementation already has useful separation of concerns:

1. NAO/NPA construction,
2. MO-in-NAO analysis,
3. NBO candidate generation,
4. primary Lewis-like assignment,
5. resonance alternatives over compatible pi matchings.

Before adding more sophisticated features, preserve these invariants:

$$
T^T S T \approx I,
$$

$$
D^{\mathrm{NAO}} = T^T S P S T,
$$

$$
\mathrm{Tr}\left(D^{\mathrm{NAO}}\right) \approx N_e,
$$

$$
\sum_A q_A = Q_{\mathrm{mol}}.
$$

Recommended priorities:

1. Keep NAO/NPA conservative and deterministic.
2. Make candidate generation broader but transparent.
3. Improve Lewis assignment scoring.
4. Keep NRA as a separate density-fitting layer and expand it incrementally.

## 2. Immediate NBO improvements

### 2.1 Improve Lewis scoring

The current primary score is mostly occupation-based:

$$
s = \sum_{v \in \mathcal{S}} n(v)
    - 2\left||\mathcal{S}| - N_{\mathrm{pair}}\right|.
$$

This should be replaced by a chemically informed score:

$$
s =
w_{\mathrm{occ}} S_{\mathrm{occ}}
- w_{\mathrm{octet}} P_{\mathrm{octet}}
- w_{\mathrm{charge}} P_{\mathrm{charge}}
- w_{\mathrm{valence}} P_{\mathrm{valence}}
- w_{\mathrm{nonclassical}} P_{\mathrm{nonclassical}}.
$$

The occupation term is

$$
S_{\mathrm{occ}} = \sum_{v \in \mathcal{S}} n(v),
\qquad
n(v) = v^T D^{\mathrm{NAO}} v.
$$

An octet penalty could be

$$
P_{\mathrm{octet}} =
\sum_A
\max\left(0, N_A^{\mathrm{Lewis}} - N_A^{\mathrm{allowed}}\right)^2.
$$

A formal-charge penalty could be

$$
P_{\mathrm{charge}} =
\sum_A
\left(q_A^{\mathrm{Lewis}} - q_A^{\mathrm{target}}\right)^2.
$$

This would activate currently reserved fields such as `formal_charges` in `NboConstraints`.

### 2.2 Separate sigma framework and pi resonance

The selected Lewis set should be decomposed into chemically meaningful parts:

$$
\mathcal{S}
=
\mathcal{S}_{\mathrm{CR}}
\cup
\mathcal{S}_{\sigma}
\cup
\mathcal{S}_{\pi}
\cup
\mathcal{S}_{\mathrm{LP}}.
$$

First determine a stable sigma framework:

$$
\mathcal{S}_{\sigma}
=
\arg\max_{\sigma} s_{\sigma}.
$$

Then enumerate resonance structures in the pi/lone-pair space:

$$
\mathcal{S}_k
=
\mathcal{S}_{\mathrm{CR}}
\cup
\mathcal{S}_{\sigma}
\cup
\mathcal{S}_{\mathrm{LP},k}
\cup
\mathcal{S}_{\pi,k}.
$$

This is important for benzene, allyl systems, carboxylates, nitro groups, nitrate, ozone, carbonyls, and related conjugated systems.

### 2.3 Improve radical and open-shell handling

The current Lewis enumeration is pair-based.  Radical systems require one-electron candidates and spin-resolved density analysis.

Useful future candidate types:

1. `SOMO`,
2. one-electron `BD(pi)`,
3. one-electron `LP`,
4. spin-resolved radical candidates.

For a singly occupied candidate,

$$
n(v) \approx 1,
$$

whereas a closed-shell pair has

$$
n(v) \approx 2.
$$

For spin-resolved analysis, use

$$
n_{\alpha}(v) = v^T D^{\alpha} v,
$$

$$
n_{\beta}(v) = v^T D^{\beta} v,
$$

and spin population

$$
m(v) = n_{\alpha}(v) - n_{\beta}(v).
$$

This would make allyl radical analysis meaningful rather than only diagnostic.

### 2.4 Add antibonding and Rydberg complements

The current candidate layer focuses on occupied candidates.  A more complete NBO implementation should generate antibonding partners.

For a bonding candidate $v_{\mathrm{BD}}$, construct a same-subspace orthogonal complement $v_{\mathrm{BD^*}}$ such that

$$
v_{\mathrm{BD}}^T v_{\mathrm{BD^*}} = 0.
$$

Its occupation is

$$
n_{\mathrm{BD^*}} =
v_{\mathrm{BD^*}}^T D^{\mathrm{NAO}} v_{\mathrm{BD^*}}.
$$

This enables donor-acceptor analysis of the form

$$
\Delta E_{i \rightarrow j^*}
\approx
q_i
\frac{|F_{ij^*}|^2}{\epsilon_{j^*} - \epsilon_i}.
$$

This should be a later stage, after the Lewis and NRA layers are stable.

## 3. Natural Resonance Analysis concept

The current resonance alternatives use a heuristic score-softmax model:

$$
w_k^{\mathrm{score}}
=
\frac{\exp\{\beta(s_k-s_{\max})\}}
     {\sum_l \exp\{\beta(s_l-s_{\max})\}}.
$$

These are ranking weights, not true Natural Resonance Analysis weights.

A real NRA-like layer should answer:

> Given a set of resonance structures, what nonnegative weights best reproduce the actual NAO density?

Let each resonance structure $k$ define an idealized Lewis density matrix

$$
D_k^{\mathrm{Lewis}}.
$$

The actual density is

$$
D^{\mathrm{NAO}}.
$$

NRA should find weights $w_k$ such that

$$
D^{\mathrm{model}}
=
\sum_k w_k D_k^{\mathrm{Lewis}},
$$

with constraints

$$
w_k \ge 0,
\qquad
\sum_k w_k = 1.
$$

The basic objective is

$$
\chi^2
=
\left\|
D^{\mathrm{NAO}}
-
\sum_k w_k D_k^{\mathrm{Lewis}}
\right\|_F^2.
$$

A chemically weighted objective is

$$
\chi^2
=
\sum_{ij} W_{ij}
\left(
D_{ij}^{\mathrm{NAO}}
-
\sum_k w_k D_{k,ij}^{\mathrm{Lewis}}
\right)^2,
$$

where $W_{ij}$ can emphasize valence, pi, or lone-pair subspaces.

## 4. Building ideal Lewis densities

For resonance structure $k$, let selected occupied NBO vectors be

$$
\{v_{k1}, v_{k2}, \ldots, v_{kN}\}.
$$

For a closed-shell structure,

$$
D_k^{\mathrm{Lewis}}
=
2\sum_p v_{kp}v_{kp}^T.
$$

For open-shell structures,

$$
D_k^{\mathrm{Lewis}}
=
2\sum_{p\in\mathrm{pairs}} v_{kp}v_{kp}^T
+
\sum_{r\in\mathrm{radicals}} v_{kr}v_{kr}^T.
$$

The first implementation can be closed-shell only.  Radical support can be added after spin-resolved candidate generation is available.

## 5. NRA fitting formulation

Vectorize the target density:

$$
d = \mathrm{vec}\left(D^{\mathrm{NAO}}\right).
$$

Build a matrix of vectorized Lewis densities:

$$
A =
\begin{bmatrix}
\mathrm{vec}(D_1^{\mathrm{Lewis}}) &
\mathrm{vec}(D_2^{\mathrm{Lewis}}) &
\cdots &
\mathrm{vec}(D_M^{\mathrm{Lewis}})
\end{bmatrix}.
$$

Then solve the constrained least-squares problem

$$
\min_w \|d - Aw\|_2^2,
$$

subject to

$$
w_k \ge 0,
\qquad
\sum_k w_k = 1.
$$

The result is a density-reconstruction weight:

$$
w_k^{\mathrm{NRA}}
=
\arg\min_{w_k \ge 0,\;\sum_k w_k=1}
\left\|
D^{\mathrm{NAO}}
-
\sum_k w_kD_k^{\mathrm{Lewis}}
\right\|^2.
$$

This should be reported separately from the current score-softmax weights.

## 6. NRA weights versus other weights

The implementation should distinguish three different meanings:

| Weight type | Symbol | Meaning |
| --- | --- | --- |
| Lewis ranking weight | $w_k^{\mathrm{score}}$ | Current heuristic softmax from Lewis score. |
| NRA density weight | $w_k^{\mathrm{NRA}}$ | Weight from density reconstruction. |
| VB/BOND wavefunction weight | $w_k^{\mathrm{VB}}$ | Future wavefunction/state-mixing contribution. |

The current `alternatives[*]["weight"]` should eventually be renamed or documented clearly as a ranking weight, not an NRA weight.

## 7. NRA implementation stages

### Stage NRA-1: structure pool implemented

Use current `alternatives` as the first resonance-structure pool:

$$
\mathcal{R} = \{R_1, R_2, \ldots, R_M\}.
$$

Each alternative already contains:

1. selected NBO list,
2. pi bonds,
3. score,
4. score-based weight,
5. electron-pair count,
6. occupation sum.

The implemented results section is:

```python
results["nra"] = {
    "subspace": "pi",
    "fit_metric": "frobenius",
    "weights": [...],
    "residual_norm": ...,
    "relative_residual": ...,
    "structures": [...],
    "warnings": [...],
}
```

### Stage NRA-2: ideal Lewis densities implemented for closed-shell alternatives

For each alternative, construct

$$
D_k^{\mathrm{Lewis}}
=
\sum_i o_i v_i v_i^T,
$$

where

$$
o_i =
\begin{cases}
2, & \text{closed-shell pair},\\
1, & \text{radical/SOMO}.
\end{cases}
$$

The current version uses $o_i=2$ for all selected pair candidates. Open-shell/SOMO support is still future work.

### Stage NRA-3: constrained least squares implemented

The implemented dependency-free constrained least-squares solver follows this approach:

1. solve equality-constrained least squares with $\sum_k w_k=1$,
2. set negative weights to zero,
3. refit on the active positive set,
4. iterate until all weights are nonnegative,
5. normalize the final weights.

This is sufficient for small resonance pools.

### Stage NRA-4: chemically restricted fitting implemented at first-pass level

Full-density fitting can be dominated by core and sigma density.  For resonance, pi-subspace fitting is often more meaningful.

Valence subspace:

$$
D^{\mathrm{val}}
=
P_{\mathrm{val}}^T D^{\mathrm{NAO}} P_{\mathrm{val}}.
$$

Pi subspace:

$$
D^{\pi}
=
P_{\pi}^T D^{\mathrm{NAO}} P_{\pi}.
$$

Recommended first option:

```python
options={
    "include_nra": True,
    "nra_subspace": "pi",
    "nra_fit_metric": "frobenius",
}
```

Useful subspace choices:

| `nra_subspace` | Meaning |
| --- | --- |
| `"full"` | Fit full NAO density. |
| `"valence"` | Fit non-core valence NAO density. |
| `"pi"` | Fit p-type/pi-capable subspace. |
| `"selected"` | Fit only NAOs used by the resonance alternatives. |

The current Frobenius implementation vectorizes symmetric blocks with $\sqrt{2}$ off-diagonal scaling so the vector norm matches the matrix Frobenius norm.

### Stage NRA-5: reporting not implemented

Add a public report method:

```python
nbo.nra_report(level="summary")
nbo.nra_report(level="full")
```

The report should show:

1. structure label,
2. pi bonds,
3. NRA weight,
4. current score/ranking weight,
5. residual contribution,
6. formal-charge pattern,
7. dominant bonds and lone pairs.

Example table:

```text
Rank  NRA weight  Score weight  Residual   Pi bonds
  1     49.8%        50.1%       1.2e-03   1-2, 3-4, 5-6
  2     49.8%        49.9%       1.2e-03   2-3, 4-5, 6-1
  3      0.2%         0.0%       8.5e-02   1-4, 2-3, 5-6
```

## 8. User-specified priors

The user may want to supply expected resonance weights.  These should be treated as priors, not as fixed answers, unless explicitly requested.

A prior-regularized objective could be

$$
\chi^2_{\mathrm{prior}}
=
\left\|
D - \sum_k w_kD_k
\right\|^2
+
\lambda\sum_k (w_k-w_k^0)^2,
$$

where $w_k^0$ are user-provided prior weights.

Constraints remain

$$
w_k \ge 0,
\qquad
\sum_k w_k = 1.
$$

This allows user guidance while still making the final weights density-aware.

## 9. Recommended roadmap

### Short term

Completed:

1. Clearly label current alternative weights as score/ranking weights.
2. Add `include_nra`, `nra_subspace`, and `nra_fit_metric` options.
3. Add `results["nra"]` for optional closed-shell density-fit NRA/NRT weights.
4. Build ideal closed-shell Lewis densities from current alternatives.
5. Fit constrained nonnegative weights in selected/pi/valence/full subspaces.
6. Add automated regression tests for NAO invariants, NBO counts, constraints, and benzene NRA weights.

Still short term:

1. Add `nra_report()` summary/full formatters.
2. Add active formal-charge and octet penalties.
3. Improve radical/open-shell handling.
4. Expand NRA regression examples beyond benzene.

### Medium term

1. Add `nra_report()`.
2. Add automatic labels for common resonance patterns:
   - Kekulé,
   - Dewar,
   - allyl terminal,
   - zwitterionic,
   - charge-separated.
3. Add user-specified resonance priors.
4. Activate `formal_charges`, `fixed_lone_pairs`, `fixed_core`, and `fragment_locks` in scoring/search.
5. Expand examples to carboxylate, nitrate, nitrobenzene, ozone, and amides.

### Long term

1. Add complete natural hybrid orbital construction.
2. Add antibonding NBOs and Rydberg complements.
3. Add donor-acceptor second-order perturbation analysis.
4. Add true VB/BOND-style state coupling if desired.
5. Keep separate:
   - Lewis ranking weights,
   - NRA density weights,
   - VB wavefunction weights.

## 10. Best next concrete implementation task

The best next implementation task is:

> Add a public NRA report and expand closed-shell NRA validation beyond benzene.

Good next targets are closed-shell pi systems where the resonance pool is chemically interpretable:

1. allyl cation,
2. allyl anion,
3. benzene Kekulé/Dewar structures,
4. carboxylate,
5. nitrate,
6. ozone.

The implemented API is:

```python
options = {
    "include_nra": True,
    "nra_subspace": "pi",
    "nra_fit_metric": "frobenius",
}
```

The first results layout is:

```python
results["nra"] = {
    "subspace": "pi",
    "fit_metric": "frobenius",
    "weights": [...],
    "residual_norm": ...,
    "relative_residual": ...,
    "structures": [...],
}
```

The next report should present these data as a compact table with NRA weight, score/ranking weight, residual contribution, pi-bond pattern, and an automatic label where available. After that, the main scientific work is improving the structure pool and scoring model: formal-charge/octet penalties, sigma/pi separation, and radical/open-shell support.
