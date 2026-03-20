# Grouped Trust-Radius Optimization (Distance-Based)

This implementation guide describes the exact code changes needed to group datapoints and reference structures using the existing Cartesian alignment distance, then optimize trust radii in subsets instead of globally.

The core idea is:
- Build `datapoint <-> reference` neighborhoods from `cartesian_distance` behavior.
- Derive connected datapoint groups via shared references.
- Skip/freeze orphan or strongly imbalanced groups.
- Run `AlphaOptimizer` per valid group.
- Write trust radii back exactly as today.

---

## 1) Files to change

1. `imdatabasepointcollecter.py` (main implementation)
2. `imforcefieldgenerator.py` (pass grouping config to collector; optional but recommended)

---

## 2) Why this fixes the current issue

Current global optimization can force unexplored PES datapoints to very small radii because they are evaluated against references that are not local to them.

Grouping prevents this by:
- only optimizing a datapoint with nearby references,
- freezing datapoints with insufficient local support,
- allowing reference reuse across cycles through persistent group cache.

---

## 3) `imdatabasepointcollecter.py` changes

### 3.1 Add grouping defaults in `__init__`

Add below the existing trust-radius settings:

```python
# Trust-radius grouping controls (distance-based)
self.trust_radius_grouping = {
    "enabled": True,
    "k_ref_neighbors": 8,          # always include k nearest references per datapoint
    "d_link_bohr": 0.9,            # include any reference within this cartesian distance (Bohr)
    "min_refs_per_group": 4,       # hard minimum references to optimize a group
    "min_ref_per_dp_ratio": 1.0,   # group balance lower bound
    "max_ref_per_dp_ratio": 25.0,  # group balance upper bound
    "freeze_orphans": True,        # keep old alpha for unsupported points
    "cache_max_age": 8,            # cycles to keep stale groups
}

# Cache persists grouping identity across optimization cycles
# Structure: {root: {"cycle": int, "groups": [group_dict, ...]}}
self._trust_radius_group_cache = {}
self._trust_radius_group_cycle = 0
```

### Explanation
- `k_ref_neighbors` prevents isolated points from seeing zero references.
- `d_link_bohr` captures local PES neighborhoods by geometric proximity.
- ratio bounds catch "very miss balanced" groups.
- cache fields make groups reusable across optimization cycles.

---

### 3.2 Add helper to pick symmetry info per root (reuses existing logic)

```python
def _get_symmetry_info_for_root(self, root):
    """Return symmetry info used by distance alignment for this electronic state root."""
    if root == 0:
        return self.non_core_symmetry_groups['gs']
    if root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
        return self.non_core_symmetry_groups['gs']
    return self.non_core_symmetry_groups['es']
```

### Explanation
This avoids duplicating state-selection logic and guarantees grouping uses the same active atom subset as interpolation distance handling.

---

### 3.3 Add pairwise distance matrix builder (based on current Cartesian alignment)

```python
def _build_dp_ref_distance_matrix(self, datapoints, molecules, symmetry_info):
    """
    D[i, j] = cartesian distance between datapoint i and reference molecule j
    using calculate_distance_to_ref (translation + rotational alignment + symmetry filtering).
    """
    n_dp = len(datapoints)
    n_ref = len(molecules)
    D = np.full((n_dp, n_ref), np.inf, dtype=np.float64)

    for i, dp in enumerate(datapoints):
        dp_xyz = dp.cartesian_coordinates
        for j, mol in enumerate(molecules):
            mol_xyz = mol.get_coordinates_in_bohr()
            _, dist, _ = self.calculate_distance_to_ref(mol_xyz, dp_xyz, symmetry_info)
            D[i, j] = float(dist)

    return D
```

### Explanation
This intentionally uses your existing `calculate_distance_to_ref(...)` instead of introducing another metric, so grouping and interpolation remain consistent.

---

### 3.4 Add group construction via shared local references

```python
def _build_trust_radius_groups(self, root, datapoints, molecules):
    """
    Build datapoint groups from shared reference neighborhoods.

    Returns dict with:
      - groups: list of {dp_idx, ref_idx, status, ratio}
      - frozen_dp_idx: set of datapoints to keep unchanged
      - distance_matrix: ndarray
    """
    cfg = self.trust_radius_grouping
    n_dp = len(datapoints)
    n_ref = len(molecules)

    out = {
        "groups": [],
        "frozen_dp_idx": set(),
        "distance_matrix": np.zeros((n_dp, n_ref), dtype=np.float64),
    }

    if n_dp == 0:
        return out

    symmetry_info = self._get_symmetry_info_for_root(root)
    D = self._build_dp_ref_distance_matrix(datapoints, molecules, symmetry_info)
    out["distance_matrix"] = D

    # 1) neighborhood refs per datapoint
    neighborhoods = []
    k_ref = int(cfg["k_ref_neighbors"])
    d_link = float(cfg["d_link_bohr"])

    for i in range(n_dp):
        if n_ref == 0:
            neighborhoods.append(set())
            continue

        nearest = np.argsort(D[i])[:min(k_ref, n_ref)]
        local = set(int(x) for x in nearest.tolist())

        within = np.where(D[i] <= d_link)[0]
        local.update(int(x) for x in within.tolist())

        neighborhoods.append(local)

    # 2) datapoint graph: connected if sharing at least one reference neighbor
    ref_to_dp = {}
    for i, refs in enumerate(neighborhoods):
        for r in refs:
            ref_to_dp.setdefault(r, []).append(i)

    adjacency = [set() for _ in range(n_dp)]
    for _, dplist in ref_to_dp.items():
        s = set(dplist)
        for i in s:
            adjacency[i].update(s)

    # 3) connected components on datapoints
    visited = set()
    components = []

    for i in range(n_dp):
        if i in visited:
            continue
        stack = [i]
        comp = []
        while stack:
            u = stack.pop()
            if u in visited:
                continue
            visited.add(u)
            comp.append(u)
            stack.extend(v for v in adjacency[u] if v not in visited)
        components.append(sorted(comp))

    # 4) per-component reference union and quality flags
    min_refs = int(cfg["min_refs_per_group"])
    min_ratio = float(cfg["min_ref_per_dp_ratio"])
    max_ratio = float(cfg["max_ref_per_dp_ratio"])

    for comp in components:
        ref_union = set()
        for i in comp:
            ref_union.update(neighborhoods[i])

        n_comp = len(comp)
        n_refs = len(ref_union)
        ratio = n_refs / max(1, n_comp)

        status = "optimize"
        if n_refs == 0:
            status = "orphan"
        elif n_refs < min_refs:
            status = "too_few_refs"
        elif ratio < min_ratio:
            status = "imbalanced_low"
        elif ratio > max_ratio:
            status = "imbalanced_high"

        group = {
            "dp_idx": comp,
            "ref_idx": sorted(ref_union),
            "status": status,
            "ratio": ratio,
            "n_dp": n_comp,
            "n_ref": n_refs,
        }
        out["groups"].append(group)

        if status != "optimize" and bool(cfg["freeze_orphans"]):
            out["frozen_dp_idx"].update(comp)

    return out
```

### Explanation
- Grouping is **purely distance-driven** and uses your requested Cartesian metric.
- Shared references create natural PES neighborhoods.
- Orphan and imbalanced groups are excluded from optimization to avoid artificial shrinkage.

---

### 3.5 Add grouped optimization orchestrator

```python
def determine_trust_radius_gradient_grouped(
    self,
    root,
    molecules,
    qm_energies,
    qm_gradients,
    im_energies,
    datapoints,
    interpolation_setting,
    sym_datapoints,
    sym_dict,
    z_matrix,
    exponent_p_q,
):
    """
    Group-aware trust-radius optimization wrapper.
    Returns full alpha list aligned with datapoints.
    """
    base_alphas = [float(dp.confidence_radius) for dp in datapoints]

    if not self.trust_radius_grouping.get("enabled", False):
        return self.determine_trust_radius_gradient(
            molecules[:], qm_energies[:], qm_gradients[:], im_energies[:],
            datapoints, interpolation_setting, sym_datapoints, sym_dict,
            z_matrix, exponent_p_q,
        )

    grouping = self._build_trust_radius_groups(root, datapoints, molecules)
    groups = grouping["groups"]

    # Persist lightweight cache for reference reuse/inspection across cycles
    self._trust_radius_group_cache[root] = {
        "cycle": self._trust_radius_group_cycle,
        "groups": groups,
    }

    final_alphas = base_alphas[:]

    for g in groups:
        if g["status"] != "optimize":
            continue

        dp_idx = g["dp_idx"]
        ref_idx = g["ref_idx"]
        if len(dp_idx) == 0 or len(ref_idx) == 0:
            continue

        sub_dps = [datapoints[i] for i in dp_idx]
        sub_mols = [molecules[j] for j in ref_idx]
        sub_qm_e = [qm_energies[j] for j in ref_idx]
        sub_qm_g = [qm_gradients[j] for j in ref_idx]
        sub_im_e = [im_energies[j] for j in ref_idx]

        # Reuse existing optimizer unchanged
        sub_alpha = self.determine_trust_radius_gradient(
            sub_mols[:], sub_qm_e[:], sub_qm_g[:], sub_im_e[:],
            sub_dps,
            interpolation_setting,
            sym_datapoints,
            sym_dict,
            z_matrix,
            exponent_p_q,
        )

        for local_i, global_i in enumerate(dp_idx):
            final_alphas[global_i] = float(sub_alpha[local_i])

    self._trust_radius_group_cycle += 1
    return final_alphas
```

### Explanation
- This function is intentionally a wrapper, so `AlphaOptimizer` remains unchanged.
- Frozen groups simply keep the previous alpha.
- Valid groups are optimized independently, which prevents unsupported points from dominating global fitting.

---

### 3.6 Integrate grouped call at the existing optimization trigger

Locate the current block in `update_forces(...)` where `multi_grad` is executed (around current line ~4163).

Replace:

```python
trust_radius = self.determine_trust_radius_gradient(...)
```

with:

```python
trust_radius = self.determine_trust_radius_gradient_grouped(
    root=self.current_state,
    molecules=chosen_structures,
    qm_energies=chosen_qm_energies,
    qm_gradients=chosen_qm_gradients,
    im_energies=chosen_im_energies,
    datapoints=self.qm_data_point_dict[self.current_state],
    interpolation_setting=self.interpolation_settings[self.current_state],
    sym_datapoints=self.qm_symmetry_datapoint_dict[self.current_state],
    sym_dict=sym_dict,
    z_matrix=self.root_z_matrix[self.current_state],
    exponent_p_q=(
        self.impes_drivers[self.current_state].exponent_p,
        self.impes_drivers[self.current_state].exponent_q,
    ),
)
```

### Explanation
This is the only runtime call-site change needed. Existing write-back (`update_confidence_radius`, in-memory updates, cache refresh) remains identical.

---

## 4) `imforcefieldgenerator.py` optional wiring

To expose grouping controls from the generator layer:

### 4.1 Add defaults in `IMForceFieldGenerator.__init__`

```python
self.trust_radius_grouping = {
    "enabled": True,
    "k_ref_neighbors": 8,
    "d_link_bohr": 0.9,
    "min_refs_per_group": 4,
    "min_ref_per_dp_ratio": 1.0,
    "max_ref_per_dp_ratio": 25.0,
    "freeze_orphans": True,
    "cache_max_age": 8,
}
```

### 4.2 Pass config to collector before `run_qmmm()`

At both places where `im_database_driver` is configured:

```python
im_database_driver.trust_radius_grouping = copy.deepcopy(self.trust_radius_grouping)
```

### Explanation
This keeps tuning in one place (generator), while collector stays execution-focused.

---

## 5) Logging highlights to add

Add one compact log right before grouped optimization loop:

```python
print(
    f"[TR-group] root={root} groups={len(groups)} "
    f"opt={sum(g['status']=='optimize' for g in groups)} "
    f"frozen={sum(g['status']!='optimize' for g in groups)}"
)
```

And per skipped group:

```python
if g["status"] != "optimize":
    print(f"[TR-group] skip status={g['status']} n_dp={g['n_dp']} n_ref={g['n_ref']} ratio={g['ratio']:.2f}")
```

This is very useful to tune `d_link_bohr` and ratio thresholds.

---

## 6) Notes and edge cases

1. `determine_trust_radius_gradient(...)` mutates its local `molecules`/`energies` list by appending datapoints as structures. That is okay here because grouped wrapper passes copies.
2. If a group has no references, it should not be optimized; freezing is the safest behavior.
3. Use a conservative initial `d_link_bohr` and inspect skip statistics before tightening.
4. Existing meaning of `self.use_opt_confidence_radius[2]` is overloaded in current code (used both as default radius and `e_x` weight). Grouping implementation does not change this behavior; consider splitting these into separate explicit keys later.

---

## 7) Minimal test checklist

1. Run with grouping disabled (`enabled=False`) and confirm identical trust-radii behavior.
2. Enable grouping and confirm logs show multiple groups.
3. Verify orphan groups are skipped/frozen (no sudden collapse of alpha).
4. Confirm HDF5 + in-memory confidence radii are still synchronized.
5. Confirm interpolation cache refresh still occurs once per optimization cycle.

---

## 8) Recommended first parameter set

```python
{
    "enabled": True,
    "k_ref_neighbors": 8,
    "d_link_bohr": 0.9,
    "min_refs_per_group": 4,
    "min_ref_per_dp_ratio": 1.0,
    "max_ref_per_dp_ratio": 20.0,
    "freeze_orphans": True,
    "cache_max_age": 8,
}
```

Then tune `d_link_bohr` and `k_ref_neighbors` from observed group skip rates.

