# XTB Sampling-Driver Integration Plan

This document details how to implement a two-layer screening workflow:

1. Main interpolation remains based on the reference ab initio database.
2. A fast first-pass check uses an XTB driver + XTB-based interpolation database.
3. The XTB interpolation database is bootstrapped from the existing ab initio database and kept synchronized when new ab initio points are added.

The goal is to reduce expensive full ab initio checks while preserving final dataset quality.

---

## 1. Target Architecture

### Current behavior
- `IMDatabasePointCollecter.point_correlation_check()` performs direct ab initio energy/gradient checks against the main interpolation model and decides whether to add a point.
- `IMForceFieldGenerator.compute()` (notably `if self.add_structures_along_rcs`) also runs direct expensive checks.

### Desired behavior
- Add an XTB-first gate:
  - Compute XTB reference energy/gradient for current structure.
  - Compare against XTB interpolation prediction (built from XTB database).
  - If agreement is good, skip expensive full ab initio check.
  - If disagreement is large, run existing full ab initio check and add point if needed.
- Maintain two interpolation DB layers:
  - `abinitio` DB (existing, authoritative).
  - `xtb_sampling` DB (new, screening-only), initialized from ab initio geometries and incrementally updated when ab initio DB changes.

---

## 2. Critical Fix Before Feature Work

In `IMForceFieldGenerator.__init__`, the sampling gradient/hessian drivers are currently initialized with `ground_state_driver` instead of the new XTB sampling driver.

Current block (problematic):

```python
qm_driver = XtbDriver()
qm_grad_driver = XtbGradientDriver(ground_state_driver)
qm_hess_driver = XtbHessianDriver(ground_state_driver)
self.sampling_driver['gs'] = (qm_driver, qm_grad_driver, qm_hess_driver)
```

Required fix:

```python
if isinstance(ground_state_driver, XtbDriver):
    qm_grad_driver = XtbGradientDriver(ground_state_driver)
    qm_hess_driver = XtbHessianDriver(ground_state_driver)
    self.drivers['gs'] = (ground_state_driver, qm_grad_driver, qm_hess_driver)
    self.sampling_driver['gs'] = self.drivers['gs']  # same backend, no need for dual mode
else:
    sampling_qm = XtbDriver()
    sampling_grad = XtbGradientDriver(sampling_qm)
    sampling_hess = XtbHessianDriver(sampling_qm)
    self.sampling_driver['gs'] = (sampling_qm, sampling_grad, sampling_hess)
```

---

## 3. New Data Model and Settings

### 3.1 New fields in `IMForceFieldGenerator`

```python
self.sampling_enabled = True
self.sampling_screen_settings = {
    "energy_kcal_per_atom": 0.8,
    "gradient_rmsd_kcal_per_mol_ang": 2.0,
    "force_orient_cos": 0.8,
}
self.sampling_imforcefieldfiles = None
self.sampling_states_interpolation_settings = {root: None for root in roots_to_follow}
```

### 3.2 Extend `dynamics_settings` passed into `IMDatabasePointCollecter`

Add keys:

```python
'sampling_enabled': self.sampling_enabled,
'sampling_driver': self.sampling_driver,
'sampling_interpolation_settings': self.sampling_states_interpolation_settings,
'sampling_screen_settings': self.sampling_screen_settings,
```

---

## 4. New XTB Interpolation Database Files

Create dedicated XTB interpolation files per root:

```python
self.sampling_imforcefieldfiles = {
    root: f"im_database_xtb_{root}.h5" for root in self.roots_to_follow
}
```

Create matching interpolation settings (same interpolation knobs, different file):

```python
self.sampling_states_interpolation_settings[root] = {
    'interpolation_type': self.interpolation_type,
    'weightfunction_type': self.weightfunction_type,
    'exponent_p': self.exponent_p,
    'exponent_q': self.exponent_q,
    'confidence_radius': self.confidence_radius,
    'imforcefield_file': self.sampling_imforcefieldfiles[root],
    'use_inverse_bond_length': self.use_inverse_bond_length,
    'use_eq_bond_length': self.use_eq_bond_length,
    'use_cosine_dihedral': self.use_cosine_dihedral,
    'use_tc_weights': self.use_tc_weights,
    'use_mpi_preload': self.use_mpi_preload,
}
```

---

## 5. Bootstrap XTB DB From Existing Ab Initio DB

Implement a generator-side helper:

```python
def _bootstrap_sampling_db_from_abinitio_db(self, root):
    ref_settings = self.states_interpolation_settings[root]
    samp_settings = self.sampling_states_interpolation_settings[root]
    ref_file = ref_settings['imforcefield_file']
    samp_file = samp_settings['imforcefield_file']

    if not os.path.exists(ref_file):
        return

    ref_drv = InterpolationDriver(self.roots_z_matrix[root])
    ref_drv.update_settings(ref_settings)
    labels, z_matrix = ref_drv.read_labels()

    # Optional incremental update:
    existing = set()
    if os.path.exists(samp_file):
        samp_drv = InterpolationDriver(self.roots_z_matrix[root])
        samp_drv.update_settings(samp_settings)
        existing, _ = samp_drv.read_labels()
        existing = set(existing)

    mol_labels = self.molecule.get_labels()
    sampling_qm, sampling_grad, sampling_hess = self.sampling_driver['gs']

    for label in labels:
        if label in existing:
            continue

        ref_dp = InterpolationDatapoint(z_matrix)
        ref_dp.update_settings(ref_settings)
        ref_dp.read_hdf5(ref_file, label)

        coords_ang = ref_dp.cartesian_coordinates * bohr_in_angstrom()
        mol = Molecule(mol_labels, coords_ang, 'angstrom')
        mol.set_charge(self.molecule.get_charge())
        mol.set_multiplicity(self.molecule.get_multiplicity())

        # XTB reference compute
        e, _, _ = self.compute_energy(sampling_qm, mol, basis=None)
        g = self.compute_gradient(sampling_grad, mol, basis=None, scf_results=None, rsp_results=None)
        h = self.compute_hessian(sampling_hess, mol, basis=None)

        # Build datapoint in sampling DB using same label and geometry metadata
        masses = mol.get_masses().copy()
        inv_sqrt = 1.0 / np.sqrt(np.repeat(masses, 3))
        grad_vec = g[0].reshape(-1)
        hess_mat = h[0].reshape(grad_vec.size, grad_vec.size)
        mw_grad = inv_sqrt * grad_vec
        mw_hess = (inv_sqrt[:, None] * hess_mat) * inv_sqrt[None, :]

        samp_dp = InterpolationDatapoint(z_matrix)
        samp_dp.update_settings(samp_settings)
        samp_dp.cartesian_coordinates = ref_dp.cartesian_coordinates
        samp_dp.eq_bond_lengths = ref_dp.eq_bond_lengths
        samp_dp.mapping_masks = getattr(ref_dp, 'mapping_masks', None)
        samp_dp.imp_int_coordinates = getattr(ref_dp, 'imp_int_coordinates', [])
        samp_dp.inv_sqrt_masses = inv_sqrt
        samp_dp.energy = e[0]
        samp_dp.gradient = mw_grad.reshape(g[0].shape)
        samp_dp.hessian = mw_hess.reshape(h[0].shape)
        samp_dp.confidence_radius = getattr(ref_dp, 'confidence_radius', self.confidence_radius)
        samp_dp.transform_gradient_and_hessian()

        samp_dp.write_hdf5(samp_file, label)
```

Call this once after `self.states_interpolation_settings` is available and before production screening.

---

## 6. Runtime Integration in `IMDatabasePointCollecter`

### 6.1 Add new runtime members

In `IMDatabasePointCollecter.__init__`:

```python
self.sampling_enabled = False
self.sampling_driver = None
self.sampling_interpolation_settings = None
self.sampling_screen_settings = None
self.sampling_impes_drivers = None
self.sampling_qm_data_point_dict = None
self.sampling_qm_symmetry_datapoint_dict = None

self.sampling_screen_stats = {
    "attempted": 0,
    "passed_skip_full_qm": 0,
    "failed_run_full_qm": 0,
}
```

### 6.2 Parse in `update_settings`

```python
self.sampling_enabled = bool(dynamics_settings.get('sampling_enabled', False))
self.sampling_driver = dynamics_settings.get('sampling_driver', None)
self.sampling_interpolation_settings = dynamics_settings.get('sampling_interpolation_settings', None)
self.sampling_screen_settings = dynamics_settings.get(
    'sampling_screen_settings',
    {"energy_kcal_per_atom": 0.8, "gradient_rmsd_kcal_per_mol_ang": 2.0, "force_orient_cos": 0.8}
)
```

### 6.3 Initialize sampling interpolation drivers in `run_qmmm`

After main `self.impes_drivers` initialization, add:

```python
if self.sampling_enabled:
    self._initialize_sampling_impes_drivers(inv_sqrt_masses)
```

Helper:

```python
def _initialize_sampling_impes_drivers(self, inv_sqrt_masses):
    self.sampling_impes_drivers = {root: None for root in self.roots_to_follow}
    self.sampling_qm_data_point_dict = {root: [] for root in self.roots_to_follow}
    self.sampling_qm_symmetry_datapoint_dict = {root: {} for root in self.roots_to_follow}

    for root in self.roots_to_follow:
        drv = InterpolationDriver(self.root_z_matrix[root])
        drv.update_settings(self.sampling_interpolation_settings[root])
        drv.symmetry_information = self.impes_drivers[root].symmetry_information
        drv.use_symmetry = self.use_symmetry
        drv.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses

        labels, _ = drv.read_labels()
        old_label = None
        for label in labels:
            dp = InterpolationDatapoint(self.root_z_matrix[root])
            dp.update_settings(self.sampling_interpolation_settings[root])
            dp.read_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)
            dp.inv_sqrt_masses = inv_sqrt_masses
            if '_symmetry' not in label:
                self.sampling_qm_data_point_dict[root].append(dp)
                old_label = dp.point_label
                self.sampling_qm_symmetry_datapoint_dict[root][old_label] = [dp]
            else:
                self.sampling_qm_symmetry_datapoint_dict[root][old_label].append(dp)

        drv.qm_data_points = self.sampling_qm_data_point_dict[root]
        drv.qm_symmetry_data_points = self.sampling_qm_symmetry_datapoint_dict[root]
        drv.prepare_runtime_data_cache(force=True)
        self.sampling_impes_drivers[root] = drv
```

---

## 7. Two-Stage Screening Logic in `point_correlation_check`

### 7.1 Stage A: XTB screen

Add helper:

```python
def _sampling_screen_check(self, molecule):
    if not self.sampling_enabled or self.sampling_driver is None or self.sampling_impes_drivers is None:
        return {"skip_full_qm": False, "details": {}}

    self.sampling_screen_stats["attempted"] += 1
    natms = len(molecule.get_labels())

    sampling_qm, sampling_grad, _ = self.sampling_driver['gs']
    xtb_energy, scf_results, rsp_results = self._compute_energy_mpi_safe(sampling_qm, molecule, basis=None)
    xtb_grad = self._compute_gradient_mpi_safe(sampling_grad, molecule, basis=None, scf_results=None, rsp_results=None)

    e_thr = float(self.sampling_screen_settings["energy_kcal_per_atom"])
    g_thr = float(self.sampling_screen_settings["gradient_rmsd_kcal_per_mol_ang"])
    c_thr = float(self.sampling_screen_settings["force_orient_cos"])

    all_ok = True
    details = {}
    for root in self.roots_to_follow:
        self._interp_compute_root_local_serial(self.sampling_impes_drivers[root], molecule)
        e_im = self.sampling_impes_drivers[root].impes_coordinate.energy
        g_im = self.sampling_impes_drivers[root].impes_coordinate.gradient

        e_diff = abs(float(xtb_energy[0]) - float(e_im)) / natms * hartree_in_kcalpermol()
        g_diff = (xtb_grad[0] - g_im) * hartree_in_kcalpermol() * bohr_in_angstrom()
        g_rmsd = float(np.sqrt((g_diff ** 2).mean()))

        gq = xtb_grad[0].ravel()
        gi = g_im.ravel()
        denom = np.linalg.norm(gq) * np.linalg.norm(gi)
        cos_theta = 1.0 if denom < 1.0e-15 else float(np.dot(gq, gi) / denom)

        details[root] = {"e_diff_kcal_per_atom": e_diff, "g_rmsd": g_rmsd, "cos": cos_theta}
        if not (e_diff <= e_thr and g_rmsd <= g_thr and cos_theta >= c_thr):
            all_ok = False

    if all_ok:
        self.sampling_screen_stats["passed_skip_full_qm"] += 1
        return {"skip_full_qm": True, "details": details}

    self.sampling_screen_stats["failed_run_full_qm"] += 1
    return {"skip_full_qm": False, "details": details}
```

### 7.2 Integrate at start of `point_correlation_check()`

```python
screen = self._sampling_screen_check(molecule)
if screen["skip_full_qm"]:
    # fast accept: no point addition, keep allowed_molecules bookkeeping
    for root in self.roots_to_follow:
        self.allowed_molecules[root]['molecules'].append(molecule)
        self.allowed_molecules[root]['im_energies'].append(self.impes_drivers[root].impes_coordinate.energy)
        self.allowed_molecules[root]['qm_energies'].append(self.impes_drivers[root].impes_coordinate.energy)
        self.allowed_molecules[root]['qm_gradients'].append(self.impes_drivers[root].impes_coordinate.gradient)
    self.add_a_point = False
    return
```

Then keep existing full ab initio path unchanged as Stage B.

---

## 8. Keep XTB DB Synced When Ab Initio Adds New Points

When `IMDatabasePointCollecter.add_point()` writes a new ab initio datapoint:

Current write point:

```python
impes_coordinate.write_hdf5(self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'], new_label)
```

Immediately mirror into XTB DB:

```python
if self.sampling_enabled:
    self._write_sampling_point_from_geometry(
        root=mol_basis[4][number],
        molecule=mol_basis[0],
        label=new_label,
        template_point=impes_coordinate
    )
```

Helper:

```python
def _write_sampling_point_from_geometry(self, root, molecule, label, template_point):
    sampling_qm, sampling_grad, sampling_hess = self.sampling_driver['gs']
    e, _, _ = self._compute_energy_mpi_safe(sampling_qm, molecule, basis=None)
    g = self._compute_gradient_mpi_safe(sampling_grad, molecule, basis=None, scf_results=None, rsp_results=None)
    h = self._compute_hessian_mpi_safe(sampling_hess, molecule, basis=None)

    masses = molecule.get_masses().copy()
    inv_sqrt = 1.0 / np.sqrt(np.repeat(masses, 3))
    grad_vec = g[0].reshape(-1)
    hess_mat = h[0].reshape(grad_vec.size, grad_vec.size)
    mw_grad = inv_sqrt * grad_vec
    mw_hess = (inv_sqrt[:, None] * hess_mat) * inv_sqrt[None, :]

    dp = InterpolationDatapoint(self.root_z_matrix[root])
    dp.update_settings(self.sampling_interpolation_settings[root])
    dp.cartesian_coordinates = template_point.cartesian_coordinates
    dp.eq_bond_lengths = template_point.eq_bond_lengths
    dp.mapping_masks = getattr(template_point, "mapping_masks", None)
    dp.imp_int_coordinates = getattr(template_point, "imp_int_coordinates", [])
    dp.inv_sqrt_masses = inv_sqrt
    dp.energy = e[0]
    dp.gradient = mw_grad.reshape(g[0].shape)
    dp.hessian = mw_hess.reshape(h[0].shape)
    dp.confidence_radius = template_point.confidence_radius
    dp.transform_gradient_and_hessian()

    if self._mpi_is_root() or (not self._mpi_is_active()):
        dp.write_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)

    # Refresh in-memory sampling interpolator
    self._reload_sampling_root_from_hdf5(root, inv_sqrt)
```

---

## 9. Integrate Into `if self.add_structures_along_rcs:` Path

In `IMForceFieldGenerator.compute()` under `if self.add_structures_along_rcs`:

- Before loop, ensure sampling DB is bootstrapped.
- Replace direct expensive check with staged check:

```python
screen_ok = self._sampling_screen_for_structure(mol, state=states)
if screen_ok:
    continue  # no expensive ab initio check needed

# Existing expensive check:
# energy_abinitio vs interpolation_abinitio
# if threshold exceeded -> self.add_point(...)
```

Suggested helper in generator:

```python
def _sampling_screen_for_structure(self, mol, state=0):
    root = int(state)
    basis = None  # XTB path
    e_xtb, _, _ = self.compute_energy(self.sampling_driver['gs'][0], mol, basis)

    drv = InterpolationDriver(self.roots_z_matrix[root])
    drv.update_settings(self.sampling_states_interpolation_settings[root])
    drv.compute(mol)
    e_im = drv.impes_coordinate.energy

    e_diff = abs(float(e_xtb[0]) - float(e_im)) * hartree_in_kcalpermol()
    return e_diff <= self.sampling_screen_settings["energy_kcal_per_atom"] * len(mol.get_labels())
```

---

## 10. MPI and Synchronization Notes

- Reuse existing collective helper pattern (`_collective_call_inline_style`) for XTB screen compute when needed.
- Keep HDF5 write ownership on root rank only.
- After root write, reload/refresh sampling interpolation caches on all ranks (same as main interpolation cache sync approach).

---

## 11. Telemetry and Safety

Add counters and periodic prints:

```python
print(
    "Sampling screen stats:",
    self.sampling_screen_stats["attempted"],
    self.sampling_screen_stats["passed_skip_full_qm"],
    self.sampling_screen_stats["failed_run_full_qm"],
)
```

Add hard fallback behavior:

```python
if self.sampling_enabled and (self.sampling_driver is None or self.sampling_driver['gs'] is None):
    print("Sampling enabled but sampling_driver unavailable. Falling back to full ab initio checks.")
    self.sampling_enabled = False
```

---

## 12. Validation Plan

1. Bootstrapping test:
- Create `im_database_xtb_0.h5` from existing `im_database_0.h5`.
- Assert label counts match.

2. Screening consistency test:
- For random points from trajectory, compare:
  - XTB vs XTB-interpolation error
  - ab initio vs ab initio-interpolation error
- Confirm XTB screen rejects most high-error ab initio points.

3. Incremental sync test:
- Add one ab initio point during simulation.
- Assert same label appears in XTB DB.

4. Regression test:
- Disable sampling and verify existing point-addition behavior is unchanged.

5. MPI test:
- 2-rank run in root-worker mode with point additions.
- Ensure no deadlocks and both DBs updated on disk.

---

## 13. Suggested Implementation Order

1. Fix `sampling_driver` initialization wiring.
2. Add sampling DB filenames/settings and bootstrap helper.
3. Pass sampling config into `IMDatabasePointCollecter`.
4. Add sampling driver/interpolator initialization in collector.
5. Add Stage-A XTB screen in `point_correlation_check`.
6. Add DB synchronization hook in `IMDatabasePointCollecter.add_point`.
7. Integrate optional screen in `add_structures_along_rcs`.
8. Add logging + tests.

---

## 14. Known Risks and Mitigations

- Risk: XTB and ab initio PES mismatch causes false negatives in screening.
  - Mitigation: conservative thresholds and mandatory full-check fallback for borderline cases.
- Risk: dual DB drift (labels out of sync).
  - Mitigation: same label reuse and post-write integrity assertions.
- Risk: MPI cache divergence.
  - Mitigation: root-only writes + explicit reload on workers.

