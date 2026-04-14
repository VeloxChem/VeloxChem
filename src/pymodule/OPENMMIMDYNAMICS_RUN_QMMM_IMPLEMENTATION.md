# OpenMMIMDynamics `run_qmmm` Migration Blueprint

## Scope

This document describes the required changes to make `OpenMMIMDynamics` consume interpolation databases produced by current `IMForceFieldGenerator` / `IMDatabasePointCollecter` logic.

Primary target:
- `openmmimdynamics.py` (`OpenMMIMDynamics`)

Reference implementation to mirror:
- `imdatabasepointcollecter.py` (`IMDatabasePointCollecter.run_qmmm` and setup helpers)

Important constraints from your request:
- This is a design + implementation blueprint only.
- No runtime code is changed in this step.

---

## Why the current setup is outdated

`OpenMMIMDynamics.run_immm` still assumes an older interpolation runtime contract:

1. Datapoint loading by label naming (`"_symmetry"`) instead of explicit `bank_role` metadata.
2. No runtime cache lifecycle (`mark_runtime_data_cache_dirty()`, `prepare_runtime_data_cache()`).
3. No MPI preload engine setup (`set_mpi_preload_engine`, dedicated interpolation communicator).
4. No separation between control-plane MPI and interpolation collective MPI.
5. No unified `update_settings(dynamics_settings, interpolation_settings, ...)` contract used by `IMForceFieldGenerator`.
6. No optimized observables path + verification path parity.
7. Missing test-hook/event contract used in integration harness.

The setup block in `IMDatabasePointCollecter.run_qmmm` now defines the canonical runtime initialization for database-backed interpolation dynamics.

---

## Current contract from `IMForceFieldGenerator`

`IMForceFieldGenerator` now drives runtime via:

```python
im_database_driver.update_settings(
    dynamics_settings,
    states_interpolation_settings,
    sampling_states_interpolation_settings,
)

stats = im_database_driver.run_qmmm(
    test_hooks={...},
    collect_step_trace=False,
    strict_test_hooks=self._test_hook_strict,
)
```

`dynamics_settings` includes at least:
- `drivers`
- `basis_set_label`
- `temperature`, `pressure`, `friction`, `timestep`, `nsteps`, `snapshots`, `ensemble`
- `trajectory_file`, `load_system`
- `desired_datapoint_density`, `converged_cycle`, `start_collect`, `roots_to_follow`
- `energy_threshold`, `grad_rmsd_thrsh`
- `sampling_drivers`, `sampling_settings`
- `metadynamics`
- `profile_runtime_timing`, `profile_interpolation_timing`
- `mpi_control_plane_enabled`, `mpi_root_worker_mode`, `mpi_reload_from_hdf5`, `mpi_debug_sync`

If `OpenMMIMDynamics` should be usable as a replacement runtime endpoint, it must accept and honor this contract.

---

## Recommended migration strategy

Use `IMDatabasePointCollecter` as backend engine for `run_qmmm` inside `OpenMMIMDynamics`.

This is the lowest-risk path because:
- setup semantics stay exactly synchronized with construction code,
- MPI/control-plane logic is reused instead of reimplemented,
- future runtime changes are inherited automatically.

If you want full in-class porting instead, see the "Full native port" section below.

---

## Drop-in section 1: imports for backend delegation

Add these imports in `openmmimdynamics.py`:

```python
from copy import deepcopy
from typing import Any, Dict, Optional

from .imdatabasepointcollecter import IMDatabasePointCollecter
```

---

## Drop-in section 2: normalize aliases + backend bridge helpers

Add the following methods to `OpenMMIMDynamics`.

```python
# Shared state references (do not deepcopy OpenMM objects)
_BACKEND_REF_FIELDS = (
    "system",
    "pdb",
    "topology",
    "molecule",
    "positions",
    "labels",
    "qm_atoms",
    "mm_subregion",
    "linking_atoms",
    "qm_force_index",
    "broken_bonds",
    "phase",
)

# Scalar/config/runtime containers
_BACKEND_CFG_FIELDS = (
    "platform",
    "openmm_precision",
    "ensemble",
    "temperature",
    "friction",
    "timestep",
    "nsteps",
    "snapshots",
    "pressure",
    "load_system",
    "out_file",
    "driver_flag",
    "distance_thrsh",
    "density_around_data_point",
    "allowed_molecule_deviation",
    "non_core_symmetry_groups",
    "all_rot_bonds",
    "identfy_relevant_int_coordinates",
    "use_symmetry",
    "cluster_run",
    "symmetry_rotors",
    "use_opt_confidence_radius",
    "bias_force_reaction_idx",
    "bias_force_reaction_prop",
)

# Values copied back to OpenMMIMDynamics after run
_BACKEND_RESULT_FIELDS = (
    "qm_potentials",
    "qm_mm_interaction_energies",
    "mm_potentials",
    "total_potentials",
    "kinetic_energies",
    "temperatures",
    "total_energies",
    "coordinates_xyz",
    "coordinates",
    "velocities",
    "velocities_np",
    "gradients",
    "state_specific_molecules",
    "point_adding_molecule",
    "allowed_molecules",
    "qm_data_point_dict",
    "qm_symmetry_datapoint_dict",
    "sorted_state_spec_im_labels",
    "impes_drivers",
    "current_state",
    "current_energy",
    "current_gradient",
)


def _normalize_runtime_aliases(self):
    """
    Normalizes old OpenMMIMDynamics attribute names to IMDatabasePointCollecter names.
    Keeps backwards compatibility with existing callers.
    """

    # Driver dictionary alias
    if not hasattr(self, "impes_drivers"):
        self.impes_drivers = getattr(self, "im_drivers", {})
    self.im_drivers = self.impes_drivers

    # Root-state molecule history alias
    if not hasattr(self, "state_specific_molecules"):
        self.state_specific_molecules = getattr(self, "root_spec_molecules", None)
    self.root_spec_molecules = self.state_specific_molecules

    # Interpolation settings alias
    if not hasattr(self, "interpolation_settings"):
        self.interpolation_settings = getattr(self, "interpolation_settings_dict", None)
    self.interpolation_settings_dict = self.interpolation_settings

    # Z-matrix alias: IMDatabasePointCollecter expects root-indexed mapping
    if not hasattr(self, "root_z_matrix"):
        if isinstance(getattr(self, "z_matrix", None), dict):
            self.root_z_matrix = self.z_matrix
        else:
            roots = list(getattr(self, "roots_to_follow", [0]))
            base_z = getattr(self, "z_matrix", None)
            self.root_z_matrix = {int(r): base_z for r in roots}


def _build_imdb_backend(self):
    """
    Build and seed an IMDatabasePointCollecter instance from current OpenMMIMDynamics state.
    """

    self._normalize_runtime_aliases()
    backend = IMDatabasePointCollecter(comm=self.comm, ostream=self.ostream)

    for name in self._BACKEND_REF_FIELDS:
        if hasattr(self, name):
            setattr(backend, name, getattr(self, name))

    for name in self._BACKEND_CFG_FIELDS:
        if hasattr(self, name):
            setattr(backend, name, deepcopy(getattr(self, name)))

    # Canonical runtime containers
    backend.roots_to_follow = list(getattr(self, "roots_to_follow", [0]))
    backend.root_z_matrix = deepcopy(getattr(self, "root_z_matrix", {}))
    backend.interpolation_settings = deepcopy(getattr(self, "interpolation_settings", None))
    backend.sampling_interpolation_settings = deepcopy(
        getattr(self, "sampling_interpolation_settings", None)
    )

    # Keep current naming in sync
    backend.impes_drivers = getattr(self, "impes_drivers", {})
    backend.state_specific_molecules = getattr(self, "state_specific_molecules", None)

    # If caller stored raw settings dicts, parse them through canonical parser.
    dyn = getattr(self, "dynamics_settings_interpolation_run", None)
    if dyn is not None:
        backend.update_settings(
            dyn,
            backend.interpolation_settings,
            backend.sampling_interpolation_settings,
        )

    return backend


def _sync_from_imdb_backend(self, backend):
    """Copy runtime results back to OpenMMIMDynamics."""

    for name in self._BACKEND_RESULT_FIELDS:
        if hasattr(backend, name):
            setattr(self, name, getattr(backend, name))

    # Keep aliases coherent
    self.im_drivers = getattr(self, "impes_drivers", {})
    self.root_spec_molecules = getattr(self, "state_specific_molecules", None)
    self.interpolation_settings_dict = getattr(self, "interpolation_settings", None)
```

---

## Drop-in section 3: canonical `update_settings` entrypoint

Add this method to `OpenMMIMDynamics` so callers can use the same contract as `IMForceFieldGenerator`:

```python
def update_settings(self, dynamics_settings, interpolation_settings=None, sampling_interpolation_settings=None):
    """
    Stores and normalizes settings using IMDatabasePointCollecter parser semantics.
    """

    parser = IMDatabasePointCollecter(comm=self.comm, ostream=self.ostream)
    parser.update_settings(
        dynamics_settings,
        interpolation_settings,
        sampling_interpolation_settings,
    )

    # Copy parsed settings onto this class.
    parsed_fields = (
        "drivers",
        "interpolation_settings",
        "sampling_interpolation_settings",
        "dynamics_settings_interpolation_run",
        "temperature",
        "starting_temperature",
        "pressure",
        "force_constant",
        "friction",
        "timestep",
        "nsteps",
        "snapshots",
        "load_system",
        "out_file",
        "reference_struc_energies_file",
        "platform",
        "openmm_precision",
        "ensemble",
        "metadynamics_settings",
        "metadynamics_enabled",
        "desired_datpoint_density",
        "unadded_cycles",
        "basis_set_label",
        "cluster_run",
        "start_collect",
        "qmc_stop",
        "nstates",
        "roots_to_follow",
        "excitation_pulse",
        "energy_threshold",
        "gradient_rmsd_thrsh",
        "sampling_enabled",
        "sampling_driver",
        "sampling_settings",
        "profile_runtime_timing",
        "profile_runtime_print_interval",
        "profile_interpolation_timing",
        "profile_interpolation_print_summary",
        "mpi_control_plane_enabled",
        "mpi_root_worker_mode",
        "mpi_reload_from_hdf5",
        "mpi_debug_sync",
    )

    for name in parsed_fields:
        if hasattr(parser, name):
            setattr(self, name, deepcopy(getattr(parser, name)))

    # Keep legacy attributes synchronized
    self.interpolation_settings_dict = self.interpolation_settings
```

---

## Drop-in section 4: add `run_qmmm` wrapper with full backend parity

Add this method to `OpenMMIMDynamics`:

```python
def run_qmmm(
    self,
    verify_optimized=False,
    verification_atol=1.0e-8,
    verification_rtol=1.0e-6,
    verification_fail_fast=False,
    use_optimized_observables=True,
    test_hooks=None,
    collect_step_trace=False,
    strict_test_hooks=False,
):
    """
    Run interpolation-based QM/MM dynamics using IMDatabasePointCollecter backend.
    """

    backend = self._build_imdb_backend()

    stats = backend.run_qmmm(
        verify_optimized=verify_optimized,
        verification_atol=verification_atol,
        verification_rtol=verification_rtol,
        verification_fail_fast=verification_fail_fast,
        use_optimized_observables=use_optimized_observables,
        test_hooks=test_hooks,
        collect_step_trace=collect_step_trace,
        strict_test_hooks=strict_test_hooks,
    )

    self._sync_from_imdb_backend(backend)
    return stats
```

---

## Drop-in section 5: keep old `run_immm` API as compatibility wrapper

Replace internals of `run_immm` with:

```python
def run_immm(
    self,
    im_driver,
    interpolation_settings,
    roots_to_follow=[0],
    restart_file=None,
    ensemble="NVE",
    temperature=298.15,
    pressure=1.0,
    friction=1.0,
    timestep=0.5,
    nsteps=1000,
    snapshots=100,
    traj_file="trajectory.pdb",
    state_file="output.xml",
    output_file="output",
):
    """
    Backwards-compatible wrapper around new settings + run_qmmm interface.
    """

    # Normalize interpolation settings to root-indexed dict
    if isinstance(interpolation_settings, dict) and all(
        isinstance(k, int) for k in interpolation_settings.keys()
    ):
        interpolation_settings_dict = interpolation_settings
    else:
        interpolation_settings_dict = {
            int(root): interpolation_settings[i]
            for i, root in enumerate(list(roots_to_follow))
        }

    # Build canonical dynamics settings expected by construction/runtime stack.
    dynamics_settings = {
        "drivers": getattr(self, "drivers", None),
        "temperature": float(temperature),
        "pressure": float(pressure),
        "friction": float(friction),
        "timestep": float(timestep),
        "nsteps": int(nsteps),
        "snapshots": int(snapshots),
        "trajectory_file": str(traj_file),
        "load_system": restart_file,
        "ensemble": str(ensemble),
        "roots_to_follow": list(roots_to_follow),
    }

    # Preserve existing per-class knobs if already set by caller.
    if hasattr(self, "dynamics_settings_interpolation_run") and self.dynamics_settings_interpolation_run:
        merged = dict(self.dynamics_settings_interpolation_run)
        merged.update(dynamics_settings)
        dynamics_settings = merged

    self.interpolation_settings = interpolation_settings_dict
    self.interpolation_settings_dict = interpolation_settings_dict
    self.dynamics_settings_interpolation_run = dynamics_settings

    self.update_settings(
        dynamics_settings,
        interpolation_settings_dict,
        getattr(self, "sampling_interpolation_settings", None),
    )

    return self.run_qmmm()
```

---

## Critical setup logic you must preserve (even if not using backend delegation)

If you choose a full native port, copy this setup pattern into `OpenMMIMDynamics.run_qmmm` before the MD loop.

### Mandatory initialization block

```python
self.allowed_molecules = {
    root: {
        "molecules": [],
        "qm_energies": [],
        "im_energies": [],
        "qm_gradients": [],
        "distances": [],
    }
    for root in self.roots_to_follow
}

self.state_specific_molecules = {root: [] for root in self.roots_to_follow}
self.sorted_state_spec_im_labels = {root: [] for root in self.roots_to_follow}
self.qm_data_point_dict = {root: [] for root in self.roots_to_follow}
self.qm_symmetry_datapoint_dict = {root: {} for root in self.roots_to_follow}
self.qm_energies_dict = {root: [] for root in self.roots_to_follow}
self.impes_drivers = {root: None for root in self.roots_to_follow}
self.sampled_molecules = {
    root: {
        "molecules": [],
        "im_energies": [],
        "qm_energies": [],
        "qm_gradients": [],
        "distances": [],
    }
    for root in self.roots_to_follow
}
self.qm_rotor_cluster_banks = {root: {} for root in self.roots_to_follow}
```

### Mandatory per-root driver bootstrap (HDF5 + runtime cache + MPI preload)

```python
masses = self.molecule.get_masses().copy()
masses_cart = np.repeat(masses, 3)
inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
self.inv_sqrt_masses = inv_sqrt_masses

for root in self.roots_to_follow:
    driver_object = InterpolationDriver(self.root_z_matrix[root])
    driver_object.update_settings(self.interpolation_settings[root])
    driver_object.enable_runtime_profiling(
        enabled=self.profile_interpolation_timing,
        reset=True,
        print_summary=False,
    )

    driver_object.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
    if root == 0:
        driver_object.symmetry_information = self.non_core_symmetry_groups["gs"]
    elif root == 1 and self.drivers["es"] is not None and self.drivers["es"][0].spin_flip:
        driver_object.symmetry_information = self.non_core_symmetry_groups["gs"]
    else:
        driver_object.symmetry_information = self.non_core_symmetry_groups["es"]

    driver_object.use_symmetry = self.use_symmetry
    self.impes_drivers[root] = driver_object

    # Rebuild datapoint banks from HDF5 using new bank_role semantics.
    self._reload_interpolation_root_from_hdf5(root, inv_sqrt_masses)

    self.prev_dens_of_points[root] = len(self.qm_data_point_dict[root])

    # Rotor-cluster registry payloads are now runtime-relevant.
    self.qm_rotor_cluster_banks[root] = self._load_rotor_cluster_bank_for_root(root)
    driver_object.qm_rotor_cluster_banks = self.qm_rotor_cluster_banks[root]

    if driver_object.qm_rotor_cluster_banks:
        first_family = next(iter(driver_object.qm_rotor_cluster_banks.values()))
        driver_object.rotor_cluster_information = first_family.get("cluster_info")
    else:
        driver_object.rotor_cluster_information = None
```

### Mandatory reload helper (use this exact style)

```python
def _reload_interpolation_root_from_hdf5(self, root, inv_sqrt_masses):
    driver_object = self.impes_drivers[root]

    im_labels, _ = driver_object.read_labels()

    self.qm_data_point_dict[root] = []
    self.qm_symmetry_datapoint_dict[root] = {}
    self.qm_energies_dict[root] = []
    self.sorted_state_spec_im_labels[root] = []

    old_label = None

    for label in im_labels:
        qm_data_point = InterpolationDatapoint(self.root_z_matrix[root])
        qm_data_point.update_settings(self.interpolation_settings[root])
        qm_data_point.read_hdf5(self.interpolation_settings[root]["imforcefield_file"], label)
        qm_data_point.inv_sqrt_masses = inv_sqrt_masses

        is_primary = (
            qm_data_point.bank_role == "core"
            or ("cluster" not in qm_data_point.point_label and "symmetry" not in qm_data_point.point_label)
        )

        if is_primary:
            self.qm_data_point_dict[root].append(qm_data_point)
            self.qm_energies_dict[root].append(qm_data_point.energy)
            self.sorted_state_spec_im_labels[root].append(label)

            old_label = qm_data_point.point_label
            self.qm_symmetry_datapoint_dict[root][old_label] = [qm_data_point]

        elif qm_data_point.bank_role == "symmetry":
            self.qm_symmetry_datapoint_dict[root][old_label].append(qm_data_point)

    driver_object.qm_symmetry_data_points = self.qm_symmetry_datapoint_dict[root]
    driver_object.qm_data_points = self.qm_data_point_dict[root]
    driver_object.labels = self.sorted_state_spec_im_labels[root]

    if len(self.qm_data_point_dict[root]) > 0:
        driver_object.impes_coordinate.eq_bond_lengths = self.qm_data_point_dict[root][0].eq_bond_lengths

    driver_object.mark_runtime_data_cache_dirty()
    driver_object.prepare_runtime_data_cache(force=True)

    use_mpi_preload = bool(self.interpolation_settings[root].get("use_mpi_preload", False))
    force_rebuild_preload = not (self._mpi_is_active() and self.mpi_root_worker_mode)

    driver_object.set_mpi_preload_engine(
        comm=self._mpi_interp_comm,
        enabled=(use_mpi_preload and self.nodes > 1),
        force_rebuild=force_rebuild_preload,
    )

    driver_object._mpi_preload_enabled_config = bool(
        getattr(driver_object, "mpi_preload_enabled", False)
    )
    self._ensure_interp_engine_comm(driver_object)

    if self._mpi_is_active() and self.mpi_root_worker_mode:
        driver_object.mpi_preload_enabled = False
```

---

## Drop-in section 6: optimized observable collection parity

Even if you do backend delegation now, keep these in a native-port branch for future decoupling.

```python
def _build_run_qmmm_runtime_cache(self):
    energy_unit = unit.kilojoules_per_mole
    gas_constant = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(energy_unit / unit.kelvin)
    return {
        "energy_unit": energy_unit,
        "qm_mm_groups": {1, 2},
        "mm_groups": {3, 4, 5, 6, 7},
        "dof": self.system.getNumParticles() * 3 - self.system.getNumConstraints(),
        "has_temp_method": hasattr(self.integrator, "computeSystemTemperature"),
        "gas_constant": gas_constant,
        "metadynamics_enabled": bool(getattr(self, "_metadynamics", None) is not None),
        "metadynamics_group": ({int(self._metadynamics_forcegroup)} if getattr(self, "_metadynamics", None) is not None else None),
    }


def _collect_step_observables_optimized(self, context, runtime_cache):
    energy_unit = runtime_cache["energy_unit"]

    qm = float(self.current_energy)
    qm_mm = context.getState(getEnergy=True, groups=runtime_cache["qm_mm_groups"]).getPotentialEnergy().value_in_unit(energy_unit)
    mm = context.getState(getEnergy=True, groups=runtime_cache["mm_groups"]).getPotentialEnergy().value_in_unit(energy_unit)
    kinetic = context.getState(getEnergy=True).getKineticEnergy().value_in_unit(energy_unit)

    meta_bias = 0.0
    if runtime_cache["metadynamics_enabled"]:
        group = runtime_cache["metadynamics_group"]
        meta_bias = context.getState(getEnergy=True, groups=group).getPotentialEnergy().value_in_unit(energy_unit)

    potential = qm + qm_mm + mm + meta_bias

    if runtime_cache["has_temp_method"]:
        temperature = self.integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
    else:
        temperature = (2.0 * kinetic) / (runtime_cache["dof"] * runtime_cache["gas_constant"])

    total = potential + kinetic

    return {
        "qm": qm,
        "qm_mm": qm_mm,
        "mm": mm,
        "metad_bias": float(meta_bias),
        "potential": potential,
        "kinetic": kinetic,
        "temperature": float(temperature),
        "total": total,
    }


def _append_step_observables(self, observables):
    self.qm_potentials.append(observables["qm"])
    self.qm_mm_interaction_energies.append(observables["qm_mm"])
    self.mm_potentials.append(observables["mm"])
    self.total_potentials.append(observables["potential"])
    self.kinetic_energies.append(observables["kinetic"])
    self.temperatures.append(observables["temperature"])
    self.total_energies.append(observables["total"])
```

---

## Drop-in section 7: force update path parity (native port mode)

```python
def _extract_qm_positions_nm(self, context):
    state = context.getState(getPositions=True)
    positions_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(positions_nm[self.qm_atoms], dtype=np.float64)


def _get_latest_qm_molecule(self, qm_positions_nm=None):
    molecule = getattr(self, "_latest_qm_molecule", None)
    if molecule is not None:
        if qm_positions_nm is None:
            return molecule
        cached_positions = getattr(self, "_latest_qm_positions_nm", None)
        current_positions = np.asarray(qm_positions_nm, dtype=np.float64)
        if cached_positions is not None and cached_positions.shape == current_positions.shape:
            if np.allclose(cached_positions, current_positions, atol=1.0e-12, rtol=0.0):
                return molecule

    if qm_positions_nm is None:
        return None

    return self._build_qm_molecule_from_positions(np.asarray(qm_positions_nm, dtype=np.float64))


def update_gradient_and_energy(self, new_positions, collective_mode=None):
    new_molecule = self._build_qm_molecule_from_positions(new_positions)
    self._latest_qm_molecule = new_molecule
    self._latest_qm_positions_nm = np.asarray(new_positions, dtype=np.float64).copy()

    self._compute_all_roots_for_molecule(new_molecule, collective_mode=bool(collective_mode))
    self._postprocess_root_state_selection()

    potential_kjmol = self.impes_drivers[self.current_state].impes_coordinate.energy * hartree_in_kjpermol()
    self.current_gradient = self.impes_drivers[self.current_state].impes_coordinate.gradient
    self.current_energy = potential_kjmol

    return self.current_gradient, potential_kjmol


def update_forces(self, context):
    if not hasattr(self, "_gradient_to_force_factor"):
        self._gradient_to_force_factor = (
            (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom())
            * unit.kilojoule_per_mole / unit.nanometer
        )

    qm_force = getattr(self, "_qm_external_force", None)
    if qm_force is None:
        qm_force = self.system.getForce(self.qm_force_index)
        self._qm_external_force = qm_force

    qm_positions = self._extract_qm_positions_nm(context)
    gradient = self.update_gradient(qm_positions)
    force = -np.array(gradient) * self._gradient_to_force_factor

    for i, atom_idx in enumerate(self.qm_atoms):
        qm_force.setParticleParameters(i, atom_idx, force[i])

    qm_force.updateParametersInContext(context)
```

---

## Full native port (if you do not delegate)

If backend delegation is not acceptable, you should copy these methods from `IMDatabasePointCollecter` into `OpenMMIMDynamics` (with name harmonization only):

1. `_set_test_hooks`, `_emit_test_hook`, `get_step_trace`
2. `update_settings`
3. `_build_run_qmmm_runtime_cache`
4. `_initialize_runtime_profilers`, `_add_runtime_timing`, `_finalize_runtime_step`
5. `_collect_step_observables_reference`, `_collect_step_observables_optimized`, `_append_step_observables`
6. `_reload_interpolation_root_from_hdf5`
7. `_initialize_sampling_impes_drivers` (if sampling is enabled)
8. `_ensure_interp_engine_comm`
9. `_build_qm_molecule_from_positions`, `_extract_qm_positions_nm`, `_get_latest_qm_molecule`
10. `_compute_all_roots_for_molecule`, `_postprocess_root_state_selection`
11. `update_gradient_and_energy`, `update_forces`, `get_qm_potential_energy`
12. Entire `run_qmmm` entrypoint including verification and print summary helpers.

For MPI parity (multi-rank runs), also port:
- `_mpi_is_active`, `_mpi_is_root`, `_mpi_bcast_control`, `_run_qmmm_worker_service`, `_mpi_worker_compute_step`, and associated command helpers.

---

## Integration checklist (must pass before enabling)

1. Single-rank, ground-state interpolation run reproduces old `run_immm` energies for same DB and initial geometry.
2. Runtime can read DBs containing `core`/`symmetry`/`cluster` bank roles without label-based assumptions.
3. `IMForceFieldGenerator`-style settings dict is accepted by `OpenMMIMDynamics.update_settings`.
4. `run_qmmm(test_hooks=...)` emits `run_start`, `step`, `run_end` payloads without crashing.
5. Optional MPI preload (`use_mpi_preload=True`) does not deadlock on multi-rank.
6. `OpenMMIMDynamics.run_immm(...)` remains callable and internally routes to `run_qmmm`.

---

## Recommended rollout order

1. Add delegation bridge (`update_settings`, `_build_imdb_backend`, `run_qmmm`, `run_immm` wrapper).
2. Validate behavior with existing interpolation integration tests.
3. Add native helper methods only if performance or architecture requires removing delegation.
4. After native parity is proven, optionally deprecate delegation path.

---

## Notes on file targets

When implementing, expected edits are primarily in:
- `openmmimdynamics.py`

No changes are required in:
- `imdatabasepointcollecter.py`
- `imforcefieldgenerator.py`

unless you decide to switch caller wiring from one runtime class to the other.

