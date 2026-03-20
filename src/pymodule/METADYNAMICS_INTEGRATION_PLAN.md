# Metadynamics Integration Plan for `IMDatabasePointCollecter`

This document outlines a detailed implementation path to integrate OpenMM metadynamics into `IMDatabasePointCollecter.run_qmmm()` so QM/MM dynamics can escape local minima and explore broader configurational space.

The approach below is designed to:
- Keep existing behavior unchanged when metadynamics is disabled.
- Reuse the current QM/MM update flow (`update_forces -> integrator step`).
- Preserve MPI root-worker synchronization behavior.
- Keep energy bookkeeping explicit (QM, QM/MM, MM, metad bias, total).

## 1. Scope and Design Goals

### Primary goal
Add optional, configurable well-tempered metadynamics on one or more collective variables (CVs) within the existing QM/MM loop in `imdatabasepointcollecter.py`.

### Non-goals for first integration
- No redesign of interpolation/QM driver logic.
- No change to existing `add_bias_force()` API.
- No replacement of old bias-force path yet; keep coexistence but enforce clear guardrails.

## 2. Configuration Schema

Add a new optional block in `dynamics_settings`:

```python
dynamics_settings["metadynamics"] = {
    "enabled": True,
    "bias_factor": 10.0,
    "hill_height_kjmol": 1.2,
    "hill_frequency": 250,
    "save_frequency": 5000,      # optional
    "bias_dir": "meta_bias",     # optional
    "force_group": 11,           # optional, default 11
    "variables": [
        {
            "type": "torsion",   # distance | angle | torsion
            "atoms": [3, 5, 7, 9],
            "min_deg": -180.0,
            "max_deg": 180.0,
            "width_deg": 12.0,
            "periodic": True
        },
        {
            "type": "distance",
            "atoms": [1, 14],
            "min_nm": 0.15,
            "max_nm": 0.45,
            "width_nm": 0.01,
            "periodic": False
        }
    ]
}
```

## 3. Imports and Class Runtime State

### 3.1 Add metadynamics imports (safe optional import)
Insert near existing OpenMM imports:

```python
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm.app.metadynamics import Metadynamics, BiasVariable
except ImportError:
    Metadynamics = None
    BiasVariable = None
    pass
```

### 3.2 Add runtime members in `__init__`
Add near other simulation/bias state members:

```python
# Metadynamics configuration/runtime
self.metadynamics_settings = None
self.metadynamics_enabled = False
self._metadynamics = None
self._metadynamics_variables = []
self._metadynamics_cv_forces = []
self._metadynamics_force_group = 11

# Per-step diagnostics
self.metad_bias_energies = []
self.metad_cv_history = []   # list[list[float]], per step
```

## 4. Settings Parsing and Validation

### 4.1 Parse in `update_settings()`
Append in `update_settings()`:

```python
if 'metadynamics' in dynamics_settings:
    self.metadynamics_settings = copy.deepcopy(dynamics_settings['metadynamics'])
else:
    self.metadynamics_settings = None

self.metadynamics_enabled = bool(
    self.metadynamics_settings is not None
    and self.metadynamics_settings.get('enabled', False)
)

if self.metadynamics_enabled:
    self._validate_metadynamics_settings()
```

### 4.2 Validation helper
Add helper method:

```python
def _validate_metadynamics_settings(self):
    cfg = self.metadynamics_settings
    if Metadynamics is None or BiasVariable is None:
        raise RuntimeError("Metadynamics requested, but openmm.app.metadynamics is unavailable.")

    required = ('variables', 'bias_factor', 'hill_height_kjmol', 'hill_frequency')
    missing = [k for k in required if k not in cfg]
    if missing:
        raise ValueError(f"metadynamics missing required keys: {missing}")

    variables = cfg['variables']
    if not isinstance(variables, list) or len(variables) == 0:
        raise ValueError("metadynamics.variables must be a non-empty list.")

    natoms = self.system.getNumParticles() if self.system is not None else None
    for idx, var in enumerate(variables):
        vtype = var.get('type', '').lower()
        atoms = var.get('atoms', [])
        if vtype not in ('distance', 'angle', 'torsion'):
            raise ValueError(f"metadynamics.variables[{idx}].type must be distance|angle|torsion.")
        expected = {'distance': 2, 'angle': 3, 'torsion': 4}[vtype]
        if len(atoms) != expected:
            raise ValueError(f"metadynamics.variables[{idx}].atoms must have {expected} entries.")
        if natoms is not None:
            for a in atoms:
                if a < 0 or a >= natoms:
                    raise ValueError(f"metadynamics.variables[{idx}] atom index {a} out of range [0, {natoms-1}].")

        # Unit-range keys
        if vtype in ('angle', 'torsion'):
            for k in ('min_deg', 'max_deg', 'width_deg'):
                if k not in var:
                    raise ValueError(f"metadynamics.variables[{idx}] missing key '{k}'.")
        else:
            for k in ('min_nm', 'max_nm', 'width_nm'):
                if k not in var:
                    raise ValueError(f"metadynamics.variables[{idx}] missing key '{k}'.")

    if float(cfg['bias_factor']) <= 1.0:
        raise ValueError("metadynamics.bias_factor must be > 1.0 for well-tempered metadynamics.")
    if float(cfg['hill_height_kjmol']) <= 0.0:
        raise ValueError("metadynamics.hill_height_kjmol must be > 0.")
    if int(cfg['hill_frequency']) <= 0:
        raise ValueError("metadynamics.hill_frequency must be > 0.")
```

## 5. CV Force and Metadynamics Object Construction

### 5.1 Build one CV force from config
Add helper:

```python
def _build_metadynamics_cv_force(self, spec):
    vtype = spec['type'].lower()
    atoms = spec['atoms']

    if vtype == 'distance':
        cv_force = mm.CustomBondForce("r")
        cv_force.addBond(int(atoms[0]), int(atoms[1]), [])
        min_val = float(spec['min_nm'])
        max_val = float(spec['max_nm'])
        width = float(spec['width_nm'])
        periodic = bool(spec.get('periodic', False))

    elif vtype == 'angle':
        cv_force = mm.CustomAngleForce("theta")
        cv_force.addAngle(int(atoms[0]), int(atoms[1]), int(atoms[2]), [])
        min_val = np.deg2rad(float(spec['min_deg']))
        max_val = np.deg2rad(float(spec['max_deg']))
        width = np.deg2rad(float(spec['width_deg']))
        periodic = bool(spec.get('periodic', False))

    elif vtype == 'torsion':
        cv_force = mm.CustomTorsionForce("theta")
        cv_force.addTorsion(int(atoms[0]), int(atoms[1]), int(atoms[2]), int(atoms[3]), [])
        min_val = np.deg2rad(float(spec['min_deg']))
        max_val = np.deg2rad(float(spec['max_deg']))
        width = np.deg2rad(float(spec['width_deg']))
        periodic = bool(spec.get('periodic', True))

    else:
        raise ValueError(f"Unsupported CV type: {vtype}")

    # Dedicated group for bias-energy diagnostics
    cv_force.setForceGroup(int(self._metadynamics_force_group))
    return cv_force, min_val, max_val, width, periodic
```

### 5.2 Build BiasVariables and Metadynamics
Add helper:

```python
def _initialize_metadynamics(self):
    if not self.metadynamics_enabled:
        self._metadynamics = None
        self._metadynamics_variables = []
        self._metadynamics_cv_forces = []
        return

    cfg = self.metadynamics_settings
    self._metadynamics_force_group = int(cfg.get('force_group', 11))

    variables = []
    cv_forces = []
    for spec in cfg['variables']:
        cv_force, min_val, max_val, width, periodic = self._build_metadynamics_cv_force(spec)
        cv_forces.append(cv_force)
        variables.append(
            BiasVariable(
                cv_force,
                min_val,
                max_val,
                width,
                periodic
            )
        )

    self._metadynamics_cv_forces = cv_forces
    self._metadynamics_variables = variables

    # NOTE: verify keyword names with your installed OpenMM version.
    self._metadynamics = Metadynamics(
        self.system,
        variables,
        self.temperature,
        float(cfg['bias_factor']),
        float(cfg['hill_height_kjmol']) * unit.kilojoules_per_mole,
        int(cfg['hill_frequency']),
        saveFrequency=(int(cfg['save_frequency']) if cfg.get('save_frequency') else None),
        biasDir=cfg.get('bias_dir', None),
    )
```

## 6. `run_qmmm()` Integration

### 6.1 Initialize metadynamics before creating `Simulation`
In `run_qmmm()`, after temperature/friction/timestep conversion and before `app.Simulation(...)`:

```python
self.metad_bias_energies = []
self.metad_cv_history = []
self._metadynamics = None

if self.metadynamics_enabled:
    # Optional guard to avoid double-biasing in early rollout
    if self.bias_force_reaction_prop is not None:
        raise RuntimeError(
            "Both legacy bias_force_reaction_prop and metadynamics are enabled. "
            "Disable one to avoid conflicting bias potentials."
        )
    self._initialize_metadynamics()
```

### 6.2 Step using metadynamics driver when enabled
Replace:

```python
self.simulation.step(1)
```

with:

```python
if self._metadynamics is not None:
    self._metadynamics.step(self.simulation, 1)
else:
    self.simulation.step(1)
```

### 6.3 Fix legacy bias `None` unsafe check
Current pattern can fail:

```python
if len(self.bias_force_reaction_prop) != 0 ...
```

Use:

```python
if (
    self.bias_force_reaction_prop is not None
    and len(self.bias_force_reaction_prop) != 0
    and len(self.global_theta_list) > self.last_force_point
    and self.point_checker != 0
    and self.step % self.bias_force_reaction_prop[3] == 0
):
    ...
```

## 7. Energy and Observable Accounting

### 7.1 Extend runtime cache
In `_build_run_qmmm_runtime_cache()`:

```python
cache = {
    'energy_unit': energy_unit,
    'qm_mm_groups': {1, 2},
    'mm_groups': {3, 4, 5, 6, 7},
    'dof': self.system.getNumParticles() * 3 - self.system.getNumConstraints(),
    'has_temp_method': hasattr(self.integrator, 'computeSystemTemperature'),
    'gas_constant': gas_constant,
    'metad_enabled': bool(self._metadynamics is not None),
    'metad_group': {int(self._metadynamics_force_group)} if self._metadynamics is not None else None,
}
return cache
```

### 7.2 Include metad bias in reference and optimized observables
In `_collect_step_observables_reference()`:

```python
meta_bias = 0.0
if runtime_cache['metad_enabled']:
    meta_bias_q = context.getState(
        getEnergy=True,
        groups=runtime_cache['metad_group']
    ).getPotentialEnergy()
    meta_bias = meta_bias_q.value_in_unit(energy_unit)
```

Then:

```python
pot = qm * energy_unit + qm_mm + mm + (meta_bias * energy_unit)
...
return {
    ...
    'metad_bias': float(meta_bias),
    ...
}
```

In `_collect_step_observables_optimized()`:

```python
meta_bias = 0.0
if runtime_cache['metad_enabled']:
    meta_bias = context.getState(
        getEnergy=True,
        groups=runtime_cache['metad_group']
    ).getPotentialEnergy().value_in_unit(energy_unit)

potential = qm + qm_mm + mm + meta_bias
...
return {
    ...
    'metad_bias': float(meta_bias),
    ...
}
```

### 7.3 Append and print
In `_append_step_observables()`:

```python
self.metad_bias_energies.append(observables.get('metad_bias', 0.0))
```

In `_print_run_qmmm_step()`:

```python
if 'metad_bias' in observables:
    print('Metadynamics Bias Energy:', observables['metad_bias'], 'kJ/mol')
```

## 8. CV Value Tracking (Optional but Recommended)

Add helper:

```python
def _collect_metad_cv_values(self, context):
    if not self._metadynamics_cv_forces:
        return []
    vals = []
    for cv_force in self._metadynamics_cv_forces:
        # Evaluate each CV force in isolation using its force group.
        # For scalar CV expressions (r/theta), this returns the CV-linked contribution.
        group = {cv_force.getForceGroup()}
        e = context.getState(getEnergy=True, groups=group).getPotentialEnergy()
        vals.append(float(e.value_in_unit(unit.kilojoules_per_mole)))
    return vals
```

Inside step loop after observable collection:

```python
if self._metadynamics is not None:
    self.metad_cv_history.append(self._collect_metad_cv_values(self.simulation.context))
else:
    self.metad_cv_history.append([])
```

## 9. HDF5 Summary Extension

In `output_file_writer()`, append a `metadynamics` group:

```python
meta_group = h5f.require_group("metadynamics")
meta_group.attrs["enabled"] = int(self._metadynamics is not None)
if self.metadynamics_settings is not None:
    meta_group.attrs["bias_factor"] = float(self.metadynamics_settings.get("bias_factor", 0.0))
    meta_group.attrs["hill_height_kjmol"] = float(self.metadynamics_settings.get("hill_height_kjmol", 0.0))
    meta_group.attrs["hill_frequency"] = int(self.metadynamics_settings.get("hill_frequency", 0))

if "bias_energy_kjmol" not in meta_group:
    ds = meta_group.create_dataset("bias_energy_kjmol", shape=(0,), maxshape=(None,), dtype="float64", chunks=True)
else:
    ds = meta_group["bias_energy_kjmol"]

new_bias = np.asarray(self.metad_bias_energies, dtype=np.float64)
if new_bias.size:
    old = ds.shape[0]
    ds.resize((old + new_bias.size,))
    ds[old:] = new_bias
```

Optional CV history dataset:

```python
if self.metad_cv_history and any(len(x) > 0 for x in self.metad_cv_history):
    ncv = max(len(x) for x in self.metad_cv_history)
    arr = np.full((len(self.metad_cv_history), ncv), np.nan, dtype=np.float64)
    for i, row in enumerate(self.metad_cv_history):
        arr[i, :len(row)] = row

    if "cv_history" not in meta_group:
        cvds = meta_group.create_dataset(
            "cv_history",
            data=arr,
            maxshape=(None, ncv),
            dtype="float64",
            chunks=True
        )
    else:
        cvds = meta_group["cv_history"]
        old = cvds.shape[0]
        cvds.resize((old + arr.shape[0], cvds.shape[1]))
        cvds[old:, :] = arr
```

## 10. Upstream Wiring in `imforcefieldgenerator.py`

Where `self.dynamics_settings` dict is created, add:

```python
'metadynamics': self.metadynamics_settings
```

and define a default in generator `__init__`:

```python
self.metadynamics_settings = None
```

## 11. MPI Considerations

- In current root-worker mode, only root rank constructs `Simulation` and performs integration.
- Worker ranks only execute interpolation compute in `_run_qmmm_worker_service()`.
- Therefore metadynamics object must remain root-only (which naturally follows current control flow).
- Do not introduce any metadynamics calls into worker service loop.

## 12. Validation Strategy

### 12.1 Baseline parity (metadynamics disabled)
- Run a short QM/MM trajectory with old settings.
- Confirm no change in `qm_potentials`, `mm_potentials`, `total_energies` except normal floating-point jitter.

### 12.2 1-CV smoke test
- Use one torsion CV.
- Validate:
  - metad bias energy increases over time,
  - sampled torsion range expands compared to unbiased run.

### 12.3 Stability checks
- Ensure no NaN in CV ranges and no out-of-range atom indices.
- Ensure restart with same `bias_dir` continues bias grid accumulation as expected.

### 12.4 MPI smoke test
- Two-rank run in root-worker mode.
- Verify no `MPI control-channel desync` and no worker deadlock.

## 13. Suggested Delivery Order

1. Add config parsing + validation + runtime members.
2. Add CV/metadynamics helper builders.
3. Integrate run-step switching (`metadynamics.step` vs `simulation.step`).
4. Add observable and HDF5 bookkeeping.
5. Add/adjust tests and one example input block.

## 14. Practical Notes

- OpenMM API details can vary by version; confirm keyword names for `Metadynamics(...)` in your installed OpenMM build.
- Keep metadynamics and legacy `bias_force_reaction_prop` mutually exclusive in first merge to reduce debugging surface.
- Preserve force groups 0-10 already used by QM/MM terms; default metadynamics to group 11.

