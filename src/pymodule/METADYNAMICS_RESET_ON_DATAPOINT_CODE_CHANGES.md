# Metadynamics Bias Reset On Datapoint Addition (Soft-Only)

This file lists the concrete code changes needed in
`imdatabasepointcollecter.py` so that when a new datapoint is added, the
metadynamics bias is reset and growth starts again from zero.

## Why `bias_dir` is needed (and when it is not)

`bias_dir` is only needed if you want metadynamics bias persistence/sharing:

- Restart continuation: keep hills on disk and continue from previous runs.
- Multi-walker sharing: multiple jobs read/write bias files in the same folder.

If you do not need either, leave `bias_dir=None` and `save_frequency=None`.
That is the cleanest setup for reset-on-datapoint behavior.

Important for your reset workflow:

- Soft reset clears in-memory hills.
- If old hill files remain in `bias_dir`, OpenMM can load them again later.
- Therefore, for a true fresh restart after each datapoint:
  set `clear_bias_dir_on_reset=True` when using `bias_dir`.

## 1) Add New State Fields In `__init__`

Insert the following block in `__init__`, right after existing metadynamics fields (`self.metadynamics_bias_energies`, `self.metadynamics_cv_history`):

```python
        # metadynamics reset controls (soft reset only)
        self.metadynamics_reset_on_datapoint_add = False
        self.metadynamics_clear_bias_dir_on_reset = False
        self._pending_metadynamics_reset = False
        self.metadynamics_reset_count = 0
```

## 2) Extend Metadynamics Settings Validation

In `_validate_metadynamics_settings`, add the following checks near the end of the function (after current `hill_frequency` validation):

```python
        clear_bias_dir = bool(cfg.get('clear_bias_dir_on_reset', False))
        if clear_bias_dir and cfg.get('bias_dir', None) is None:
            raise ValueError(
                "metadynamics.clear_bias_dir_on_reset=True requires metadynamics.bias_dir to be set."
            )
```

## 3) Parse New Settings In `update_settings`

In `update_settings`, right after:

```python
        self.metadynamics_enabled = bool(self.metadynamics_settings is not None and
                                         self.metadynamics_settings.get('enabled', False))
```

add:

```python
        if self.metadynamics_settings is not None:
            self.metadynamics_reset_on_datapoint_add = bool(
                self.metadynamics_settings.get('reset_on_datapoint_add', False)
            )
            self.metadynamics_clear_bias_dir_on_reset = bool(
                self.metadynamics_settings.get('clear_bias_dir_on_reset', False)
            )
        else:
            self.metadynamics_reset_on_datapoint_add = False
            self.metadynamics_clear_bias_dir_on_reset = False
```

## 4) Add Helper Methods For Soft Bias Reset

Add the following methods near `_initialize_metadynamics` (same class scope):

```python
    def _clear_metadynamics_bias_dir(self):
        if self.metadynamics_settings is None:
            return

        bias_dir = self.metadynamics_settings.get('bias_dir', None)
        if not bias_dir:
            return

        if not os.path.isdir(bias_dir):
            return

        for fname in os.listdir(bias_dir):
            if not (fname.startswith('bias_') or fname.startswith('temp_')):
                continue
            if not fname.endswith('.npy'):
                continue
            fpath = os.path.join(bias_dir, fname)
            try:
                os.remove(fpath)
            except OSError:
                # Non-fatal: keep simulation running even if cleanup is partial
                pass

    def _reset_metadynamics_bias(self, context, reason='datapoint_added'):
        """
        Soft reset only:
        clear accumulated metadynamics hills and restart bias growth from zero.
        """
        if self._metadynamics is None or not self.metadynamics_enabled:
            return False

        if self.metadynamics_clear_bias_dir_on_reset:
            self._clear_metadynamics_bias_dir()

        mtd = self._metadynamics
        mtd._selfBias.fill(0.0)
        mtd._totalBias.fill(0.0)
        mtd._loadedBiases = {}

        if len(mtd.variables) == 1:
            mtd._table.setFunctionParameters(mtd._totalBias.flatten(), *mtd._limits)
        else:
            mtd._table.setFunctionParameters(
                *mtd._widths,
                mtd._totalBias.flatten(),
                *mtd._limits
            )
        mtd._force.updateParametersInContext(context)

        # Reset local history containers as well.
        self.metad_bias_energies = []
        self.metad_cv_history = []
        self.metadynamics_bias_energies = []
        self.metadynamics_cv_history = []
        self.metadynamics_reset_count += 1

        print(
            f"[metadynamics] soft bias reset done (reason={reason}, "
            f"count={self.metadynamics_reset_count})"
        )
        return True
```

Note: this does not remove the metadynamics force object from the `System`.
It removes the accumulated bias surface (hills), which is what you want for
restart-from-zero behavior without context reinitialization.

## 5) Mark Reset As Pending After Successful `add_point(...)`

In `point_correlation_check`, after each successful `self.add_point(...)` call (there are two locations), add:

```python
                if self.metadynamics_enabled and self.metadynamics_reset_on_datapoint_add:
                    self._pending_metadynamics_reset = True
```

Applied in both branches:

```python
                self.add_point(state_specific_molecules, self.non_core_symmetry_groups)
                if self.metadynamics_enabled and self.metadynamics_reset_on_datapoint_add:
                    self._pending_metadynamics_reset = True
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0
```

and:

```python
                self.add_point(state_specific_molecules, self.non_core_symmetry_groups)
                if self.metadynamics_enabled and self.metadynamics_reset_on_datapoint_add:
                    self._pending_metadynamics_reset = True
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0
```

## 6) Execute Reset In Existing Dynamics Reset Branch

In `update_forces`, inside:

```python
            if self.point_checker == 0:
```

add this at the top of that block:

```python
                if self._pending_metadynamics_reset and self._metadynamics is not None:
                    self._reset_metadynamics_bias(context, reason='datapoint_added')
                    self._pending_metadynamics_reset = False
```

That keeps reset timing aligned with your current logic where coordinates and velocities are already reset after new datapoints.

## 7) Optional: Clear Pending Flag At Run Start

In `run_qmmm`, right before entering the main timestep loop, initialize:

```python
        self._pending_metadynamics_reset = False
        self.metadynamics_reset_count = 0
```

This avoids carrying any stale flags across restarts or repeated calls.

## 8) Example Input Settings (Soft-Only)

Use this in your dynamics input to activate behavior:

```python
"metadynamics": {
    "enabled": True,
    "variables": [
        {
            "type": "torsion",
            "atoms": [3, 4, 7, 9],
            "min_deg": -180.0,
            "max_deg": 180.0,
            "width_deg": 10.0,
            "periodic": True
        }
    ],
    "bias_factor": 10.0,
    "hill_height_kjmol": 1.2,
    "hill_frequency": 100,
    "save_frequency": 1000,
    "bias_dir": "meta_bias",

    "reset_on_datapoint_add": True,
    "clear_bias_dir_on_reset": True
}
```

If you do not need persisted hills:

```python
"metadynamics": {
    "enabled": True,
    "...": "...",
    "save_frequency": None,
    "bias_dir": None,
    "reset_on_datapoint_add": True
}
```

## 9) Notes On MPI Behavior

- In root-worker mode, this reset runs on the root rank in the same branch where the dynamics is reset (`point_checker == 0`), so it does not add new control messages.
- Because this is soft-only, there is no `context.reinitialize()` and no extra force graph mutation.
