# 1. Executive Summary

## Purpose of the Workflow
The `IMForceFieldGenerator`-centered workflow builds interpolation databases for QM/MM dynamics by:
- generating initial QM datapoints,
- running adaptive QM/MM trajectories,
- adding new datapoints when interpolation error exceeds thresholds,
- storing all model data in HDF5 for runtime interpolation (`InterpolationDriver`).

Core participants:
- `IMForceFieldGenerator`: construction orchestrator.
- `IMDatabasePointCollecter`: QM/MM runtime + adaptive datapoint growth.
- `InterpolationDriver`: interpolation kernel (energy/gradient) with symmetry/cluster modes.
- `InterpolationDatapoint`: coordinate transforms + persisted datapoint schema.

## Current Strengths
- Rich feature set: symmetry-expanded datapoints, rotor-cluster banking, trust-radius optimization, MPI preload path, sampling screen, runtime profiling hooks.
- Scientifically meaningful interpolation core: local Taylor expansion with normalized Shepard-style blending and derivative-consistent gradient assembly.
- HDF5 persistence schema is relatively complete and includes both scientific arrays and cluster metadata.
- A non-trivial integration test exists (`tests/integration/interpolation_test.py`) covering an end-to-end ground-state workflow.

## Main Merge-Readiness Blockers
The implementation is currently blocked primarily by correctness and maintainability issues:
1. Multiple confirmed logic bugs in control/configuration paths (some critical).
2. Multi-state scaffolding is present but effectively non-functional/inconsistent.
3. Severe duplication between `IMForceFieldGenerator` and `IMDatabasePointCollecter` for datapoint creation and rotor-cluster registry logic.
4. Extremely large functions with hidden side effects and poor workflow boundaries.
5. Inconsistent configuration schemas and typo-driven mode drift.

## Highest-Priority Changes (Pre-Merge)
1. Fix critical correctness bugs:
   - `all(roots_to_follow) == 0` assertion logic.
   - desired-density gate in `IMForceFieldGenerator.compute`.
   - Bayesian trust-radius call/signature mismatch.
   - sampling cluster-bank load using only the last `root`.
   - `InterpolationDriver.calculate_translation_coordinates` variable bug.
2. Remove/guard dead or broken runtime branches (`simple_interpolation` path).
3. Unify duplicated datapoint/cluster registry logic into shared utilities.
4. Stabilize configuration contracts (`use_opt_confidence_radius`, sampling keys, state selection semantics).
5. Split mega-functions into explicit stages and add minimal API-level documentation/tests for each stage.

---

# 2. High-Level Architecture

## Class Responsibilities

### `IMForceFieldGenerator`
Primary orchestration layer.
- Prepares driver bundles (`gs`/`es`) and defaults.
- Builds Z-matrix and symmetry metadata (`set_up_the_system`).
- Seeds database (`add_point` / `add_point_rotor`) before trajectory expansion.
- Defines scan structures and target datapoint densities.
- Spawns one `IMDatabasePointCollecter` per scan structure and runs adaptive QM/MM construction.

### `IMDatabasePointCollecter`
Adaptive trajectory-expansion runtime.
- Builds OpenMM system and QM/MM partition (`system_from_molecule`).
- Loads existing HDF5 interpolation datapoints into `InterpolationDriver` objects (per root).
- Executes QM/MM dynamics (`run_qmmm`) and updates forces each step using interpolation.
- Decides whether to trigger expensive QM verification and datapoint insertion (`point_correlation_check`).
- Adds new core/symmetry/cluster datapoints and updates runtime caches.

### `InterpolationDriver`
Runtime interpolation engine.
- Reads datapoints from HDF5.
- Computes interpolated energy/gradient for a geometry.
- Supports symmetry-expanded weighting and rotor-cluster assembly.
- Supports Cartesian/internal distance metrics and trust-radius derivatives.

### `InterpolationDatapoint`
Datapoint data model + coordinate transforms.
- Stores Cartesian/internal quantities.
- Computes B/B2 and transforms gradient/Hessian to internal space.
- Handles bond/dihedral coordinate transforms.
- Serializes/deserializes complete datapoint payload to/from HDF5.

## Interaction Topology
`IMForceFieldGenerator.compute` is the top-level entry point. It configures per-root interpolation settings and calls `IMDatabasePointCollecter.run_qmmm`, which repeatedly calls `InterpolationDriver.compute` and occasionally `IMDatabasePointCollecter.add_point*` to grow the HDF5 database. `InterpolationDatapoint` is the storage/transform primitive used by both the generator and collector.

## Where Separation Is Good
- Interpolation math is mostly localized to `InterpolationDriver`/`InterpolationDatapoint`.
- OpenMM simulation machinery is mostly localized to `IMDatabasePointCollecter`.
- HDF5 datapoint schema is centralized in `InterpolationDatapoint`.

## Where Responsibilities Are Entangled
- `IMForceFieldGenerator` duplicates substantial datapoint-creation and rotor-bank serialization logic that also exists in `IMDatabasePointCollecter`.
- Both generator and collector perform near-identical QM-evaluate/transform/write datapoint flows.
- Multi-state policy is partially in generator, partially in collector, partially implicit in list index conventions.
- Symmetry metadata is passed around as fragile positional lists instead of typed structures.

## Intended Lifecycle of a Construction Run
1. User builds `IMForceFieldGenerator` with drivers and options.
2. `compute(...)` sets output files/settings and calls `set_up_the_system(...)`.
3. Initial datapoints are added via `add_point` or `add_point_rotor`.
4. Candidate structures are generated along configured dihedral scans.
5. For each candidate structure:
   - instantiate/configure `IMDatabasePointCollecter`,
   - initialize OpenMM/QM/MM and interpolation drivers,
   - run QM/MM (`run_qmmm`) with adaptive point checks.
6. New points are written to root-specific HDF5 databases.
7. Optional trust-radius optimization and sampling DB updates occur.
8. Summary/trace artifacts are written.

---

# 3. End-to-End Workflow Reconstruction

## Step-by-Step Runtime Flow (Actual Code Path)

### Stage A: `IMForceFieldGenerator.compute`
Location: `imforcefieldgenerator.py:1316+`

1. Test hook plumbing is installed.
2. Database filenames are resolved (`im_database_{root}.h5`, sampling mirrors).
3. `set_up_the_system(molecule, extract_z_matrix=...)` is called.
4. Single-root branch (currently dominant) prepares:
   - `states_interpolation_settings[root]`,
   - `dynamics_settings` payload for collector.
5. If no database exists, seed datapoints are generated:
   - optionally from minimized conformers,
   - then `add_point(...)` (or cluster mode path).
6. Build scan structures with `determine_molecules_along_dihedral_scan`.
7. For each scan structure, instantiate `IMDatabasePointCollecter` and run `run_qmmm`.

Observed branch hazard:
- Multi-root path (`len(self.roots_to_follow) > 1`) currently does not execute the main setup/run logic that exists in the `else` block.

### Stage B: `set_up_the_system`
Location: `imforcefieldgenerator.py:774+`

1. Build/restore root-specific Z-matrix.
2. Build CH3-focused symmetry groups.
3. Build rotatable bond definitions and promote selected ring torsions to impropers.
4. Build packed `symmetry_information` per root:
   - atom map,
   - rot groups,
   - core/non-core atoms,
   - rotatable bond indices,
   - dihedral map and partition indices.
5. Build `self.symmetry_rotors` for rotor-cluster logic.
6. Optionally generate conformers and transition-like structures.

### Stage C: Initial Datapoint Seeding (`add_point` / `add_point_rotor`)
Locations:
- `imforcefieldgenerator.py:3501+` (`add_point`)
- `imforcefieldgenerator.py:2792+` (`add_point_rotor`)

Flow:
1. Build per-state molecule batches (`gs` / `es` entries).
2. Optional symmetry rotations and constrained optimization.
3. Evaluate E/G/H with configured drivers.
4. Mass-weight gradients/Hessians.
5. Build `InterpolationDatapoint` objects.
6. Assign labels, family labels, bank roles (`core`, `symmetry`, `cluster`).
7. Transform to internal-space tensors and write to HDF5.
8. For cluster mode, also write rotor-cluster registry JSON blocks.

### Stage D: `IMDatabasePointCollecter.run_qmmm`
Location: `imdatabasepointcollecter.py:2038+`

1. Validate system and initialize OpenMM simulation.
2. Build per-root interpolation drivers (`InterpolationDriver`) and load HDF5 datapoints.
3. Optionally initialize sampling interpolation drivers.
4. Main step loop:
   - `update_forces(context)`:
     - extract QM coordinates,
     - `update_gradient_and_energy` -> interpolation for all roots,
     - select current state,
     - run point-add heuristics and optional `point_correlation_check`.
   - integrate 1 step.
   - stop based on density/unadded cycles or step limits.
5. Emit summaries and write output artifacts.

### Stage E: Adaptive Expansion (`point_correlation_check`)
Location: `imdatabasepointcollecter.py:3650+`

1. Optional cheap sampling screen (`_sampling_screen_check`).
2. Compute QM energies/gradients for relevant roots.
3. Compare against interpolation predictions (energy per atom + gradient RMSD + cosine alignment).
4. Decide to add point or store as allowed molecule.
5. If configured, optimize trust radii (`determine_trust_radius_gradient` or Bayesian path).
6. If adding point:
   - optionally identify important internal coordinates,
   - constrained optimize,
   - call `add_point`/`add_point_rotor`,
   - refresh interpolation caches.

## Data and State Propagation
- Root-indexed dictionaries are central (`root -> datapoints/settings/drivers`).
- HDF5 files are both persistence and active communication medium between stages.
- `symmetry_information` is passed as positional array-like tuples, causing implicit coupling and fragile indexing.
- Label conventions (`point_N`, `_symmetry_k`, `_core`, `_cluster_*`) drive many grouping operations.

## Workflow Complexity / Clarity Problems
1. Huge methods hide stage boundaries (`compute`, `run_qmmm`, `add_point*`, `point_correlation_check`).
2. Control flags have mixed types and overloaded meanings.
3. Multi-state logic is present but partly stubbed/implicit.
4. Duplication causes diverging behavior (e.g., cluster threshold mismatch).
5. Runtime behavior depends on many mutable object attributes set across classes.

---

# 4. Merge-Blocking Issues

## Issue 1: Ground-State Restriction Assertion Is Logically Wrong
- Severity: `critical`
- Affected: `imforcefieldgenerator.py::__init__`, `imforcefieldgenerator.py::compute`
- Location:
  - `imforcefieldgenerator.py:218`
  - `imforcefieldgenerator.py:1411`
- Problem:
  - Uses `all(roots_to_follow) == 0` and `all(self.roots_to_follow) == 0`.
  - For `[0,1]`, `all(...)` is `False`, and `False == 0` evaluates to `True`.
- Why it matters:
  - Multi-root unsupported cases can pass assertion gates silently.
  - Scientific mode gating becomes unreliable.
- Proposed resolution:
  - Replace with explicit condition: `all(int(r) == 0 for r in roots_to_follow)`.

## Issue 2: Multi-Root Branch in `compute` Is Effectively Incomplete
- Severity: `critical`
- Affected: `imforcefieldgenerator.py::compute`
- Location: `imforcefieldgenerator.py:1408-1417` and subsequent large `else` block.
- Problem:
  - The single-root execution body is under the `else`; the `len>1` branch does not perform equivalent setup/run.
- Why it matters:
  - Multi-root runs can appear accepted but not execute intended workflow.
- Proposed resolution:
  - Either hard-fail multi-root immediately pre-merge, or implement full branch parity.

## Issue 3: Desired Datapoint Density Gate Is Broken
- Severity: `critical`
- Affected: `imforcefieldgenerator.py::compute`
- Location: `imforcefieldgenerator.py:1756-1773`
- Problem:
  - `desired_density` is never set `True`.
  - typo variable `desiered_point_density` changes type from `int` to `bool`.
- Why it matters:
  - Structural scans continue even when target density is reached; expensive QM/MM overrun.
- Proposed resolution:
  - Track threshold in a separate boolean (`reached_target_density`).

## Issue 4: Bayesian Trust-Radius Path Has Signature Mismatches
- Severity: `critical`
- Affected:
  - `imdatabasepointcollecter.py::point_correlation_check`
  - `imdatabasepointcollecter.py::determine_beysian_trust_radius`
- Location:
  - call: `imdatabasepointcollecter.py:3845-3850`
  - definition: `imdatabasepointcollecter.py:5491`
  - internal call: `imdatabasepointcollecter.py:5499`
- Problem:
  - Caller omits required `z_matrix` argument.
  - Function calls `calculate_distance_to_ref` without required `symmetry_info` argument.
- Why it matters:
  - Bayesian mode is non-functional and will fail at runtime.
- Proposed resolution:
  - Fix call signature and pass `sym_dict`/`z_matrix` consistently.

## Issue 5: Sampling Rotor Cluster Banks Loaded Only for Last Root
- Severity: `high`
- Affected: `imdatabasepointcollecter.py::run_qmmm`
- Location: `imdatabasepointcollecter.py:2261-2263`
- Problem:
  - Sampling load occurs after `for root in self.roots_to_follow`, using loop variable `root` from last iteration.
- Why it matters:
  - Multi-root sampling uses incomplete cluster bank state.
- Proposed resolution:
  - Move sampling load inside root loop or loop explicitly in sampling block.

## Issue 6: `states_basis` Argument Is Ignored
- Severity: `high`
- Affected: `imforcefieldgenerator.py::compute`
- Location: `imforcefieldgenerator.py:1316, 1357`
- Problem:
  - Method accepts `states_basis` but unconditionally overwrites it.
- Why it matters:
  - API contract is misleading; user-provided basis routing is silently ignored.
- Proposed resolution:
  - Respect input when provided; fallback to defaults only when `None`.

## Issue 7: Default `use_minimized_structures` Is a Tuple by Accident
- Severity: `high`
- Affected: `imforcefieldgenerator.py::__init__`
- Location: `imforcefieldgenerator.py:387`
- Problem:
  - Trailing comma creates tuple: `([True, [], []],)`.
- Why it matters:
  - Later list indexing semantics are brittle and can fail unexpectedly.
- Proposed resolution:
  - Remove trailing comma.

## Issue 8: Missing `simple_interpolation` Implementation
- Severity: `high`
- Affected: `interpolationdriver.py::compute`
- Location: `interpolationdriver.py:1099-1102`
- Problem:
  - Branch calls `self.simple_interpolation()` but no such method exists.
- Why it matters:
  - Selecting `interpolation_type='simple'` crashes immediately.
- Proposed resolution:
  - Implement method or explicitly disallow `simple` with clear error.

## Issue 9: `InterpolationDriver.calculate_translation_coordinates` Uses Undefined Variable
- Severity: `high`
- Affected: `interpolationdriver.py::calculate_translation_coordinates`
- Location: `interpolationdriver.py:2056-2061`
- Problem:
  - `w = np.asarray(w, ...)` references undefined `w`; should be `weights`.
- Why it matters:
  - Weighted-centering mode is broken.
- Proposed resolution:
  - Replace with `w = np.asarray(weights, dtype=float)` and validate shape/sum.

## Issue 10: Divergent Rotor Cluster Coupling Thresholds
- Severity: `high`
- Affected:
  - `imforcefieldgenerator.py::add_point_rotor`
  - `imdatabasepointcollecter.py::add_point_rotor`
- Location:
  - `imforcefieldgenerator.py:3246` (`threshold=0.01`)
  - `imdatabasepointcollecter.py:4627` (`threshold=0.04`)
- Problem:
  - Same conceptual clustering step uses different hard-coded thresholds.
- Why it matters:
  - Initial vs runtime-added cluster families can diverge for identical chemistry.
- Proposed resolution:
  - Centralize threshold as shared constant/config key.

## Issue 11: `use_opt_confidence_radius` Contract Is Inconsistent and Typo-Prone
- Severity: `high`
- Affected:
  - `imforcefieldgenerator.py::__init__`
  - `imdatabasepointcollecter.py::__init__`, `point_correlation_check`, `determine_trust_radius_gradient`
- Location:
  - generator default typo: `imforcefieldgenerator.py:392` (`'mutli_grad'`)
  - collector default bool: `imdatabasepointcollecter.py:381`
- Problem:
  - Sometimes treated as bool, elsewhere as list with positional indices.
  - Mode string typo in default (`mutli_grad` vs `multi_grad`).
- Why it matters:
  - Runtime branching can silently skip optimization paths or crash in non-generator usage.
- Proposed resolution:
  - Replace with typed dict/dataclass (`enabled`, `mode`, `initial_radius`, `energy_weight`).

## Issue 12: File Existence Checks Use `os.listdir(os.getcwd())`
- Severity: `medium`
- Affected:
  - `imforcefieldgenerator.py`
  - `imdatabasepointcollecter.py`
- Location:
  - examples: `imforcefieldgenerator.py:2006, 3081, 3768`; `imdatabasepointcollecter.py:4424`
- Problem:
  - Existence logic fails for absolute paths and non-CWD paths.
- Why it matters:
  - Incorrect label loading and datapoint counts depending on execution directory.
- Proposed resolution:
  - Replace with `os.path.exists(path)`.

## Issue 13: `run_qmmm` Return Value Mismatch with Caller Expectation
- Severity: `medium`
- Affected:
  - caller: `imforcefieldgenerator.py::compute`
  - callee: `imdatabasepointcollecter.py::run_qmmm`
- Location:
  - caller assignment: `imforcefieldgenerator.py:1788`
  - callee return: `imdatabasepointcollecter.py:2435`
- Problem:
  - Caller stores `stats`, but callee always returns `None`.
- Why it matters:
  - Misleading telemetry and dead integration channel.
- Proposed resolution:
  - Either return a structured stats object or remove expectation.

## Issue 14: State Selection Is Hardcoded to First Root
- Severity: `high`
- Affected: `imdatabasepointcollecter.py::_postprocess_root_state_selection`
- Location: `imdatabasepointcollecter.py:3341-3344`
- Problem:
  - `self.current_state = self.roots_to_follow[0]` each update.
- Why it matters:
  - Multi-state dynamics machinery does not actually perform state selection.
- Proposed resolution:
  - Implement explicit selection policy (adiabatic/diabatic/min-energy or configured rule).

## Issue 15: Root Variable Leakage Causes Wrong Symmetry Branch Choice
- Severity: `medium`
- Affected: `imdatabasepointcollecter.py::point_correlation_check`
- Location: `imdatabasepointcollecter.py:3819`
- Problem:
  - Condition uses `root` after loop scope; resolves to last iterated root.
- Why it matters:
  - Wrong symmetry dictionary can be used for trust-radius optimization.
- Proposed resolution:
  - Use `self.current_state` directly in condition.

## Issue 16: Driver-Type Check Typo in Excited-State Branch
- Severity: `medium`
- Affected: `imdatabasepointcollecter.py::add_point`
- Location: `imdatabasepointcollecter.py:4982`
- Problem:
  - `isinstance(drivers, TdaEigenSolver)` should be `drivers[0]`.
- Why it matters:
  - Excited-state derivative index routing may be skipped incorrectly.
- Proposed resolution:
  - Fix `isinstance(drivers[0], TdaEigenSolver)`.

---

# 5. Detailed Refactoring Recommendations

## A. Architecture

### Recommendation A1: Extract Shared Datapoint/Registry Utilities
- Current problem:
  - `add_point*`, `_write_cluster_registry_for_family`, `_read_cluster_registry_for_family`, `_load_rotor_cluster_bank_for_root`, `_build_opt_constraint_list` duplicated in both generator and collector.
- Target design:
  - New shared module (for example `imdb_construction_shared.py`) containing pure helpers.
- Expected benefit:
  - Eliminates divergence bugs and reduces maintenance load.
- Risks:
  - Refactor touches many call sites; do staged extraction with parity tests.

### Recommendation A2: Split Orchestration from Execution
- Current problem:
  - `IMForceFieldGenerator.compute` and `IMDatabasePointCollecter.run_qmmm` are monolithic state machines.
- Target design:
  - Extract explicit stages (setup, seed, scan loop, run step, point decision, writeback).
- Expected benefit:
  - Clear workflow boundaries and testable units.
- Risks:
  - Must preserve side-effect order (especially MPI and HDF5 I/O).

## B. API Clarity

### Recommendation B1: Replace Positional Config Lists with Typed Configs
- Current problem:
  - `use_opt_confidence_radius` is positional and inconsistent (`bool` vs list).
- Target design:
  - Use dataclass or dict with named fields.
- Expected benefit:
  - Eliminates index/typo bugs, improves readability.
- Risks:
  - Requires migration at every call site.

### Recommendation B2: Clarify Public API Contracts
- Current problem:
  - `states_basis` argument ignored.
  - `run_qmmm` caller expects stats but receives `None`.
- Target design:
  - Harmonize signatures and returns.
- Expected benefit:
  - More predictable integration behavior.
- Risks:
  - May affect existing scripts relying on implicit behavior.

## C. Naming and Terminology

### Recommendation C1: Normalize Names and Fix Typos
- Current problem:
  - `IMDatabasePointCollecter` misspelling, `determine_beysian_trust_radius`, `identfy_*`, `mutli_grad`, `desiered_*`.
- Target design:
  - Keep backward-compatible aliases but migrate internals/docs to corrected names.
- Expected benefit:
  - Reduces cognitive friction and typo-driven logic errors.
- Risks:
  - Potential breakage if external code imports renamed members directly; provide deprecation window.

## D. Configuration Handling

### Recommendation D1: Central Validation at Entry
- Current problem:
  - Hard-coded defaults spread across classes; key mismatches (`g_thrsh_...` vs `g_rmsd_...`).
- Target design:
  - Validate and normalize `dynamics_settings` and interpolation settings once.
- Expected benefit:
  - Immediate error reporting; fewer latent runtime failures.
- Risks:
  - Existing flexible scripts may need updates.

## E. Data Flow and State

### Recommendation E1: Replace `symmetry_information` Positional Tuple with Structured Object
- Current problem:
  - Access via magic indices (`[3]`, `[4]`, `[7]`...) across files.
- Target design:
  - Dataclass with named fields (`core_atoms`, `excluded_atoms`, `rot_bonds`, `dihedral_map`, ...).
- Expected benefit:
  - Safer refactors and fewer index mistakes.
- Risks:
  - Medium migration effort.

### Recommendation E2: Make State-Selection Policy Explicit
- Current problem:
  - `current_state` reset to first root each step.
- Target design:
  - Introduce policy method (`select_current_state(...)`) and unit tests.
- Expected benefit:
  - Multi-state support becomes real and testable.
- Risks:
  - Scientific behavior change if policy differs from implicit old behavior.

## F. Error Handling and Logging

### Recommendation F1: Reduce Unstructured Prints
- Current problem:
  - Large volume of `print(...)` debug traces in production path.
- Target design:
  - Use structured logging levels and per-rank control.
- Expected benefit:
  - Better diagnostics, cleaner MPI runs.
- Risks:
  - Need log-level defaults aligned with current user expectations.

## G. Scientific Transparency

### Recommendation G1: Document Implemented Equations Adjacent to Code
- Current problem:
  - Core formulas are present but spread out; intent not consistently stated.
- Target design:
  - Add concise docstrings/equations for Taylor model, weight normalization, trust-radius objectives.
- Expected benefit:
  - Easier peer review and scientific reproducibility.
- Risks:
  - Low.

## H. Testability

### Recommendation H1: Add Focused Regression Tests for Merge Blockers
- Current problem:
  - Integration test exists but does not systematically cover all risky branches.
- Target design:
  - Add small tests for each bugfix path:
    - assertion semantics,
    - density gate,
    - bayesian mode call path,
    - sampling root loading,
    - simple interpolation guard,
    - config key normalization.
- Expected benefit:
  - Prevents immediate regressions during cleanup.
- Risks:
  - Low.

---

# 6. Drop-In Code Snippets and Replacement Instructions

## Snippet 1: Fix Ground-State Assertion Logic
Type: direct replacement

Where:
- `imforcefieldgenerator.py` at `__init__` and `compute` assertions.

Replace with:
```python
assert_msg_critical(
    all(int(root) == 0 for root in roots_to_follow),
    "The current version is restricted to ground-state potentials. Later version will allow multi-state construction!"
)
```
And similarly in `compute`:
```python
assert_msg_critical(
    all(int(root) == 0 for root in self.roots_to_follow),
    "Only ground-state potential construction is currently supported! Will be updated to multi-state in the near future!"
)
```
Follow-up edits:
- If multi-state remains unsupported, fail early regardless of list length to avoid silent partial execution.

## Snippet 2: Fix `use_minimized_structures` Default Tuple Bug
Type: direct replacement

Where:
- `imforcefieldgenerator.py:__init__`

Replace:
```python
self.use_minimized_structures = [True, [], []],
```
With:
```python
self.use_minimized_structures = [True, [], []]
```

## Snippet 3: Respect `states_basis` Argument in `compute`
Type: direct replacement

Where:
- `imforcefieldgenerator.py::compute`

Replace:
```python
states_basis = {'gs':self.gs_basis_set_label, 'es':self.es_basis_set_label}
```
With:
```python
if states_basis is None:
    states_basis = {'gs': self.gs_basis_set_label, 'es': self.es_basis_set_label}
else:
    states_basis = {
        'gs': states_basis.get('gs', self.gs_basis_set_label),
        'es': states_basis.get('es', self.es_basis_set_label),
    }
```

## Snippet 4: Fix Desired-Density Gate
Type: direct replacement

Where:
- `imforcefieldgenerator.py::compute` in scan loop.

Replace block around `desired_density/desiered_point_density` with:
```python
desired_point_density = int(self.dynamics_settings['desired_datapoint_density'])
reached_target_density = False
current_structure_density = {}

for root in density_of_datapoints.keys():
    value = density_of_datapoints[root][key][current_dihedral_angle]
    current_structure_density[root] = value
    if value >= desired_point_density:
        reached_target_density = True

if not reached_target_density:
    # existing run_qmmm body
    ...
```
Follow-up edits:
- Remove misspelled variable `desiered_point_density` throughout this function.

## Snippet 5: Fix Sampling Rotor-Bank Load for All Roots
Type: helper addition / local restructure

Where:
- `imdatabasepointcollecter.py::run_qmmm`

Replace:
```python
if self.sampling_enabled:
    self.qm_sampling_rotor_cluster_banks[root] = self._load_rotor_cluster_bank_for_root(root, sampling_mode=True)
    self._initialize_sampling_impes_drivers(inv_sqrt_masses)
```
With:
```python
if self.sampling_enabled:
    for root in self.roots_to_follow:
        self.qm_sampling_rotor_cluster_banks[root] = self._load_rotor_cluster_bank_for_root(
            root, sampling_mode=True
        )
    self._initialize_sampling_impes_drivers(inv_sqrt_masses)
```

## Snippet 6: Fix Bayesian Trust-Radius Signature Usage
Type: direct replacement + interface cleanup

Where:
- call site: `imdatabasepointcollecter.py::point_correlation_check`
- callee: `imdatabasepointcollecter.py::determine_beysian_trust_radius`

Call site replacement:
```python
trust_radius = self.determine_beysian_trust_radius(
    self.allowed_molecules[self.current_state]['molecules'],
    self.allowed_molecules[self.current_state]['qm_energies'],
    self.qm_data_point_dict[self.current_state],
    self.interpolation_settings[self.current_state],
    self.qm_symmetry_datapoint_dict[self.current_state],
    sym_dict,
    self.root_z_matrix[self.current_state],
)
```
Inside function replacement:
```python
_, distance, _ = self.calculate_distance_to_ref(
    mol.get_coordinates_in_bohr(),
    dp.cartesian_coordinates,
    sym_dict,
)
```

## Snippet 7: Fix Weighted-Centering Variable Bug
Type: direct replacement

Where:
- `interpolationdriver.py::calculate_translation_coordinates`

Replace:
```python
if weights is not None:
    w = np.asarray(w, dtype=float)
    W = np.sum(w)
    return given_coordinates - np.sum(given_coordinates * w[:, None], axis=0) / W
```
With:
```python
if weights is not None:
    w = np.asarray(weights, dtype=float).reshape(-1)
    if w.shape[0] != given_coordinates.shape[0]:
        raise ValueError("weights length must match number of coordinates")
    W = float(np.sum(w))
    if abs(W) < 1.0e-15:
        raise ValueError("sum(weights) must be non-zero")
    return given_coordinates - np.sum(given_coordinates * w[:, None], axis=0) / W
```

## Snippet 8: Guard Unsupported `simple` Interpolation Mode
Type: direct replacement

Where:
- `interpolationdriver.py::compute`

Replace simple-branch body with:
```python
if self.interpolation_type == 'simple':
    raise NotImplementedError(
        "InterpolationDriver: interpolation_type='simple' is currently not implemented. "
        "Use 'shepard'."
    )
```

## Snippet 9: Replace CWD-Listdir File Checks
Type: direct replacement

Where:
- all occurrences like `if target_file in os.listdir(os.getcwd()):`

Replace with:
```python
if os.path.exists(target_file):
```
Follow-up edits:
- keep behavior consistent for relative and absolute paths.

## Snippet 10: Fix Excited-State Driver Type Check Typo
Type: direct replacement

Where:
- `imdatabasepointcollecter.py:4982`

Replace:
```python
if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers, TdaEigenSolver):
```
With:
```python
if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
```

## Snippet 11: Fix Root Leakage in Symmetry Dictionary Choice
Type: direct replacement

Where:
- `imdatabasepointcollecter.py::point_correlation_check`

Replace:
```python
if self.current_state > 0 or root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
```
With:
```python
if self.current_state > 0 or (
    self.current_state == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip
):
```

## Snippet 12: Unify Rotor Cluster Coupling Threshold
Type: helper addition

Where:
- module-level constant in both generator and collector (or shared module).

Add:
```python
ROTOR_CLUSTER_COUPLING_THRESHOLD = 0.02  # choose validated value
```
Replace both call sites with:
```python
clusters = build_rotor_clusters(self.symmetry_rotors, coupling_map, threshold=ROTOR_CLUSTER_COUPLING_THRESHOLD)
```
Follow-up edits:
- Validate scientific impact on a representative benchmark set before final value freeze.

---

# 7. Function-by-Function Review

## `IMForceFieldGenerator.__init__` (`imforcefieldgenerator.py:216`)
- Purpose:
  - Prepare defaults, driver bundles, and feature flags.
- Inputs/Outputs:
  - Inputs: ground/excited drivers, roots list.
  - Output: configured generator object.
- Side effects:
  - Initializes many mutable global-run attributes.
- State dependencies:
  - Assumes root policy and driver types.
- Scientific meaning:
  - Sets interpolation and trust-radius defaults that directly affect PES approximation.
- Weaknesses:
  - Wrong ground-state assertion logic.
  - Mutable default arg (`roots_to_follow=[0]`).
  - tuple bug in `use_minimized_structures`.
- Documentation improvements:
  - Document supported driver combinations and current multi-state limitations.

## `IMForceFieldGenerator.set_up_the_system` (`imforcefieldgenerator.py:774`)
- Purpose:
  - Build internal-coordinate/symmetry/rotor metadata for construction.
- Inputs/Outputs:
  - Input: `molecule`, optional existing Z-matrix extraction flags.
  - Output via state mutation: `roots_z_matrix`, `symmetry_information`, `symmetry_rotors`, `conformal_structures`.
- Side effects:
  - Heavy mutation of multiple shared attributes.
- Dependencies:
  - MM topology generation and symmetry detection heuristics.
- Scientific meaning:
  - Defines coordinate manifold and symmetry handling used by interpolation and datapoint generation.
- Weaknesses:
  - Complex nested local functions.
  - Uses loop variable `root` outside loop when building `self.symmetry_rotors`.
- Documentation improvements:
  - Replace index-based `symmetry_information` with named schema.

## `IMForceFieldGenerator.compute` (`imforcefieldgenerator.py:1316`)
- Purpose:
  - Top-level database construction workflow.
- Inputs/Outputs:
  - Input: molecule, optional basis map.
  - Output: `self.im_results` dictionary.
- Side effects:
  - Writes databases, trajectories, summary files; spawns collector runs.
- Dependencies:
  - Every other major class.
- Scientific meaning:
  - Orchestrates dataset generation that defines interpolation model fidelity.
- Weaknesses:
  - Monolithic.
  - ignores `states_basis` input.
  - desired-density gating bug.
  - multi-root branch incomplete.
- Documentation improvements:
  - Explicitly document phase sequence and expected artifacts.

## `IMForceFieldGenerator.determine_molecules_along_dihedral_scan` (`imforcefieldgenerator.py:1835`)
- Purpose:
  - Build scan structures and density bins.
- Inputs/Outputs:
  - Inputs: initial molecules and dihedral spec.
  - Outputs: `point_densities`, `sampled_molecules`, `allowed_deviation`.
- Side effects:
  - none beyond return values.
- Scientific meaning:
  - Defines sampling coverage in torsional space.
- Weaknesses:
  - No schema object for dihedral spec; tuple indices are implicit.
- Documentation improvements:
  - Provide formal input schema (`[(dihedral, n_sampling, state, start), ...]`).

## `IMForceFieldGenerator.determine_datapoint_density` (`imforcefieldgenerator.py:1959`)
- Purpose:
  - Recount datapoints near scan bins from HDF5.
- Side effects:
  - Reads database files.
- Weaknesses:
  - Uses CWD `listdir` existence checks.
- Documentation improvements:
  - Clarify angle-vector metric and label filtering assumptions.

## `IMForceFieldGenerator.add_point` / `add_point_rotor` (`imforcefieldgenerator.py:3501`, `2792`)
- Purpose:
  - Generate and persist new datapoints.
- Inputs/Outputs:
  - Inputs: molecule/basis/root payload + interpolation settings.
  - Output: written HDF5 datapoints and optional cluster registry.
- Side effects:
  - Driver execution + file writes + symmetry metadata mutation.
- Scientific meaning:
  - Adds local Taylor anchors for interpolation.
- Weaknesses:
  - Very large and duplicated with collector versions.
  - hard-coded branches and debug prints.
- Documentation improvements:
  - Separate function docs for core, symmetry, and cluster-bank write paths.

## `IMDatabasePointCollecter.update_settings` (`imdatabasepointcollecter.py:1476`)
- Purpose:
  - Apply dynamics/interpolation runtime settings.
- Weaknesses:
  - Key mismatch fallback (`g_thrsh_...` vs access uses `g_rmsd_...`).
  - `qmc_stop` parsed but unused.
- Documentation improvements:
  - Provide strict settings schema with defaults.

## `IMDatabasePointCollecter.system_from_molecule` (`imdatabasepointcollecter.py:507`)
- Purpose:
  - Build OpenMM system and QM/MM partition.
- Side effects:
  - Writes XML/PDB files; mutates forcefield/system state.
- Weaknesses:
  - Very broad responsibility; includes IO, partitioning, solvent setup.
- Documentation improvements:
  - Split into topology build, system build, and serialization phases.

## `IMDatabasePointCollecter.run_qmmm` (`imdatabasepointcollecter.py:2038`)
- Purpose:
  - Execute QM/MM loop with adaptive interpolation-driven force updates.
- Inputs/Outputs:
  - No explicit scientific return (`None`).
- Side effects:
  - Trajectory/output writes, datapoint growth, trust-radius optimization.
- Weaknesses:
  - Large and complex control flow.
  - sampling bank load bug for multi-root.
- Documentation improvements:
  - Document stop criteria and runtime mode flags.

## `IMDatabasePointCollecter.update_gradient_and_energy` / `update_forces` (`imdatabasepointcollecter.py:3203`, `3370`)
- Purpose:
  - Build QM molecule, evaluate interpolation on all roots, inject forces.
- Scientific meaning:
  - Couples interpolation model to MD propagation.
- Weaknesses:
  - `current_state` hard-set by `_postprocess_root_state_selection`.
  - point-add heuristics and force update tightly coupled.
- Documentation improvements:
  - Clarify units and force conversion equation.

## `IMDatabasePointCollecter.point_correlation_check` (`imdatabasepointcollecter.py:3650`)
- Purpose:
  - Decide if QM validation implies adding datapoint.
- Scientific meaning:
  - Governs adaptive enrichment based on model error indicators.
- Weaknesses:
  - Very complex conditionals and state mutation.
  - root leakage bug, missing argument in bayesian call.
  - `self.skipping_value` reset neutralizes adaptive behavior.
- Documentation improvements:
  - Formalize thresholds and decision tree.

## `IMDatabasePointCollecter.determine_trust_radius_gradient` (`imdatabasepointcollecter.py:5764`)
- Purpose:
  - Optimize per-datapoint confidence radii using reference structures.
- Inputs/Outputs:
  - Returns optimized alpha vector.
- Scientific meaning:
  - Calibrates trust radii to balance interpolation errors.
- Weaknesses:
  - High complexity; implicit config from positional list.
- Documentation improvements:
  - Document objective terms (`energy` vs `force`) and hyperparameters.

## `InterpolationDriver.update_settings` (`interpolationdriver.py:232`)
- Purpose:
  - Parse interpolation config and synchronize datapoint settings.
- Weaknesses:
  - No explicit validation for unsupported interpolation types beyond runtime branch.
- Documentation improvements:
  - Enumerate accepted values and defaults.

## `InterpolationDriver.read_labels` / `read_qm_data_points` (`interpolationdriver.py:528`, `1267`)
- Purpose:
  - Discover datapoint labels and load objects from HDF5.
- Weaknesses:
  - Implicit label parsing by suffix; brittle across naming variants.
- Documentation improvements:
  - Specify label naming invariants.

## `InterpolationDriver.compute` (`interpolationdriver.py:1056`)
- Purpose:
  - Reset current coordinate and run interpolation algorithm.
- Weaknesses:
  - dead `simple` branch call.
- Documentation improvements:
  - Clarify active algorithm choices.

## `InterpolationDriver.shepard_interpolation` (`interpolationdriver.py:1125`)
- Purpose:
  - Compute normalized weighted energy and gradient over datapoints.
- Scientific meaning:
  - Main interpolation equation implementation.
- Weaknesses:
  - Intermediate variables (`hessian_error`) are semantically inconsistent.
- Documentation improvements:
  - Include final equations in docstring.

## `InterpolationDriver.compute_potential` (`interpolationdriver.py:1660`)
- Purpose:
  - Evaluate one datapoint-local Taylor model (with optional symmetry/cluster handling).
- Scientific meaning:
  - Local surrogate model around each reference point.
- Weaknesses:
  - Third return value semantics vary by branch.
- Documentation improvements:
  - Define return tuple contract clearly.

## `InterpolationDriver.cartesian_distance` / `internal_distance` (`interpolationdriver.py:2229`, `2071`)
- Purpose:
  - Build distance denominators and weight derivatives.
- Scientific meaning:
  - Determines interpolation influence and gradients.
- Weaknesses:
  - Complex metric-specific formulas spread across methods.
- Documentation improvements:
  - Add formula references and unit conventions.

## `InterpolationDriver.determine_important_internal_coordinates` (`interpolationdriver.py:2340`)
- Purpose:
  - Identify coordinate constraints tied to interpolation error sources.
- Scientific meaning:
  - Attempts targeted constrained optimization for new datapoints.
- Weaknesses:
  - Large and difficult to verify; currently under-tested branch.
- Documentation improvements:
  - Document ranking criteria and tie-breaking.

## `InterpolationDatapoint.transform_gradient_to_internal_coordinates` / `transform_hessian_to_internal_coordinates` (`interpolationdatapoint.py:448`, `498`)
- Purpose:
  - Convert Cartesian derivatives to internal-coordinate representation.
- Scientific meaning:
  - Essential for Taylor interpolation in internal space.
- Weaknesses:
  - SVD rank handling logic prints warnings but lacks explicit failure policy.
- Documentation improvements:
  - State assumptions on rank and singular value thresholds.

## `InterpolationDatapoint.write_hdf5` / `read_hdf5` (`interpolationdatapoint.py:1132`, `1297`)
- Purpose:
  - Persist/load datapoint tensors + metadata.
- Side effects:
  - Extensive HDF5 dataset writes.
- Weaknesses:
  - Contract depends on label suffix conventions.
- Documentation improvements:
  - Publish dataset schema table and required fields.

---

# 8. Mathematical Framework Documentation

This section distinguishes:
- Confirmed implemented mathematics.
- Inferred intended mathematics.
- Ambiguous/partially implemented mathematics.

## 8.1 Confirmed Implemented Mathematics

### A. Local Taylor Model per Datapoint
From `interpolationdriver.py::compute_potential` and `interpolation_preload_mpi.py::evaluate_candidate_taylor`.

For datapoint \(i\):
\[
U_i(q) = E_i + \Delta q_i^\top g_i + \frac{1}{2}\Delta q_i^\top H_i\Delta q_i
\]
where:
- \(E_i\): reference datapoint energy,
- \(g_i\): internal gradient,
- \(H_i\): internal Hessian,
- \(\Delta q_i = q - q_i\), with dihedral handling below.

When `use_cosine_dihedral=False`, dihedral components use:
\[
\Delta q_{\mathrm{eff}} = \sin(\Delta q_{\mathrm{raw}}), \quad
\text{chain factor } c = \cos(\Delta q_{\mathrm{raw}})
\]
and chain-rule correction is applied before Cartesian back-projection.

### B. Interpolated Energy and Gradient Assembly
From `interpolationdriver.py::shepard_interpolation`.

Given raw weights \(w_i\), normalized weights:
\[
W_i = \frac{w_i}{\sum_j w_j}
\]
Interpolated energy:
\[
E(q) = \sum_i W_i U_i
\]
Interpolated gradient:
\[
\nabla E(q) = \sum_i W_i \nabla U_i + \sum_i U_i \nabla W_i
\]
which is exactly implemented through tensor contractions.

### C. Trust-Radius / Shepard Denominator
From `interpolationdriver.py::shepard_weight_gradient` and related functions.

Using distance \(d\), confidence radius \(R\), exponents \(p,q\):
\[
D = d^2,\quad u = \frac{D}{R^2},\quad \text{denom} = u^p + u^q
\]
Raw weight is \(w = 1/\text{denom}\) (with clipping for numerical stability).

### D. Target-Customized Weight Modifier
From `interpolationdriver.py::VL_target_customized_shepard_weight_gradient` and `internal_distance`.

A second measure \(D_{\mathrm{imp}}\) is used via:
\[
\text{denom}_{\mathrm{tc}} = (u^p+u^q)\exp(\alpha D_{\mathrm{imp}})
\]
thus:
\[
w_{\mathrm{tc}} = \frac{\exp(-\alpha D_{\mathrm{imp}})}{u^p+u^q}
\]
with explicit derivative terms in code.

### E. Symmetry Signature Distance Metric
From `interpolationdriver.py` and `interpolation_preload_mpi.py` grouped signature routines.

For periodic phase differences \(\Delta\theta\):
\[
d^2 = \operatorname{mean}\left(2(1-\cos \Delta\theta)\right)
\]
Gradient term uses:
\[
\nabla d^2 \sim \operatorname{mean}\left(2\sin(\Delta\theta)\,B_{\theta}\right)
\]
where \(B_{\theta}\) are selected Wilson-B rows.

### F. Rotor Coupling Score for Cluster Construction
From `rotorclass.py::rotor_coupling_score`:
\[
C_{ab} = \frac{\|H_{ab}\|_F}{\sqrt{\|H_{aa}\|_F\|H_{bb}\|_F + \varepsilon}}
\]
Clusters are connected components with thresholded \(C_{ab}\).

### G. Internal Coordinate Transform Extensions
From `interpolationdatapoint.py`.

1. Inverse bond mode:
\[
q=1/r,\quad \frac{dq}{dr}=-1/r^2
\]

2. Switched equilibrium-bond transform (`_switched_bond_transform`):
- blends local displacement \(L=r-r_{eq}\) and reciprocal-like term \(R=-r_{eq}^2(1/r-1/r_{eq})\),
- blending weight uses smoothstep polynomial in \(\log(r/r_{eq})^2\) between inner/outer windows.

### H. Trust-Radius Optimization Objective (AlphaOptimizer)
From `alphaoptimizer.py`.

For each reference structure:
- Energy residual term: \(L_E = (E_{interp}-E_{QM})^2\).
- Force residual term: anisotropic blend of parallel/perpendicular components against reference gradient direction.
- Combined loss per structure:
\[
L = e_x L_E + \frac{1}{2}(1-e_x)L_F
\]
Averaged across structures and optimized over confidence radii vector \(\alpha\).

## 8.2 Inferred Intended Mathematics

### A. Bayesian Trust Radius
From `determine_beysian_trust_radius`:
\[
S_i = \sum_k \frac{(\Delta E_{ik})^2}{\sigma^2 d_{ik}^6},\quad
R_i = S_i^{-1/6}
\]
with \(\sigma=0.1\) hard-coded in code.

Status: intended but currently call-signature bugs make this path unreliable.

## 8.3 Ambiguous or Partially Implemented Mathematics

1. `compute_potential` third return value (`hessian_error`) is inconsistent by branch.
2. `simple` interpolation mode is referenced but not implemented.
3. Multi-state state-switching logic is scaffolded but effectively disabled (`current_state` fixed to first root).

---

# 9. Inconsistencies and Unnecessary Code

## Item 1: Dead/Incorrect Always-True Branch
- Location:
  - `imforcefieldgenerator.py:3792`
  - `imdatabasepointcollecter.py:4442`
- Quote:
```python
if 1 == 1 or len(element) == 2 and self.use_minimized_structures[0]:
```
- Why inconsistent:
  - `1 == 1` makes condition always true.
- Recommendation:
  - Delete `1 == 1` and keep intended condition only.

## Item 2: Unused `point_adding_molecule`
- Location:
  - declared: `imdatabasepointcollecter.py:306`
  - consumed externally: `imforcefieldgenerator.py:1812`
- Why inconsistent:
  - Never populated in collector, but collected by generator.
- Recommendation:
  - Either implement population or remove from outputs.

## Item 3: Parsed but Unused `qmc_stop`
- Location: `imdatabasepointcollecter.py:1577-1578`
- Why inconsistent:
  - No subsequent use in runtime stop logic.
- Recommendation:
  - Implement stop check or remove setting.

## Item 4: Unused Reference/Verification Helpers
- Location:
  - `_collect_step_observables_reference`: `imdatabasepointcollecter.py:1764`
  - `_update_qmmm_verification_stats`: `imdatabasepointcollecter.py:1873`
- Why inconsistent:
  - Not called in active run path.
- Recommendation:
  - Remove or integrate behind explicit verification mode.

## Item 5: Unused `determine_conformal_structures`
- Location: `imforcefieldgenerator.py:1911`
- Why inconsistent:
  - No call sites found in workflow.
- Recommendation:
  - Remove or route usage through `set_up_the_system`.

## Item 6: Broken `simple_interpolation` Branch
- Location: `interpolationdriver.py:1099-1102`
- Quote:
```python
if self.interpolation_type == 'simple':
    self.simple_interpolation()
```
- Why inconsistent:
  - No method implementation exists.
- Recommendation:
  - Remove branch or implement method.

## Item 7: Duplicate Initialization Fields
- Location: `imdatabasepointcollecter.py::__init__`
- Examples:
  - `self.cluster_run` set at lines `336` and `365`.
  - `self.coordinates_xyz` set at lines `372` and `388`.
- Why inconsistent:
  - Redundant assignment obscures intended defaults.
- Recommendation:
  - Keep single initialization per field.

## Item 8: CWD-Based File Existence Checks
- Location examples:
  - `imforcefieldgenerator.py:2006`
  - `imforcefieldgenerator.py:3081`
  - `imdatabasepointcollecter.py:4424`
- Why inconsistent:
  - Fails for absolute paths/non-CWD execution.
- Recommendation:
  - Use `os.path.exists(path)` consistently.

## Item 9: Adaptive Skip Logic Neutralized
- Location: `imdatabasepointcollecter.py:3901`
- Quote:
```python
self.skipping_value = 0
```
- Why inconsistent:
  - Immediately erases skip value computed above.
- Recommendation:
  - Remove unconditional reset or gate it explicitly.

## Item 10: Inconsistent Sampling Key Names
- Location:
  - fallback key: `imdatabasepointcollecter.py:1597` (`g_thrsh_kcal_ang_per_atom`)
  - consumption key: `imdatabasepointcollecter.py:3602` (`g_rmsd_thrsh_kcal_ang_per_atom`)
- Why inconsistent:
  - Default fallback does not match read key.
- Recommendation:
  - Standardize one key and add compatibility mapping.

## Item 11: Output Stream Rank Logic Overwritten
- Location: `interpolationdriver.py:94-105`
- Quote:
```python
if ostream is None:
    ...
ostream = OutputStream(sys.stdout)
self.ostream = ostream
```
- Why inconsistent:
  - Reassigns stream for all ranks, ignoring root/non-root setup.
- Recommendation:
  - Remove overriding assignment at line 104.

## Item 12: Misleading Return Variable `hessian_error`
- Location: `interpolationdriver.py::compute_potential`
- Why inconsistent:
  - Returns different semantics by branch (`0.0`, list, or gradient-like vector).
- Recommendation:
  - Rename and enforce single type contract.

## Item 13: Legacy Type/Mode Drift in `use_opt_confidence_radius`
- Location:
  - generator default list with typo: `imforcefieldgenerator.py:392`
  - collector default bool: `imdatabasepointcollecter.py:381`
- Why inconsistent:
  - Different type expectations in different layers.
- Recommendation:
  - Replace with validated structured config.

## Item 14: `run_qmmm` Stats Channel Is Dead
- Location:
  - caller expects stats: `imforcefieldgenerator.py:1788`
  - callee returns `None`: `imdatabasepointcollecter.py:2435`
- Why inconsistent:
  - Hook payload includes `qmmm_stats` but always null.
- Recommendation:
  - Return real stats dictionary.

---

# 10. Developer-Facing Usage Documentation

## Entry Point and Object Creation

Typical workflow entry point:
1. Construct QM drivers (`ScfRestrictedDriver`, `XtbDriver`, optionally excited-state drivers).
2. Instantiate `IMForceFieldGenerator`.
3. Configure key fields.
4. Call `compute(molecule)`.

Example pattern (simplified):
```python
ffg = IMForceFieldGenerator(ground_state_driver=gs_driver, roots_to_follow=[0])
ffg.imforcefieldfiles = {0: "im_database_0.h5"}
ffg.use_symmetry = True
ffg.cluster_run = False
ffg.use_opt_confidence_radius = [False, "multi_grad", 0.5, 0.5]
results = ffg.compute(molecule)
```

## Configuration Surfaces

### Generator-Level (`IMForceFieldGenerator`)
Most important:
- `imforcefieldfiles`: per-root HDF5 output.
- interpolation controls: `interpolation_type`, `weightfunction_type`, `exponent_p`, `exponent_q`, `confidence_radius`, transform flags.
- construction controls: `desired_point_density`, `energy_threshold`, `gradient_rmsd_thrsh`, `start_collect`, `converged_cycle`.
- symmetry/cluster controls: `use_symmetry`, `cluster_run`, `exclude_non_core`.
- datapoint optimization controls: `use_minimized_structures`, `identfy_relevant_int_coordinates`, `use_opt_confidence_radius`.

### Collector-Level (`IMDatabasePointCollecter`)
Configured indirectly via `dynamics_settings` and `interpolation_settings`.
Critical runtime controls:
- MD settings: `ensemble`, `temperature`, `timestep`, `nsteps`, `snapshots`.
- adaptive growth: thresholds and density targets.
- MPI mode flags: `mpi_control_plane_enabled`, `mpi_root_worker_mode`, `mpi_reload_from_hdf5`.
- sampling mode: `sampling_drivers`, `sampling_settings`.

### Interpolation-Level (`InterpolationDriver`)
Important settings:
- `interpolation_type` (currently effectively `shepard` only),
- `weightfunction_type` (`cartesian` or `internal`),
- `exponent_p/q`,
- `use_tc_weights`,
- coordinate transform flags (`use_inverse_bond_length`, `use_eq_bond_length`, `use_cosine_dihedral`).

## Typical Construction Run Lifecycle
1. `compute` prepares system and database paths.
2. Initial datapoints are seeded into HDF5.
3. Scan structures are generated (dihedral bins or default single structure).
4. For each structure, collector runs adaptive QM/MM.
5. New datapoints are appended to HDF5.
6. Optional trust-radius optimization updates confidence radii.
7. Summary outputs are written.

## Outputs Generated
- Primary interpolation DB: `im_database_<root>.h5`.
- Optional sampling DB: `im_database_sampling_<root>.h5`.
- Trajectory PDBs: `trajectory_<state>_<dihedral>_<idx>.pdb`.
- Collector summary HDF5: typically `summary_output.h5`.
- Optional reference structure HDF5: `QM_ref_along_traj.h5`.

## Parameters That Are Easy to Misuse
1. `use_opt_confidence_radius` (type/shape assumptions are fragile).
2. Multi-root settings (`roots_to_follow`) due partial support and hardcoded state selection.
3. Sampling settings key names (`g_rmsd...` vs `g_thrsh...`).
4. File paths if relying on CWD-sensitive checks.

## Extension Points for Developers
1. New interpolation metrics:
   - add in `InterpolationDriver.cartesian_distance/internal_distance` family.
2. New trust-radius objectives:
   - extend `AlphaOptimizer` and `determine_trust_radius_gradient`.
3. New symmetry policies:
   - refactor `set_up_the_system` symmetry builders and mask construction.
4. New state-selection rules:
   - implement policy in `_postprocess_root_state_selection`.

---

# 11. Suggested Merge Preparation Checklist

## Required Pre-Merge Fixes
- [ ] Fix ground-state assertion logic (`all(root == 0 ...)`).
- [ ] Fix desired-density gating bug.
- [ ] Fix Bayesian trust-radius call/signature mismatches.
- [ ] Fix sampling rotor-bank load for all roots.
- [ ] Fix `states_basis` argument handling.
- [ ] Remove tuple bug in `use_minimized_structures` default.
- [ ] Guard/remove broken `simple_interpolation` branch.
- [ ] Fix weighted-centering variable bug (`weights` vs `w`).
- [ ] Replace CWD-based existence checks with `os.path.exists`.
- [ ] Fix `isinstance(drivers, TdaEigenSolver)` typo.
- [ ] Resolve multi-root dead-end behavior in `compute` (hard-fail or full implementation).

## Strongly Recommended Refactors (Near-Term)
- [ ] Extract shared datapoint/cluster registry helpers from generator + collector.
- [ ] Unify rotor cluster threshold into one config/constant.
- [ ] Convert `symmetry_information` to named structured object.
- [ ] Convert `use_opt_confidence_radius` to named config schema.
- [ ] Split mega-functions into stage-specific helpers.
- [ ] Remove dead/legacy helper methods and stale branches.

## Documentation and Scientific Transparency
- [ ] Add explicit equations for interpolation and trust-radius logic to docstrings.
- [ ] Document exact HDF5 schema and label conventions.
- [ ] Document supported/unsupported multi-state behavior clearly.

## Testing
- [ ] Add targeted regression tests for each pre-merge bug fix.
- [ ] Add tests for Bayesian mode and sampling mode.
- [ ] Add tests for multi-root guard behavior.
- [ ] Validate parity before/after refactor on representative molecules.

## Final Validation Before Merge
- [ ] Re-run integration test (`tests/integration/interpolation_test.py`).
- [ ] Verify no unintended change in scientific outputs for baseline systems.
- [ ] Confirm MPI mode still behaves correctly in root-worker setup.
- [ ] Confirm logs/outputs remain interpretable and reproducible.

