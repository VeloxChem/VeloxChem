# IMForceFieldGenerator MPI Rewrite Plan

## Scope

This document covers:

1. Coherence and MPI-risk findings in the current implementation.
2. A concrete rewrite plan to move remaining legacy energy/optimization calls into the existing MPI framework.
3. Root-cause analysis and fix plan for `confirm_database_quality` hanging during energy evaluation.

---

## Critical Findings (Current Code)

### 1) `confirm_database_quality` can deadlock because ranks diverge before collective calls

File: `imforcefieldgenerator.py`

- `all_structures = given_molecular_strucutres[root]` is rank-local (`line 2264`).
- `if len(all_structures) == 0: continue` (`line 2269`) allows some ranks to skip the loop.
- Other ranks enter `_collective_call_inline_style(...)` for energy (`line 2398`) which contains barriers/allgather.

If worker ranks have empty `state_specific_molecules` while root has data, root waits forever in the collective phase.

### 2) Non-deterministic random selection can desynchronize rank loop counts

File: `imforcefieldgenerator.py`

- Rank-local `random.sample(...)` is used in selection logic (`line 2336`) before collectives.

Even with equal list sizes, different random outcomes per rank can lead to different iteration lengths and collective mismatch.

### 3) Debug/process-termination statements still present in production path

File: `imforcefieldgenerator.py`

- `exit()` (`line 2277`)
- `raise SystemExit(0)` (`line 2406`)

These can terminate one rank while others continue, causing MPI hangs or hard aborts.

### 4) Invalid driver access in excited-state branch

File: `imforcefieldgenerator.py`

- `drivers` is a tuple but code uses `drivers['es']` (`line 2458`, `line 2465`).

This is a runtime type error path.

### 5) Variable name mismatch in confirm-db write path

File: `imforcefieldgenerator.py`

- Code uses `reference_energy` (`line 2485`) while computed variable is `reference_energies` (`line 2399`).

This is a runtime NameError when branch is hit.

### 6) Sampling screen currently does not short-circuit QM call path

File: `imdatabasepointcollecter.py`

- `_sampling_screen_check(...)` is called (`line 3796`) and may return `skip_full_qm=True`.
- Full QM energy/gradient still runs afterward (`line 3818`, `line 3824`) due missing `continue/return`.

Result: sampling pre-screen does not actually prevent expensive QM evaluation.

### 7) Sampling settings key mismatch

File: `imdatabasepointcollecter.py`

- Default key in `update_settings`: `'g_thrsh_kcal_ang_per_atom'` (`line 1283`)
- Reader in `_sampling_screen_check`: `'g_rmsd_thrsh_kcal_ang_per_atom'` (`line 3700`)

This can silently misconfigure thresholds depending on which dictionary path is used.

### 8) Sampling labels reload mismatch

File: `imdatabasepointcollecter.py`

- Sampling reload sets `driver_object.labels = self.sorted_state_spec_im_labels[root]` (`line 4943`), but this is the main DB label list, not sampling labels.

Potential downstream mismatch between in-memory sampling data points and labels.

### 9) Metadynamics torsion indexing inconsistent with validator

File: `imdatabasepointcollecter.py`

- Validator enforces 0-based atoms (`line 1129`).
- Torsion force subtracts 1 (`line 926`), while distance/angle do not.

This likely shifts torsion atoms incorrectly for metadynamics CVs.

---

## Root Cause: Why `confirm_database_quality` Hangs

The hang is most likely a **collective desynchronization**:

1. `run_qmmm()` in root-worker mode only accumulates full `state_specific_molecules` on root.
2. Workers return with empty per-root molecule lists.
3. `confirm_database_quality()` runs on all ranks and branches by rank-local list length.
4. Some ranks skip, others call `_collective_call_inline_style` (barriers + allgather).
5. Barrier/allgather never completes.

Secondary contributors:

- rank-local randomness before collectives,
- leftover `exit()`/`SystemExit`,
- potential tuple/dict errors in ES path.

---

## Fix Plan for `confirm_database_quality` (Immediate)

### A. Make structure input rank-consistent

1. Build `all_structures` payload on root only.
2. Broadcast serialized structures (XYZ + charge + multiplicity) to all ranks.
3. Reconstruct molecule objects identically on each rank.

### B. Make random selection deterministic across ranks

1. Perform random selection on root only.
2. Broadcast selected indices (or serialized selected molecules) to all ranks.
3. Ensure each rank runs identical loop counts before collectives.

### C. Remove unsafe termination

1. Replace `exit()` and `raise SystemExit(0)` with proper control flow (`raise RuntimeError` or clean `break/continue`).

### D. Fix variable/type errors in ES branch

1. Replace `drivers['es']` with tuple-index access.
2. Replace `reference_energy` with `reference_energies`.

### E. Add MPI safety assertions

Before each collective phase, assert:

- same number of structures,
- same random-selection length,
- same loop counter expectations across ranks.

---

## MPI Rewrite Plan for `IMForceFieldGenerator`

### Phase 1: Introduce one MPI-safe execution API (single entry point)

Create helper wrappers in `IMForceFieldGenerator`:

- `_energy_mpi_safe(...)`
- `_gradient_mpi_safe(...)`
- `_hessian_mpi_safe(...)`
- `_opt_mpi_safe(...)`

Behavior:

- In pure collective phase: call all-rank collectives.
- In root-worker/root-local phase: temporarily enforce local communicator (`COMM_SELF`) on participating drivers.

This mirrors the robust approach already used in `IMDatabasePointCollecter`.

### Phase 2: Replace legacy direct compute calls

Replace direct calls in bootstrap/setup/confirm/add-point paths:

- `self.compute_energy(...)`
- `self._run_optimization(...)`

with the safe wrappers above. Prioritize these hot spots:

1. `confirm_database_quality`
2. `add_point`
3. setup/initial-geometry optimization branches in `compute`
4. `_bootstrap_sampling_db_from_abinito_db`

### Phase 3: Separate root-only logic from collective logic

Refactor methods into explicit sections:

- `root_local_preparation` (selection, bookkeeping, file decisions)
- `collective_qm_phase` (energy/gradient/hessian collectives)
- `root_only_io_phase` (HDF5 writes, logs, plotting)

Rule: no rank-local branch is allowed to skip a collective call in collective sections.

### Phase 4: Standardize MPI data transport

For every method that may run on all ranks:

1. Serialize molecule payloads explicitly for broadcasts.
2. Reconstruct canonical molecule objects on all ranks.
3. Broadcast deterministic control parameters (selected indices, roots, thresholds).

### Phase 5: Harden HDF5 update model

Use root-only writes for interpolation/sampling files and enforce barriers around reload points.

For sampling bootstrap:

- guard `_bootstrap_sampling_db_from_abinito_db` to root write only,
- then trigger synchronized reload of sampling interpolation drivers.

### Phase 6: Fix sampling/metadynamics coherence issues

1. In `point_correlation_check`, actually short-circuit full QM when sampling screen passes.
2. Unify sampling threshold keys to `g_rmsd_thrsh_kcal_ang_per_atom`.
3. Fix sampling reload labels to use sampling label list.
4. Fix metadynamics torsion indexing to match 0-based validation or change validation if 1-based is intended.

### Phase 7: Verification and regression tests

Minimum test matrix:

1. Single-rank run.
2. Multi-rank collective mode (`mpi_root_worker_mode=False`).
3. Multi-rank root-worker mode (`mpi_root_worker_mode=True`).
4. `confirm_database_quality` with empty/non-empty `state_specific_molecules`.
5. Sampling enabled/disabled.
6. Metadynamics enabled/disabled.

Required assertions:

- no rank deadlock,
- identical number/order of collective calls per rank,
- consistent point counts across ranks after add/reload,
- deterministic behavior under fixed RNG seed.

---

## Practical Implementation Order

1. Remove `SystemExit`/`exit` and fix ES tuple/dict + variable-name bugs.
2. Make `confirm_database_quality` structure and random-selection broadcasts deterministic.
3. Introduce MPI-safe wrappers in `IMForceFieldGenerator`.
4. Migrate direct compute/optimization calls to wrappers.
5. Fix sampling short-circuit and settings-key mismatch.
6. Fix sampling label reload and metadynamics torsion indexing.
7. Run multi-rank validation scenarios.

