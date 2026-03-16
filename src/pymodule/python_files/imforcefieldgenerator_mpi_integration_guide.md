# IMForceFieldGenerator MPI Integration Guide (Minimal-Helper Style)

This document rewrites the MPI integration strategy so it follows the style of `ConformerGenerator`: explicit MPI operations in-place (barrier, bcast, gather), minimal abstraction, and clear phase boundaries.

Goal:

1. Keep `IMForceFieldGenerator` readable for learning.
2. Use the same MPI rhythm as `ConformerGenerator`.
3. Avoid desync while allowing MPI-enabled subcomponents to use workers correctly.

This is a practical implementation guide, not only a conceptual overview.

## 1. Key Idea: You Do Not Need Many Helpers

The previous guide used multiple helper functions to reduce repeated code and enforce correctness. That is valid, but not required.

You can implement a robust MPI integration with:

1. Inline MPI code in each phase.
2. One optional tiny helper for synchronized error handling.
3. Strict discipline: all ranks enter collectives in the same order.

That is exactly how `ConformerGenerator` is written.

## 2. What `ConformerGenerator` Does (And Why It Works)

The `ConformerGenerator` class is a good template because it follows a consistent sequence.

### 2.1 It stores communicator/rank once

Reference:
- `conformergenerator.py:62-75`

Pattern:

```python
if comm is None:
    comm = MPI.COMM_WORLD
self._comm = comm
self._rank = comm.Get_rank()
self._size = comm.Get_size()
```

Why it helps:
- Every MPI call in the class uses the same communicator.
- No accidental mixing of communicators.

### 2.2 It uses barrier after shared file generation

Reference:
- `conformergenerator.py:208-213`

Pattern:

```python
mmff_gen.write_openmm_files(filename=top_file_name)
self._comm.barrier()
```

Why it helps:
- Prevents one rank from reading files while another rank is still writing.

### 2.3 Root builds data, then broadcasts

Reference:
- `conformergenerator.py:523-524`

Pattern:

```python
conformation_dih_arr = comm.bcast(conformation_dih_arr, root=mpi_master())
```

Why it helps:
- One authoritative source creates global work description.
- All ranks receive identical work metadata.

### 2.4 It gathers rank-local performance/results to root

Reference:
- `conformergenerator.py:569-570`, `595-596`, `610-611`

Pattern:

```python
dt_list = comm.gather(conf_dt, root=mpi_master())
gathered_energy_coords = comm.gather(sorted_energy_coords, root=mpi_master())
```

Why it helps:
- Root performs final global decisions, workers only contribute local pieces.

### 2.5 Root-only final postprocessing

Reference:
- `conformergenerator.py:613-736`

Pattern:

```python
if rank == mpi_master():
    # sort global list, filter duplicates, return dict
    return conformers_dict
else:
    return None
```

Why it helps:
- Final ownership is clear.
- No duplicate global side effects.

## 3. Why `IMForceFieldGenerator` Is Harder

`IMForceFieldGenerator` orchestrates many branches and many classes:

1. Conformer generation.
2. Atom-transfer path guessing.
3. Optimization tasks.
4. Energy/gradient/hessian tasks.
5. Database feeding and dynamics setup.

Some of these should be root-only. Some must be collective because the called class itself uses MPI.

If those are mixed without explicit phases, you get desync.

## 4. Minimal-Helper Rules for `IMForceFieldGenerator`

Follow these four rules everywhere.

1. Before every intrinsic-MPI call: `comm.barrier()`.
2. All ranks must call that intrinsic-MPI function.
3. Immediately after: error synchronization (`allgather`) and `comm.barrier()`.
4. Root does final interpretation, then broadcasts compact payload if others need state.

If you always do this, you can keep code mostly inline.

## 5. Current Anchors in `imforcefieldgenerator.py`

Use these locations for integration:

- `_run_optimization(...)` near line 368
- `_generate_conformers_mpi_synchronized(...)` near line 391
- `set_up_the_system(...)` near line 457
- `compute(...)` near line 743
- `determine_atom_transfer_reaction_path(...)` near line 1487
- `compute_energy(...)` near line 2768
- `compute_gradient(...)` near line 2835
- `compute_hessian(...)` near line 2877

## 6. Side-by-Side Mapping: Conformer vs IMForceField

### 6.1 Shared-file safety

Conformer style:

```python
write_openmm_files(...)
comm.barrier()
```

IMForceField equivalent:

```python
# after any root/global setup that produces shared artifacts
comm.barrier()
```

Use this after:
- topology/forcefield files used by all ranks,
- any global setup step whose outputs are read by next phase.

### 6.2 Root creates, all consume

Conformer style:

```python
if rank == root:
    payload = build_work_array(...)
payload = comm.bcast(payload, root=root)
```

IMForceField equivalent:
- Use this for root-only metadata decisions:
  - selected molecules,
  - selected roots,
  - ts/path routing data,
  - state transition flags.

### 6.3 Gather then root decides

Conformer style:

```python
local_result = ...
all_result = comm.gather(local_result, root=root)
if rank == root:
    finalize(all_result)
```

IMForceField equivalent:
- Use for distributed metrics/logging and distributed candidate results.
- Keep final policy/root updates on rank 0.

## 7. Inline MPI Phase Template (No Big Helper Layer)

Use this exact inline block around intrinsic-MPI calls.

```python
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
root = mpi_master()

comm.barrier()  # phase enter
local_err = None
result = None
try:
    result = intrinsic_mpi_call(...)
except Exception as exc:
    local_err = f"rank {rank}: {exc}"

all_err = comm.allgather(local_err)
all_err = [e for e in all_err if e is not None]
if all_err:
    if rank == root:
        raise RuntimeError("\n".join(all_err))
    raise RuntimeError("MPI phase failed on another rank")

comm.barrier()  # phase exit

# root-only postprocess
payload = None
if rank == root:
    payload = build_payload_from(result)

payload = comm.bcast(payload, root=root)
```

This pattern is the same operationally as in `ConformerGenerator`, just generalized to IMForceField phases.

## 8. Integrating `determine_atom_transfer_reaction_path` (Root-Only)

You asked specifically that this task should be head-rank only.

### 8.1 Why root-only here

`determine_atom_transfer_reaction_path(...)` is orchestration and final routing logic, not a distributed kernel where each rank should independently produce separate global paths.

### 8.2 Inline pattern

```python
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
root = mpi_master()

comm.barrier()
if rank == root:
    path_mols, ts_mol = self.determine_atom_transfer_reaction_path(reactants, products, scf=scf)
    payload = {
        "path_xyz": [m.get_xyz_string() for m in path_mols],
        "ts_xyz": ts_mol.get_xyz_string() if ts_mol is not None else None,
        "charge": self.molecule.get_charge(),
        "multiplicity": self.molecule.get_multiplicity(),
    }
else:
    payload = None

payload = comm.bcast(payload, root=root)
comm.barrier()

path_mols = []
for xyz in payload["path_xyz"]:
    m = Molecule.from_xyz_string(xyz)
    m.set_charge(payload["charge"])
    m.set_multiplicity(payload["multiplicity"])
    path_mols.append(m)
```

Comparison to ConformerGenerator:
- same root-generate + bcast structure as `generate()` does for work arrays.


### 8.3 Loop-Level Walkthrough for Your Exact Conformer Block

This subsection maps directly to your `set_up_the_system` code block where you call:

```python
conformal_structures, dihedral_canditates = self._generate_conformers_mpi_synchronized(molecule)
```

and then loop over `dihedral_canditates` while adding transition structures.

#### Critical correction first

In your snippet, this line is wrong:

```python
if rank == rank:
```

It is always `True`. It must be:

```python
if rank == root:
```

Otherwise all ranks execute root-only logic, including `determine_atom_transfer_reaction_path`, which is exactly what we want to avoid.

#### Why this pattern is needed

- `_generate_conformers_mpi_synchronized(...)` already uses all ranks collectively (Conformer-style MPI phase).
- `determine_atom_transfer_reaction_path(...)` should be root-only in this loop.
- After root builds `conformers_plus_ts`, broadcast one canonical payload to all ranks.

This mirrors `ConformerGenerator` design:
- collective compute phase,
- root postprocess,
- broadcast/gather for shared state.

#### Copy/Paste example for this exact loop context

```python
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
root = mpi_master()

# 1) Enter conformer phase together (collective)
comm.barrier()
conformal_structures, dihedral_canditates = self._generate_conformers_mpi_synchronized(molecule)
comm.barrier()

# 2) Root-only loop: build conformers_plus_ts and call transfer-path
conformers_payload = None
if rank == root:
    conformers_plus_ts = {0: {}}

    if len(dihedral_canditates) > 0:
        global_counter = 0
        for entry_idx, entries in enumerate(dihedral_canditates[:]):
            dihedral = entries[0]
            dih_key = tuple([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])

            conformers_plus_ts[0][dih_key] = []
            for i in range(len(entries[1])):
                global_counter += i
                conf_mol = conformal_structures['molecules'][entry_idx + global_counter]
                conformers_plus_ts[0][dih_key].append((conf_mol, 'normal'))

            # root-only call
            _, ts_molecule = self.determine_atom_transfer_reaction_path(
                [conformers_plus_ts[0][dih_key][-2][0]],
                [conformers_plus_ts[0][dih_key][-1][0]],
                scf=False,
            )
            ts_molecule.set_dihedral((3, 4, 6, 10), np.pi / 2, 'radian')
            conformers_plus_ts[0][dih_key].append((ts_molecule, 'transition'))

    else:
        if len(conformal_structures['molecules']) == 0:
            raise RuntimeError('ConformerGenerator returned no conformers.')
        conformers_plus_ts[0][None] = [(conformal_structures['molecules'][0], 'normal')]

    # 3) Root serializes to broadcast-safe payload
    conformers_payload = {}
    for state, dih_map in conformers_plus_ts.items():
        conformers_payload[state] = {}
        for dih_key, mol_entries in dih_map.items():
            key_payload = '__NONE__' if dih_key is None else tuple(dih_key)
            conformers_payload[state][key_payload] = [
                {
                    'xyz': mol_obj.get_xyz_string(),
                    'tag': tag,
                }
                for mol_obj, tag in mol_entries
            ]

# 4) Broadcast built conformer set to all ranks
conformers_payload = comm.bcast(conformers_payload, root=root)
comm.barrier()

# 5) Reconstruct self.conformal_structures identically on all ranks
rebuilt = {}
for state, dih_map in conformers_payload.items():
    rebuilt[state] = {}
    for key_payload, entries in dih_map.items():
        dih_key = None if key_payload == '__NONE__' else tuple(key_payload)
        rebuilt_entries = []
        for item in entries:
            mol_obj = Molecule.from_xyz_string(item['xyz'])
            mol_obj.set_charge(molecule.get_charge())
            mol_obj.set_multiplicity(molecule.get_multiplicity())
            rebuilt_entries.append((mol_obj, item['tag']))
        rebuilt[state][dih_key] = rebuilt_entries

self.conformal_structures = rebuilt
```

#### Step-by-step interpretation of MPI behavior

1. All ranks enter `_generate_conformers_mpi_synchronized` together.
2. That function runs ConformerGenerator MPI-distributed work safely.
3. After return, all ranks have synchronized conformer base data.
4. Only root enters the expensive transfer-path loop and creates TS additions.
5. Root serializes the result into plain payload (XYZ + tags).
6. Payload is broadcast once.
7. All ranks reconstruct identical `self.conformal_structures` for downstream logic.

#### Why this is better than per-rank transfer-path calls

- Avoids duplicated TS search and possible filesystem races.
- Keeps branch decisions deterministic.
- Preserves worker alignment before next collective phase.
- Same high-level rhythm as `ConformerGenerator` (`collective -> root finalize -> shared data sync`).

#### Optional robustness improvement

Add rank-synchronized error handling around root-only loop by broadcasting an error payload if root fails; then all ranks raise the same `RuntimeError` message. This prevents one rank from hanging while another throws.

## 9. Integrating Optimization and Energy/Gradient/Hessian Calls

Your main concern is correct worker participation for classes that already do MPI internally.

### 9.1 Principle

If the called class internally uses MPI collectives, all ranks must call it in lockstep.

This includes likely branches involving:
- `OptimizationDriver.compute(...)`
- SCF/LR solver calls inside `compute_energy(...)`
- Hessian branches in `compute_hessian(...)`

### 9.2 Inline collective call for optimization

```python
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
root = mpi_master()

comm.barrier()
local_err = None
opt_result = None
try:
    opt_result = self._run_optimization(...)
except Exception as exc:
    local_err = f"rank {rank}: {exc}"

errs = comm.allgather(local_err)
errs = [e for e in errs if e is not None]
if errs:
    if rank == root:
        raise RuntimeError("Optimization phase failed\n" + "\n".join(errs))
    raise RuntimeError("Optimization phase failed on another rank")

comm.barrier()

# root final decision
if rank == root:
    optimized_molecule, opt_info = opt_result
```

### 9.3 Inline collective call for compute_energy

```python
comm.barrier()
local_err = None
energy_result = None
try:
    energy_result = self.compute_energy(qm_driver, molecule, basis)
except Exception as exc:
    local_err = f"rank {rank}: {exc}"

errs = comm.allgather(local_err)
errs = [e for e in errs if e is not None]
if errs:
    if rank == root:
        raise RuntimeError("compute_energy failed\n" + "\n".join(errs))
    raise RuntimeError("compute_energy failed on another rank")
comm.barrier()

if rank == root:
    qm_energy, scf_results, rsp_results = energy_result
    payload = {
        "qm_energy": np.asarray(qm_energy).tolist(),
        "has_scf": scf_results is not None,
        "has_rsp": rsp_results is not None,
    }
else:
    payload = None

payload = comm.bcast(payload, root=root)
```

### 9.4 Print Collective Results and Exit the Whole Program (Debug Mode)

You asked how to print values from:

```python
scf_results_mpi = _collective_call_inline_style(comm, rank, root, 'Energy calculation', self.compute_energy, self.drivers['gs'][0], molecule, current_basis)
```

and then terminate the full MPI job.

#### 9.4.1 Root-only debug print after collective return

```python
scf_results_mpi = _collective_call_inline_style(
    comm, rank, root,
    'Energy calculation',
    self.compute_energy,
    self.drivers['gs'][0], molecule, current_basis,
)

if rank == root:
    qm_energy, scf_results, rsp_results = scf_results_mpi
    qm_energy = np.asarray(qm_energy)
    print(f"[DBG] qm_energy shape={qm_energy.shape} first={qm_energy.ravel()[0]:.12f}", flush=True)
    print(f"[DBG] has_scf={scf_results is not None} has_rsp={rsp_results is not None}", flush=True)
```

Notes:
- Print only on root to avoid massive multi-rank output noise.
- Print compact summaries (shape/first value/flags), not full matrices.

#### 9.4.2 Coordinated clean stop (recommended first)

Use this when you want all ranks to exit together without error stack spam.

```python
stop_after_debug = True  # temporary debug flag

# root decides, all ranks receive same decision
stop_after_debug = comm.bcast(stop_after_debug if rank == root else None, root=root)

if stop_after_debug:
    comm.barrier()
    # if inside compute(), prefer returning a valid result object
    return self.im_results
```

If you are not in a function where `return` is possible, use:

```python
if stop_after_debug:
    comm.barrier()
    raise SystemExit(0)
```

#### 9.4.3 Immediate hard stop (debug emergency)

Use this when ranks are hanging and you need guaranteed full termination.

```python
if rank == root:
    print('[DBG] Hard-stopping MPI job after energy checkpoint', flush=True)
comm.Abort(0)
```

Caution:
- `Abort` kills all ranks immediately and may skip clean file/DB finalization.
- Use only as temporary debugging tool.

#### 9.4.4 Why this matches ConformerGenerator style

This is the same control rhythm used in `ConformerGenerator`:

1. collective compute,
2. root inspects/prints global result summary,
3. synchronized continuation decision for all ranks.

### 9.5 Same pattern for gradient/hessian

Use identical phase structure. Keep final interpretation on root.

Comparison to ConformerGenerator:
- mirrors `gather -> root aggregate` and `root decision` sections.


### 9.6 `add_point` After Hessian: Root-Only Database Writing Pattern

This subsection addresses your exact pain point: after collective energy/gradient/hessian calculation in `add_point`, how to write safely to HDF5 and still keep all ranks consistent.

Relevant code area in your file:
- `imforcefieldgenerator.py`, function `add_point(...)`.
- Compute loop around lines where you call:
  - `_collective_call_inline_style(... self.compute_energy ...)`
  - `_collective_call_inline_style(... self.compute_gradient ...)`
  - `self.compute_hessian(...)`
- HDF5 write at `impes_coordinate.write_hdf5(...)`.

#### 9.6.1 Principle

1. Collective compute phases: all ranks participate.
2. File write phases: root rank only.
3. Synchronization metadata: root broadcasts a small receipt (file, label, state), not full hessian arrays.

This avoids file corruption and keeps worker state deterministic.

#### 9.6.2 Recommended loop choreography

For each `mol_basis` entry:

1. All ranks: collective energy call.
2. All ranks: collective gradient call.
3. All ranks: collective hessian call.
4. Root only:
- build `InterpolationDatapoint`,
- assign labels,
- write HDF5 (`write_hdf5`),
- collect `write_receipts`.
5. Root broadcasts `write_receipts`.
6. Barrier before next molecule iteration.

#### 9.6.3 Copy/paste skeleton for your `add_point` loop

```python
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
root = mpi_master()

# --- collective compute ---
scf_results_mpi = _collective_call_inline_style(
    comm, rank, root,
    'Energy calculation',
    self.compute_energy,
    drivers[0], mol_basis[0], mol_basis[1]
)
energies, scf_results, rsp_results = scf_results_mpi

gradient_results_mpi = _collective_call_inline_style(
    comm, rank, root,
    'Gradient calculation',
    self.compute_gradient,
    drivers[1], mol_basis[0], mol_basis[1], scf_results, rsp_results
)
gradients = gradient_results_mpi

hessian_results_mpi = _collective_call_inline_style(
    comm, rank, root,
    'Hessian calculation',
    self.compute_hessian,
    drivers[2], mol_basis[0], mol_basis[1]
)
hessians = hessian_results_mpi

# --- root-only write phase ---
write_receipts = None
if rank == root:
    write_receipts = []

    for number in range(len(energies)):
        target_root = mol_basis[4][number]
        target_file = interpolation_settings[target_root]['imforcefield_file']

        # build datapoint object (same logic as current code)
        impes_coordinate = InterpolationDatapoint(self.roots_z_matrix[target_root])
        impes_coordinate.update_settings(interpolation_settings[target_root])
        impes_coordinate.energy = energies[number]
        impes_coordinate.gradient = gradients[number]
        impes_coordinate.hessian = hessians[number]
        impes_coordinate.transform_gradient_and_hessian()

        # label policy remains root-owned
        label = determine_next_label_somehow(target_file)  # your existing label logic

        # single-writer operation
        impes_coordinate.write_hdf5(target_file, label)

        write_receipts.append({
            'root': int(target_root),
            'file': str(target_file),
            'label': str(label),
            'symmetry': bool(mol_basis[5]),
        })

# broadcast what was written (small payload)
write_receipts = comm.bcast(write_receipts if rank == root else None, root=root)
comm.barrier()

# Optional: worker-side bookkeeping only, no file writes
# for item in write_receipts:
#     update_local_metadata(item)
```

#### 9.6.4 What to broadcast vs what not to broadcast

Broadcast:
- point labels,
- target root/state,
- target file name,
- symmetry flag,
- any control-flow decision needed by workers.

Do NOT broadcast (usually):
- full hessian tensors,
- full gradient arrays,
- full `InterpolationDatapoint` objects.

Reason: workers do not need full raw data if root already committed it to database.

#### 9.6.5 Important correction in your current code

Inside `add_point`, you currently compute gradient twice:

```python
gradient_results_mpi = _collective_call_inline_style(... self.compute_gradient ...)
gradients = gradient_results_mpi[0]

gradients = self.compute_gradient(...)  # second call
```

The second direct call re-executes gradient outside your collective wrapper and can desync or waste time.
Keep only one path (the collective one).

Similarly, hessian should follow your collective pattern as well if the hessian driver has internal MPI behavior.

#### 9.6.6 Why root-only writing is enough

Yes, your understanding is correct: head rank can write the database, then pass only required metadata to all ranks.

This is exactly the ConformerGenerator philosophy:
1. all ranks compute in distributed phase,
2. root consolidates/finalizes,
3. root shares compact state needed for next phase.

#### 9.6.7 Debug checkpoint pattern for this section

```python
if rank == root:
    print(f"[DBG:add_point] wrote {len(write_receipts)} points", flush=True)
    for item in write_receipts[:3]:
        print(f"  -> {item['file']} :: {item['label']} (state {item['root']})", flush=True)

stop_after_write_debug = False
stop_after_write_debug = comm.bcast(stop_after_write_debug if rank == root else None, root=root)
if stop_after_write_debug:
    comm.barrier()
    raise SystemExit(0)
```

This gives you a controlled way to verify label progression and root-only writing before full dynamics continuation.

## 10. Keep `compute()` Readable Without Many Helpers

You can still keep a clear structure with comments and repeated inline blocks.

Suggested phase comments in `compute()`:

```python
# PHASE 0: setup (root-only metadata + bcast)
# PHASE 1: system setup sync
# PHASE 2: conformer phase (collective call)
# PHASE 3: root-only transfer-path decisions + bcast
# PHASE 4: collective optimization/QM calls
# PHASE 5: root-only data routing / add_point
# PHASE 6: sync before handing off to IMDatabasePointCollecter
```

This gives Conformer-like readability without large helper abstraction.

## 11. One Optional Tiny Helper (If You Want)

If you want minimal abstraction but less duplication, keep only this one helper:

```python
def _collective_call_inline_style(comm, rank, root, phase_name, fn, *args, **kwargs):
    comm.barrier()
    local_err = None
    result = None
    try:
        result = fn(*args, **kwargs)
    except Exception as exc:
        local_err = f"rank {rank}: {exc}"
    errs = [e for e in comm.allgather(local_err) if e is not None]
    if errs:
        if rank == root:
            raise RuntimeError(f"{phase_name} failed\n" + "\n".join(errs))
        raise RuntimeError(f"{phase_name} failed on another rank")
    comm.barrier()
    return result
```

This keeps the class close to Conformer style while avoiding repeated try/gather boilerplate.

## 12. Practical "Do/Don’t" for This Class

Do:
- Keep one MPI communicator (`MPI.COMM_WORLD`) for outer orchestration.
- Enter intrinsic-MPI subcalls collectively.
- Use root-only finalization + broadcast compact payload.
- Use barriers at phase edges.

Don’t:
- Let only root call a method that internally uses MPI collectives.
- Let each rank independently run root-style global decisions.
- Mix unrelated collectives without clear phase boundaries.

## 13. Incremental Integration Plan (Learning-Oriented)

### Step 1
Add explicit `comm/rank/root` declarations and phase comments in `compute()`.

### Step 2
Apply inline root-only + bcast pattern to `determine_atom_transfer_reaction_path` call sites.

### Step 3
Wrap optimization call sites in inline collective blocks.

### Step 4
Wrap `compute_energy`, `compute_gradient`, `compute_hessian` call sites in inline collective blocks.

### Step 5
Before entering `IMDatabasePointCollecter.run_qmmm()`, add a final synchronization barrier.

### Step 6
Enable debug sync printouts temporarily and verify all ranks progress through same phases.

## 14. Validation Checklist

1. Serial run:
- `python3 vlx_input.py`
- Must behave exactly as before.

2. MPI small:
- `mpirun -np 2 python3 vlx_input.py`
- No desync and no hangs.

3. Conformer path:
- Ensure conformer generation phase enters/exits on all ranks.

4. Optimization + Hessian branch:
- Ensure all ranks enter these phases together.

5. End-to-end:
- `mpirun -np N python3 vlx_input.py`
- Confirm run continues beyond point addition and back into dynamics.

## 15. Short Conclusion

Yes, you can keep the orchestration in `ConformerGenerator` style.

You do not need many helper functions if you enforce the same four MPI habits everywhere:

1. barrier at phase entry,
2. collective call on all ranks,
3. synchronized error handling,
4. root finalize + broadcast + barrier.

That gives you both clarity for learning and robust MPI behavior in a complex orchestrator class.
