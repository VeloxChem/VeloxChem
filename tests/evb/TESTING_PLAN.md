# EVB driver testing plan

Status: proposed (2026-07-07). No test code written yet — this file is the plan.

## Goal

The EVB driver has grown many code paths (multiple system configurations,
several FEP recalculation modes, reporters, MPI, post-processing). The current
suite (`tests/test_evb.py`) covers only three slices. This plan sets up robust,
tiered testing across the driver, distinguishing **cheap smoke tests** from
**costlier tests that run (cheap) sampling**. Full production runs and
multi-rank MPI are explicitly out of scope for now (Tier C, deferred).

## Decisions (locked)

- **Layout:** new `tests/evb/` package (not one growing file).
- **CI:** all EVB tests stay excluded from CI (keep `timeconsuming`). The
  cheap/costly split is for *local fast-feedback*, not CI gating.
- **Data:** reference the existing `tests/data/evb_*` files in place; do not
  duplicate them.

## Key lever

The committed JSON force fields (`tests/data/evb_ethanol_ff_data.json`,
`evb_ethene_H2O_ff_data.json`) let everything *below* the FF builder run
**without QM**. Load them and you can build systems, evaluate energies, and run
tiny sampling in seconds. QM is therefore quarantined to a single layer; the
rest is cheap.

## Tiers

- **Tier A — Smoke (no MD, ~seconds).** Construct-and-inspect. Build systems for
  every config branch; create an OpenMM `Context` and assert potential energy is
  finite (not NaN) at lambda = 0 / 0.5 / 1. Pure I/O round-trips. Config
  validation errors. Data processing on canned CSVs.
- **Tier B — Cheap sampling (tiny MD, ~seconds-2min each).** Run a real but
  minimal FEP: 9-atom vacuum reaction, ~3 lambda points, `sample_steps` ~20,
  `write_step` = 1. One run per `recalc_mode` and per reporter option, plus the
  full `build -> run -> compute` integration test.
- **Tier C — Deferred (not now).** Production-length runs; multi-rank MPI
  (`mpirun -n 2` for `offload` / async reporter). Stub/document only.

## Markers

All three registered in `tests/conftest.py`. Applied per-file via module-level
`pytestmark`.

- Every EVB test keeps `@pytest.mark.timeconsuming` -> stays out of CI as today.
- Plus one tier marker: `evb_smoke` (Tier A) or `evb_sampling` (Tier B).
- The QM FF-builder test also gets `evb_qm` so it can be excluded from quick runs.

Local usage:

- `pytest tests/evb -m "evb_smoke"` — the inner-loop tier (seconds).
- `pytest tests/evb -m "evb_sampling"` — the MD tier, before committing.
- `pytest tests/evb` — everything (still auto-skipped in CI via `timeconsuming`).

## Package layout

```
tests/evb/
  __init__.py
  conftest.py            # fixtures + assert_energy_finite helper
  test_evb_ffbuilder.py  # Tier A mapping (no-QM) + the one evb_qm test
  test_evb_sysbuilder.py # Tier A: config-branch matrix + energy-finite oracle
  test_evb_driver_io.py  # Tier A: h5/json round-trips, folders, load_initialisation
  test_evb_dataproc.py   # Tier A: canned-CSV compute + dGevb path variants
  test_evb_sampling.py   # Tier B: FEP per recalc_mode, reporters, restart, integration
```

Data files stay in `tests/data/`; fixtures reference them via
`Path(__file__).parents[1] / "data"`.

## Shared infrastructure (build first)

1. Fixtures in `tests/evb/conftest.py`:
   - `evb_ff_pair` — session-scoped; loads the two committed JSON FFs + their
     xyz molecules. Backbone of every non-QM test.
   - `chdir_tmp` — `monkeypatch.chdir(tmp_path)`. `build_systems` /
     `load_initialisation` write folders into cwd, so driver-level tests must run
     in a temp dir.
   - `tiny_fep_config` — config dict with microscopic step counts for Tier B.
   - `built_systems` — runs the system builder once (vacuum), returns
     `(systems, topology, positions)` for energy-eval tests.
2. `assert_energy_finite(system, positions, lambdas)` — creates a `Context` on
   the Reference platform, sets positions, asserts `getState(getEnergy=True)` is
   finite. The cheap correctness oracle for Tier A.
3. Register markers in `tests/conftest.py`.

## Untested surface (what the matrix must reach)

| Layer | File | Cost | Untested paths |
|---|---|---|---|
| FF builder | `reaffbuilder.py` | QM | list vs single input, precomputed charges/Hessians, reparameterize/optimize_mol/optimize_ff flags, forced breaking_bonds, multiplicity overrides, atom mapping |
| System builder | `reactionsystembuilder.py` | cheap | CNT, graphene, E_field, no_reactant, arbitrary solvents, NPT/barostat, constraints, decompose_nb/decompose_bonded, soft-core variants, frozen_atoms, PDB input, config type-validation errors |
| FEP driver | `evbfepdriver.py` | MD | *entirely untested*: recalc_mode {auto, inline, offload, deferred}, n_replicas, fwd/back sweep, isobaric, constrain_H, minimize_every_lambda, restart/skip_initial_equil, MPI async reporter |
| Reporter | `evbreporter.py` | via FEP | *entirely untested*: energy/force/forcegroup rows, GPU recalculator, MPI client/server |
| Data processing | `evbdataprocessing.py` | cheap | partial: only top-level compute; discretised vs analytical dGevb, fitting, plotting |
| Orchestration | `evbdriver.py` | glue+I/O | build_systems folder/file creation, load_initialisation (+restart, OLD_ renaming), h5 round-trip, options.json update, default_system_configurations for every name |

## Tier A matrix (smoke)

- `default_system_configurations`: parametrize every name (`vacuum`,
  `vacuum_NVE`, `water`, `water_NPT`, `E_field`, `no_reactant`, `ts_guesser`,
  `debug`, a real solvent, an unknown name -> `ValueError`). Assert required keys.
- System builder, parametrized over config branches (vacuum, water NVT, water
  NPT, implicit gbn, E_field, no_reactant, CNT, graphene, decompose_nb,
  decompose_bonded, soft-core on/off, constraints, frozen_atoms) -> build +
  `assert_energy_finite`. Keep the serialized-XML reference compare only for the
  2-3 canonical configs already covered (brittle; don't scale it to every branch).
- Config type-validation: wrong-typed option -> expect `TypeError` / `ValueError`.
- Driver I/O: `_save_dict_as_h5` / `_load_dict_from_h5` round-trip;
  `update_options_json` create + merge; `build_systems` produces the expected
  folder/file set; `load_initialisation` reads it back, incl. `OLD_` renaming and
  `restart=` -> `skip_initial_equil` branch.
- Data processing on canned CSVs (keep existing test); add sub-tests for
  discretised vs analytical dGevb.

## Tier B matrix (cheap sampling)

- `test_fep_smoke[recalc_mode]` — parametrize `inline`, `deferred` (single-rank);
  tiny vacuum FEP; assert `*_Energies.csv` / `*_Data_combined.csv` exist with
  expected row counts and finite energies.
- `test_reporter_options` — one tiny run each with `report_forces`,
  `report_velocities`, `report_forcegroups`; assert extra columns/files appear
  and are well-formed.
- `test_fep_then_dataprocessing` — the real integration test: tiny FEP ->
  `compute_energy_profiles` -> assert a finite barrier / dG. Only test that
  exercises the full `build_ff -> build_systems -> run_FEP -> compute` chain
  (start from committed JSON FFs to skip QM).
- `test_restart` — build, run, reload with `restart=` pdb, run again.
- Cheap variants: `n_replicas=2`, `minimize_every_lambda=True`.

Tier B assertions are **structural + sanity-band** (finite, correct shape/row
counts, energies within a loose band), not golden-value — a Langevin integrator
won't reproduce bit-for-bit across platforms.

## Build order

1. Fixtures + `assert_energy_finite`; re-tier the existing three tests into the
   package (add tier markers).
2. Tier A system-builder branch matrix (biggest coverage-per-effort win).
3. Tier A driver I/O tests.
4. Tier B `inline` / `deferred` FEP smoke + FEP->dataprocessing integration.
5. Tier B reporter-option and restart tests.
6. Leave MPI (`offload` / async) and production as documented Tier-C stubs.

## Open questions

1. **QM scope.** Is the single ethanol->ethene+H2O RESP build enough QM coverage,
   or also add a no-QM FF-builder test that feeds precomputed
   `reactant_partial_charges` / `reactant_hessians` to exercise mapping/combination
   cheaply? Recommendation: add the latter — moves most FF-builder coverage into
   Tier A.
2. **CNT/graphene cost.** If even construction is slow, move those from Tier A to
   Tier B. Decide by measuring during implementation.
3. **Existing tests.** Migrate `tests/test_evb.py` into the package, or leave it
   and only add new files? Recommendation: migrate for a single source of truth.
