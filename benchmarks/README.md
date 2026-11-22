# Benchmarks

Cases that are run again to check that performance and results have not
regressed. Each case is an ordinary VeloxChem input file, run through the
command line (`python -m veloxchem`) twice: with `ri: yes` and with `ri: no`.
The four-center row is the reference for the RI row's results and the
baseline for its speedup, and both are taken in the same run on the same
machine.

```
benchmarks/
  cases/<suite>/*.inp   the cases; the first comment line describes each
  geometries/           the molecules the cases read with xyzfile:
  run.py                runs cases and writes results/<run id>.json
  compare.py            summarizes one run, or compares two
  results/              not tracked
```

## Suites

| suite | cases | four-center wall time, 12 cores (rough) |
| --- | --- | --- |
| `quick` | caffeine and a nitroxide radical in def2-SVP: SCF (HF, B3LYP, BLYP, CAM-B3LYP, UKS), SCF gradient (RI-JK, RI-J, UKS), TDA, RPA | 10 min |
| `standard` | diffuse and triple-zeta SCF and gradient, TDA of a 70-atom molecule, TPA transition, reduced TPA, 3PA transition, optimization (closed and open shell) | 1 to 1.5 h |

The RI rows take a fraction of that. `python benchmarks/run.py --list` prints
every case.

## Running

```sh
export OMP_NUM_THREADS=12
python benchmarks/run.py quick                       # both variants
python benchmarks/run.py quick standard --ri yes     # RI rows only
python benchmarks/run.py benchmarks/cases/quick/scf_hf.inp
python benchmarks/run.py quick --launcher 'mpirun -n 2'
python benchmarks/run.py quick --repeat 3            # smallest wall time counts
```

Run on an otherwise idle machine and set `OMP_NUM_THREADS`; the runner warns
when it is unset and when a BLAS library in the process serves fewer threads
than it asks for. Reinstall after changing the code: the runner warns when the
build commit printed in the output is not the checked-out commit (only for a
build configured with `ENABLE_COMMIT_HASH`).

Each variant runs as its own process in a fresh directory,
`results/<run id>/<suite>/<case>/<variant>/`, so no checkpoint file of an
earlier run is read as a restart, and a killed case costs one row and not the
run. The input and the output are kept there; checkpoint and h5 files are
removed unless `--keep-files` is given.

## Comparing

```sh
python benchmarks/compare.py results/RUN.json              # RI against four-center
python benchmarks/compare.py results/BASE.json results/NEW.json
```

With one run, each case's RI row is set against its four-center row: speedup,
peak memory of each, and the largest difference of each result. The
differences are the error of the fitting and are reported, not judged.

With two runs, each row of the new run is set against the same case and variant
of the base run. Results outside the tolerances in `compare.py` fail; wall time
and peak memory are reported as new/base ratios and marked above 10 % growth
(`--time-tol`, `--mem-tol`). Ratios are only meaningful between runs with the
same machine, threads, ranks and launcher, and the differences are printed
first. `compare.py` exits with 1 when a result fails or a row has a problem.

## What a row records

- the exit code, the wall time of each repetition, and the peak resident set
  size (the largest process on this host: with `mpirun`, the largest rank);
- the results, from the final h5 file where it holds them (SCF energy,
  excitation energies, oscillator and TPA strengths, TPA cross sections,
  optimization energies and steps) and from the output where it does not (the
  last molecular gradient, the 3PA strengths). Transition moments are not
  recorded: their sign follows the phase of the response vectors;
- the solver stages, iterations and times the output prints, and its warnings;
- the RI flags of the SCF (`ri`, `ri_jk`, `ri_coulomb`, `ri_fallback`).

A row has a problem when the process failed, no results were read, or an RI row
did not use RI: the SCF reports a fallback in the h5 file, and the gradient and
response drivers in a warning that names the four-center integrals. Such a row
measured the four-center path and its timings say nothing about RI.

Every run records the commit, whether tracked files were modified, the build
commit, the veloxchem version and path, the machine and CPU, the launcher, the
`OMP_*`, `OPENBLAS_*`, `MKL_*` and `VLX_*` environment, and the BLAS libraries
with their thread counts (when threadpoolctl is installed).

## Rules

- **A run is never appended to.** One file is one set of cases, one machine,
  one commit. The runner refuses a run id that exists.
- **Ratios are computed when compared, never stored.** A stored ratio goes
  stale when either of its rows is measured again.
- **Results are local.** Timings belong to a machine, so runs are not tracked;
  keep a base run of your machine and compare against it.

## Adding a case

Put an input file in a suite directory, with a first comment line that says
what it measures. Use `xyzfile:` with a name from `geometries/` (or a file
beside the case), and do not set `ri`: the runner adds it. A case should stay
meaningful as the code changes. A measurement that answers one design question
belongs with that question, not here.

## Geometries

`caffeine.xyz` (24 atoms), `nitroxide.xyz` (a nitronyl nitroxide radical,
C7H13N2O2, 24 atoms, doublet) and `tagrisso.xyz` (osimertinib, 70 atoms). Their
sources are not recorded.

## History

The scripts, data and notes that measured the SIMD integral and RI drivers
during their development (`BENCHMARKS.md`, `benchmarks/scripts/`,
`benchmarks/data/`) were removed after commit `150f4f5a26`. They describe the
code as it was in September 2026, compare paths that no longer exist
(conventional and direct RI-JK), and set driver attributes that were renamed
or removed. Read them with `git show 150f4f5a26:BENCHMARKS.md`.
