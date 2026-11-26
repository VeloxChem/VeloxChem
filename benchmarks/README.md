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

| suite | cases | both variants | RI only |
| --- | --- | --- | --- |
| `quick` | caffeine and a nitroxide radical in def2-SVP: SCF (HF, B3LYP, BLYP, CAM-B3LYP, TPSSh, UKS), SCF gradient (RI-JK, RI-J, UKS), TDA, RPA, polarizability | 15 min | 3 min |
| `standard` | diffuse and triple-zeta SCF and gradient, TDA of a 70-atom molecule, unrestricted TDA, TDDFT gradient, CPP, SHG, THG, TPA transition, TPA, 3PA transition, optimization (closed and open shell) | 2 h | 25 min |
| `extended` | Hessian and IR, polarizability gradient, excited-state optimization; noradrenaline in def2-TZVPD (RPA, TPA, gradient); a triplet guanine dimer (UHF, UKS); a TDDFT gradient with RI-J | 7.5 h | 45 min |

Times measured on a 12th Gen Intel Core i9-12900K, one rank, with 8 threads
bound to the 8 performance cores (2026-11-25), as in the binding example under
"Running". `--ri yes` runs the RI rows alone. `python benchmarks/run.py
--list` prints every case.

Whether an RI row fits depends on the memory of the machine. RI-JK for a
60-atom water cluster in def2-SVPD (780 basis functions) was predicted to need
7.4 GB a rank with B3LYP and 14.6 GB with CAM-B3LYP (two ranks, 12.8 GB
available to each), and the CAM-B3LYP run fell back to the four-center
integrals. A row that fell back is reported as a problem, not as a measurement
of RI.

## Cases kept elsewhere

`run.py` also takes a directory of input files from anywhere, for cases that
are not kept in the repository (larger molecules, cluster layouts). A geometry
named by `xyzfile:` is looked up beside the case, then in a `geometries/`
folder beside it or beside one of its two parent folders, then in
`benchmarks/geometries/`. Any other value in the input that names a file
beside the case (point charges, parameter files) is copied with it.

```sh
python benchmarks/run.py path/to/suite --ri yes --launcher srun
python benchmarks/run.py path/to/suite --ri yes --memory-limit 8G
```

`--memory-limit` runs each process under a memory limit without swap, through
a transient systemd scope (Linux with systemd); the program reads the limit as
it reads one set by a batch system. Under a limit a fallback to the
four-center integrals is what the limit is meant to bring about, so it is
reported as a note, and a killed process is the problem.

## Running

```sh
export OMP_NUM_THREADS=8
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

To bind the threads of each case, pass `OMP_PLACES` and `OMP_PROC_BIND`
through the launcher, not the environment of the runner:

```sh
export OMP_NUM_THREADS=8
python benchmarks/run.py quick --launcher 'env OMP_PLACES=cores OMP_PROC_BIND=close'
```

Where numpy uses a BLAS library that runs on OpenMP (the OpenMP build of
OpenBLAS in a conda installation, for example), these variables in the
environment of the runner bind the runner itself to one place when it imports
numpy, and every case it starts inherits that place. The runner stops before
the first case when it may run on fewer logical cores than `OMP_NUM_THREADS`.

With MPI, `OMP_NUM_THREADS` is the number of threads of each rank, and the MPI
launcher has to bind each rank to its own `OMP_NUM_THREADS` cores; without
that, every rank places its threads from the first core, and the ranks share
those cores. The runner does not see the ranks and cannot check this. With
MPICH (Hydra) and 4 threads a rank:

```sh
export OMP_NUM_THREADS=4
python benchmarks/run.py quick --launcher 'env OMP_PLACES=cores OMP_PROC_BIND=close mpirun -n 2 -bind-to core:4'
```

With Open MPI the binding is `--map-by slot:PE=4 --bind-to core`.

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
  linear response functions, C6, harmonic frequencies and IR intensities,
  optimization energies and steps) and from the output where it does not (the
  last molecular gradient, the 3PA strengths, the SHG beta, the THG gamma, the
  polarizability gradient). Transition moments are not recorded: their sign
  follows the phase of the response vectors;
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
belongs with that question, not here. A case is kept in the repository only
when its molecule can be published; others are run from a directory elsewhere
(see "Cases kept elsewhere").

## Geometries

| file | molecule | atoms |
| --- | --- | --- |
| `caffeine.xyz` | caffeine | 24 |
| `nitroxide.xyz` | a nitronyl nitroxide radical, C7H13N2O2 (doublet) | 24 |
| `noradrenaline.xyz` | noradrenaline, C8H11NO3 | 23 |
| `guanine-2.xyz` | a guanine dimer, C12H14N10O2 | 38 |
| `tagrisso.xyz` | osimertinib | 70 |

Their sources are not recorded.

## History

The scripts, data and notes that measured the SIMD integral and RI drivers
during their development (`BENCHMARKS.md`, `benchmarks/scripts/`,
`benchmarks/data/`) were removed after commit `150f4f5a26`. They describe the
code as it was in September 2026, compare paths that no longer exist
(conventional and direct RI-JK), and set driver attributes that were renamed
or removed. Read them with `git show 150f4f5a26:BENCHMARKS.md`.
