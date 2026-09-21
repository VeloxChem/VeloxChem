# Benchmarks

The measurements this project makes, the data they produced, and the scripts which
produce it. A table in `BENCHMARKS.md` should be rendered from a file here rather
than typed by hand.

```
benchmarks/
  geometries/   the molecules, read by name
  scripts/      scfbench.py, the suites, render.py
  data/scf/     one file per run, rendered beside itself
```

## The rules this exists to enforce

**Never invent a benchmark.** A suite says which molecules, bases, methods and
functionals are measured. Re-running a benchmark means running an existing suite
again, not choosing a grid. Ask before starting any run: they take the machine for
an hour or two and its memory with it.

**One file per run, never appended to.** A run is one suite, one machine, one
commit. That is what makes every ratio inside a file comparable: the `full` row a
speedup is measured against was taken minutes away from the row it is compared
with, on the same machine, at the same commit.

`BENCHMARKS.md` had three vintages of the same caffeine row in three sections, with
nothing to say which was current, and two of them disagreed by a factor of three.
Comparing against the stale one turned a 1.4x into a 3.9x. Hence the format.

**Speedups are computed when a table is rendered, never stored.** A stored ratio
goes stale the moment one of its two columns is measured again.

## Running a suite

```sh
source ~/Environments/vlxenv/bin/activate
source veloxchem.sh
OMP_NUM_THREADS=14 python benchmarks/scripts/scf_laptop.py --machine m4max
```

`--molecule` picks the geometry, `--bases` a subset of the pairs, `--quick` the two
smallest for checking the wiring. The json is rewritten after every row, so a run
stopped early leaves a file which is still valid.

## Running on a node

`scf_node.py` measures one molecule in one basis by the four methods, because on a
node a grid of eight bases is a day. The ranks and the threads are whatever the
launcher gave and are recorded rather than chosen.

```sh
srun --mpi=pmix -N 2 --ntasks-per-node 8 --cpus-per-task 32 \
     --export=ALL,OMP_NUM_THREADS=32,OMP_PROC_BIND=spread,OMP_PLACES=cores,\
LD_PRELOAD=<...>/openblas/lib/libopenblas.so \
     python benchmarks/scripts/scf_node.py --machine epyc9755 \
            --molecule tagrisso --basis def2-tzvp
```

`--molecule` takes a name in `geometries/` or a path to an xyz, so a molecule on the
cluster needs no copy into the repo. `--basis` picks the fitting set by the pairing
rule unless `--aux` overrides it. B3LYP is the default functional; `--functional HF`
for Hartree-Fock. One rank per NUMA domain is the design point: 8 x 32 on one node,
16 x 32 on two.

**VeloxChem's own output is kept**, one file per calculation under `--outdir`, with
`timing` on. The iteration table and the per-iteration breakdown are in there and
nowhere else. They stay on the machine the job ran on and are not tracked: they are
large and cannot be regenerated.

Three things the script refuses or complains about, because each fails looking like
something else:

| | |
| --- | --- |
| the launcher asked for more tasks than the communicator has | srun without `--mpi=pmix` starts one singleton per task instead of failing, and the job measures a fraction of what it was asked to |
| a BLAS serves fewer threads than the rank was given | it does not refuse them; it warns once per thread and dies later somewhere unrelated. The module's OpenBLAS is `MAX_THREADS=48` |
| threadpoolctl is missing | not fatal, but then no BLAS is recorded **and numpy's pool is not resized**, which on the node is the difference between 105 and 3500 Gflop/s |

`--force` measures anyway. Never cap numpy with `OPENBLAS_NUM_THREADS`: every
OpenBLAS in the process reads it, the driver's included, and setting it above a
library's compiled ceiling is what causes the crash above.

`scripts/test_node.py` covers that branch with a stub, because the machine the code
is written on has no threadpoolctl and the machine it matters on does. Run it
directly or under pytest; the repository's own pytest run does not collect it.

Then render it, which writes the markdown and the pdf beside the data:

```sh
python benchmarks/scripts/render.py benchmarks/data/scf/<run>.json
```

`--functional HF` narrows it, `--out -` prints instead of writing, `--no-pdf` skips
the pdf. The pdfs are not tracked; they are regenerated from the json in a second.

## What an SCF record holds

Per calculation: molecule, atoms, basis, nao, fitting set, naux, occupied orbitals,
method, the ri mode **read back from the driver**, functional, convergence
threshold, energy, iterations, and the time split four ways --

| field | what it is |
| --- | --- |
| `ri_setup` | the metric and the B vectors, formed once before the first build |
| `fock_2e_total` | the two-electron build, summed over the iterations |
| `fock_xc_total` | the quadrature, summed the same way |
| `remainder` | the wall clock less the three above |

The split matters. On tagrisso at def2-svpd, 199 seconds of the conventional
route's 357 are B vectors, and the simd route forms the same thing in 17 -- while
the two differ by only 2.9x on the builds. Before the setup had a column of its own
it sat inside the remainder and neither number could be seen.

The methods are `full` (four centre, no approximation), `ri_jk_conventional`,
`ri_jk_simd` (which asks for `automatic` and records what it got), and
`ri_jk_simd_direct`.

## Provenance

Every file records the commit, whether the tree was dirty, the veloxchem version,
the machine and its cpu, the ranks and threads, and every BLAS the process loaded
with its thread count. A number whose origin is not identifiable is not evidence:
two of this year's wrong conclusions came from comparing runs whose configurations
differed in a way nobody had written down.

## The geometries

`caffeine.xyz` and `tagrisso.xyz` for now. Their provenance is not recorded -- they
came from a working directory, not a citable source. Anything added here should
carry where it came from.
