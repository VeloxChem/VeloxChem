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
