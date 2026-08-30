# Benchmarking notes

Notes on the performance of the sparse matrix and SIMD integrals path. Each
entry records what was measured, on what, and what the numbers do and do not
cover.

## Machine

| | |
|---|---|
| CPU | Apple M4 Max, 14 cores (10 performance) |
| Memory | 36 GB |
| OS | macOS 26.6.2 |
| Build | `make -j 12 release` from `src/` |
| Threads | `OMP_NUM_THREADS` unset, OpenMP default |

## Molecules

Geometries are the usual benchmark set.

| molecule | atoms |
|---|---|
| tagrisso | 70 |
| taxol | 110 |
| crambin | 642 |
| ubiquitin | 1231 |

## Dense reconstruction of a sparse matrix

`CSparseMatrix::to_dense`, through the `SparseMatrix.to_numpy` binding, on a
sparse matrix built with the overlap screener at a threshold of 1.0e-14, then
allocated and zeroed. Timings are the best of the stated number of runs.

| molecule | basis | nao | dense GB | sparse GB | to_numpy s | GB/s | runs |
|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.004 | 0.001 | 0.0004 | 10.4 | 5 |
| tagrisso | def2-tzvp | 1345 | 0.013 | 0.004 | 0.0005 | 29.0 | 5 |
| tagrisso | def2-qzvp | 3099 | 0.072 | 0.019 | 0.002 | 41.7 | 5 |
| taxol | def2-svp | 1099 | 0.009 | 0.002 | 0.0003 | 34.8 | 5 |
| taxol | def2-tzvp | 2185 | 0.036 | 0.009 | 0.001 | 36.2 | 5 |
| taxol | def2-qzvp | 4947 | 0.182 | 0.041 | 0.004 | 42.6 | 5 |
| crambin | def2-svp | 6177 | 0.28 | 0.02 | 0.004 | 73.0 | 5 |
| crambin | def2-tzvp | 12063 | 1.08 | 0.10 | 0.017 | 63.5 | 3 |
| crambin | def2-qzvp | 28167 | 5.91 | 0.41 | 0.258 | 22.9 | 1 |
| ubiquitin | def2-svp | 11577 | 1.00 | 0.05 | 0.010 | 103.4 | 3 |
| ubiquitin | def2-tzvp | 22442 | 3.75 | 0.20 | 0.047 | 79.3 | 3 |
| ubiquitin | def2-qzvp | 53197 | 21.08 | 0.87 | 2.20 | 9.6 | 1 |

The cost tracks the size of the dense matrix, not the sparsity: the sparse
values are 0.05 to 0.87 GB against 1 to 21 GB of output. The rates of 60 to 100
GB/s in the middle of the range are consistent with writing the dense array at
memory bandwidth. The fall to 9.6 GB/s at 21 GB is page faulting and memory
pressure on a 36 GB machine, not the traversal.

Two independent runs agreed to within a few percent everywhere except the two
largest cases, which run once and have no minimum to take.

An attempt to split the cost into the zero fill and the sparse scatter gave
inconsistent results, as the allocation and page fault behaviour of the numpy
allocator and of pybind11 differ too much to compare. No split is reported.

## Dense reconstruction against the CMatrix chain

The same job through `CMatrix.to_numpy`, which calls `full_matrix()` and copies
the result into a numpy array. Both sides start from a zeroed container, so only
the reconstruction is measured.

| molecule | basis | nao | dense GB | CMatrix s | SparseMatrix s | speedup | runs |
|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.004 | 0.0004 | 0.0004 | 1.0x | 5 |
| tagrisso | def2-tzvp | 1345 | 0.013 | 0.001 | 0.0005 | 2.1x | 5 |
| tagrisso | def2-qzvp | 3099 | 0.072 | 0.006 | 0.002 | 3.5x | 5 |
| taxol | def2-svp | 1099 | 0.009 | 0.001 | 0.0003 | 2.6x | 5 |
| taxol | def2-tzvp | 2185 | 0.036 | 0.003 | 0.001 | 3.0x | 5 |
| taxol | def2-qzvp | 4947 | 0.182 | 0.015 | 0.004 | 3.5x | 5 |
| crambin | def2-svp | 6177 | 0.28 | 0.024 | 0.004 | 6.3x | 5 |
| crambin | def2-tzvp | 12063 | 1.08 | 0.103 | 0.017 | 6.2x | 3 |
| crambin | def2-qzvp | 28167 | 5.91 | 1.094 | 0.269 | 4.1x | 1 |
| ubiquitin | def2-svp | 11577 | 1.00 | 0.100 | 0.010 | 10.3x | 3 |
| ubiquitin | def2-tzvp | 22442 | 3.75 | 0.425 | 0.047 | 9.0x | 3 |
| ubiquitin | def2-qzvp | 53197 | 21.08 | out of memory | 2.780 | | 1 |

The advantage grows with the size of the molecule and shrinks with the richness
of the basis, as high angular momentum fills the matrix in and leaves less
sparsity to skip. Below a few thousand basis functions the two are
indistinguishable.

Peak resident set size, each path in its own process:

| molecule | basis | dense GB | SparseMatrix GB | CMatrix GB |
|---|---|---|---|---|
| crambin | def2-tzvp | 1.08 | 1.39 | 3.05 |
| ubiquitin | def2-tzvp | 3.75 | 4.18 | 10.07 |
| crambin | def2-qzvp | 5.91 | 6.53 | 15.65 |

The memory is the more consequential difference. The CMatrix chain holds the
blocked matrix, the full submatrix and the numpy copy at once, which measures
2.2 to 2.7 times the dense matrix. The sparse path holds the dense array and the
sparse values, which measures 1.10 to 1.29 times. This is why ubiquitin in
def2-qzvp reconstructs in 2.8 s through the sparse path and does not run at all
through the CMatrix chain on this machine, where it would need about 53 GB.

The default limit of `SparseMatrix.to_numpy` is 8 GB, so crambin in def2-qzvp
passes by default and ubiquitin in def2-qzvp requires an explicit `max_memory`.

## What the timings of the dense reconstruction do not cover

At the time the two sections above were measured the kernel implemented `(s|s)`
alone, so the sparse matrices timed there are built with the real screening
structure of the basis, allocated and zeroed, rather than computed. The traversal
visits every block, combination of basis functions, angular component and atom
pair exactly as it would for a computed matrix, and both sides of the comparison
are handicapped equally, but the values are zeros. The timings therefore stand,
as the work is the same. The note which followed them, that the angular component
handling had never been checked against reference values, does not: every
combination of angular momenta up to six has since been implemented and checked
against the reference driver.

## Profile of the (s|s) overlap kernel

Measured when `(s|s)` was the only kernel, single threaded, on hydrogen and
hydrogen/helium lattices, which were the largest systems it could then reach. The
kernel has since been rewritten, so the shares below describe that kernel and not
the present one, but the conclusion which mattered, that the loop is bound by the
throughput of the exponential and not by the memory it touches, still holds.

### Where the time of the driver goes

Phases of `CSimdOverlapDriver::compute`, with the sparsity construction timed
separately through the exported `SparseMatrix` constructor.

| case | atoms | total s | sparsity s | alloc+zero s | integrals s | integrals % |
|---|---|---|---|---|---|---|
| H512 sto-3g | 512 | 0.0097 | 0.0088 | 0.0000 | 0.0009 | 8.9% |
| H1331 sto-3g | 1331 | 0.0780 | 0.0750 | 0.0000 | 0.0029 | 3.7% |
| H1331 sto-6g | 1331 | 0.0920 | 0.0769 | 0.0000 | 0.0150 | 16.3% |
| H2197 sto-6g | 2197 | 0.2536 | 0.2329 | 0.0001 | 0.0207 | 8.1% |

The integrals are 4 to 16 percent of the call. The rest is building the
sparsity pattern, which a sampling profiler attributes almost entirely to the
stable sort of atom pairs by distance in `CAtomBasisPairGroup::sort_by_distance`:
11380 leaf samples against 1518 for the vector exponential and 81 for the body
of the kernel itself.

### Where the time of the kernel goes

Phases of `compute_ss_overlap`, measured with temporary timers.

| case | calls | mean nvalues | dimensions | alloc+zero | accumulate | copy+fill |
|---|---|---|---|---|---|---|
| H1331 sto-3g | 5 | 496903 | 0.3% | 0.7% | 97.3% | 1.7% |
| H1331 sto-6g | 5 | 758887 | 0.3% | 0.3% | 98.7% | 0.7% |
| H/He729 6-311g | 81 | 35022 | 1.7% | 4.0% | 90.1% | 4.2% |
| H/He125 6-31g | 36 | 2222 | 11.2% | 2.3% | 84.7% | 1.9% |

The accumulation loop is 85 to 99 percent of the kernel. The screening of the
pairs of primitives costs under two percent except on the smallest calls, where
the kernel itself takes a fraction of a millisecond, so precomputing the reach of
a pair of primitives to take the square root and the exponential out of the
bisection would gain nothing.

### What limits the accumulation loop

A standalone benchmark of the loop, built with the flags of the project, in
nanoseconds per atom pair.

| n | array MB | loop as written | same loop without exp | copy only | loop on all threads | threading |
|---|---|---|---|---|---|---|
| 8192 | 0.1 | 3.372 | 0.208 | 0.208 | 6.022 | 0.56x |
| 131072 | 1.0 | 2.623 | 0.279 | 0.279 | 0.603 | 4.35x |
| 500000 | 3.8 | 1.888 | 0.201 | 0.190 | 0.293 | 6.45x |
| 4000000 | 30.5 | 1.403 | 0.147 | 0.146 | 0.206 | 6.81x |

The loop is bound by the throughput of the exponential and not by the memory it
touches. Removing the exponential and keeping the same loads, multiplications
and store makes the loop nine times faster, and the loop sustains 12 to 16 GB/s
where the same loop without the exponential reaches 111 to 152 GB/s. The vector
exponential of the platform is two wide and there is no four wide version, as
linking against `_simd_exp_d4` fails, so the cost of a single evaluation is not
ours to improve. Only evaluating the exponential fewer times would help, and the
number of evaluations is the number of pairs of primitives times the number of
atom pairs they reach, which the screening already minimizes.

### What threading the loop would buy

Measured with a parallel region around the loop, entered above a threshold on
the number of atom pairs, in nanoseconds per atom pair.

| threads | H1331 sto-6g | speedup | H2197 sto-6g | speedup |
|---|---|---|---|---|
| 1 | 1.449 | 1.00x | 1.474 | 1.00x |
| 2 | 0.785 | 1.85x | 0.759 | 1.94x |
| 4 | 0.512 | 2.83x | 0.454 | 3.25x |
| 6 | 0.395 | 3.67x | 0.339 | 4.35x |
| 8 | 0.382 | 3.79x | 0.318 | 4.64x |
| 10 | 0.383 | 3.78x | 0.298 | 4.95x |
| 12 | 0.418 | 3.46x | 0.332 | 4.44x |
| 14 | 0.428 | 3.39x | 0.327 | 4.51x |

Scaling is close to linear to four threads, peaks at ten, which is the number of
performance cores, and degrades beyond it as the remaining threads land on
efficiency cores and the static schedule waits for them. The peak of 4.95x falls
short of the 6.81x of the standalone loop because the pairs of primitives which
reach few atom pairs stay below the threshold and run serially.

The cost of forking and joining the threads is repaid between sixteen thousand
and thirty three thousand atom pairs, measured as 0.86x at 16384 and 1.45x at
32768.

### Why the parallelization was removed again

The kernel and the coordinates carried parallel regions and no longer do. The
runtime of the platform reports `max_active_levels` of one, so a parallel region
inside an active parallel region runs with a single thread. A parallel region in
the kernel and a parallel region over the combinations of basis functions of a
block therefore do not combine, and the outer one silently disables the inner
one. Which of the two levels is right depends on how the work divides between
combinations of basis functions and atom pairs, and with only the `(s|s)` case
implemented every system available is either a single very large call or too few
calls to be representative. The decision was postponed until the kernels of higher
angular momenta existed and the choice could be measured on a real molecule.

That decision has since been taken, and the level is the combinations of basis
functions of a block. See the section below.

## The SIMD overlap driver against the reference overlap driver

Measured once every combination of angular momenta up to six was implemented, so
that the driver runs an arbitrary basis. Each case runs in its own process, as
the largest of them hold twenty one gigabytes and the memory pressure of one case
would otherwise distort the next. Times are the best of three runs, or of a
single run above twelve thousand basis functions, in seconds.

The phases are separated. `sparsity` is the construction of the sparsity pattern,
which the exported constructor of the sparse matrix performs on its own;
`integrals` is the remainder of `compute`, i.e. the coordinates, the solid
harmonics and the kernels; `compute` is their sum, which is what the driver does.
`to dense` is the reconstruction of the dense matrix as a numpy array, and the
two reference columns are the same two steps of the reference driver.

### Single thread

| molecule | basis | nao | dense GB | sparse GB | sparsity | integrals | compute | ref compute | total | integrals | to dense | ref dense |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.00 | 0.001 | 0.0001 | 0.0004 | 0.0005 | 0.0012 | 2.3x | 3.0x | 0.0002 | 0.0003 |
| tagrisso | def2-tzvp | 1345 | 0.01 | 0.004 | 0.0002 | 0.0010 | 0.0012 | 0.0042 | 3.6x | 4.2x | 0.0008 | 0.0010 |
| tagrisso | def2-qzvp | 3099 | 0.07 | 0.017 | 0.0003 | 0.0029 | 0.0033 | 0.0193 | 5.9x | 6.6x | 0.0042 | 0.0062 |
| taxol | def2-svp | 1099 | 0.01 | 0.002 | 0.0003 | 0.0006 | 0.0008 | 0.0027 | 3.4x | 4.9x | 0.0004 | 0.0007 |
| taxol | def2-tzvp | 2185 | 0.04 | 0.009 | 0.0003 | 0.0016 | 0.0019 | 0.0093 | 4.9x | 5.9x | 0.0020 | 0.0032 |
| taxol | def2-qzvp | 4947 | 0.18 | 0.037 | 0.0005 | 0.0050 | 0.0054 | 0.0486 | 8.9x | 9.8x | 0.0102 | 0.0172 |
| crambin | def2-svp | 6177 | 0.28 | 0.025 | 0.0105 | 0.0044 | 0.0149 | 0.0943 | 6.4x | 21.5x | 0.0087 | 0.0268 |
| crambin | def2-tzvp | 12063 | 1.08 | 0.097 | 0.0107 | 0.0170 | 0.0277 | 0.3483 | 12.6x | 20.5x | 0.0849 | 0.1536 |
| crambin | def2-qzvp | 28167 | 5.91 | 0.409 | 0.0112 | 0.0605 | 0.0717 | 1.8501 | 25.8x | 30.6x | 0.4500 | 1.1526 |
| ubiquitin | def2-svp | 11577 | 1.00 | 0.052 | 0.0481 | 0.0081 | 0.0562 | 0.3383 | 6.0x | 41.9x | 0.0220 | 0.1080 |
| ubiquitin | def2-tzvp | 22442 | 3.75 | 0.200 | 0.0505 | 0.0383 | 0.0888 | 1.2537 | 14.1x | 32.7x | 0.2673 | 0.5746 |
| ubiquitin | def2-qzvp | 53197 | 21.09 | 0.872 | 0.0487 | 0.1289 | 0.1777 | too large | | | 3.9617 | |

### All fourteen threads

| molecule | basis | sparsity | integrals | compute | ref compute | total | integrals | to dense | ref dense |
|---|---|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 0.0001 | 0.0004 | 0.0006 | 0.0005 | 0.9x | 1.2x | 0.0002 | 0.0003 |
| tagrisso | def2-tzvp | 0.0002 | 0.0010 | 0.0012 | 0.0014 | 1.1x | 1.3x | 0.0005 | 0.0010 |
| tagrisso | def2-qzvp | 0.0004 | 0.0029 | 0.0033 | 0.0033 | 1.0x | 1.1x | 0.0018 | 0.0062 |
| taxol | def2-svp | 0.0002 | 0.0006 | 0.0008 | 0.0017 | 2.2x | 2.9x | 0.0003 | 0.0007 |
| taxol | def2-tzvp | 0.0003 | 0.0015 | 0.0018 | 0.0024 | 1.3x | 1.6x | 0.0010 | 0.0031 |
| taxol | def2-qzvp | 0.0004 | 0.0050 | 0.0054 | 0.0068 | 1.3x | 1.4x | 0.0043 | 0.0152 |
| crambin | def2-svp | 0.0044 | 0.0023 | 0.0067 | 0.0125 | 1.9x | 5.5x | 0.0039 | 0.0249 |
| crambin | def2-tzvp | 0.0050 | 0.0096 | 0.0145 | 0.0689 | 4.7x | 7.2x | 0.0496 | 0.1351 |
| crambin | def2-qzvp | 0.0051 | 0.0298 | 0.0349 | 0.3223 | 9.2x | 10.8x | 0.2648 | 1.1073 |
| ubiquitin | def2-svp | 0.0191 | 0.0034 | 0.0225 | 0.0445 | 2.0x | 13.2x | 0.0097 | 0.1012 |
| ubiquitin | def2-tzvp | 0.0200 | 0.0170 | 0.0370 | 0.2530 | 6.8x | 14.9x | 0.1633 | 0.6791 |
| ubiquitin | def2-qzvp | 0.0201 | 0.0572 | 0.0773 | too large | | | 1.6281 | |

### What the threads buy, fourteen against one

| molecule | basis | sparsity | integrals | compute | to dense | reference |
|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 0.93x | 0.91x | 0.93x | 0.92x | 2.31x |
| tagrisso | def2-tzvp | 0.90x | 0.94x | 0.94x | 1.56x | 3.03x |
| tagrisso | def2-qzvp | 0.92x | 1.00x | 0.99x | 2.31x | 5.78x |
| taxol | def2-svp | 1.30x | 0.96x | 1.05x | 1.41x | 1.63x |
| taxol | def2-tzvp | 1.14x | 1.02x | 1.04x | 1.93x | 3.90x |
| taxol | def2-qzvp | 1.07x | 0.99x | 1.00x | 2.39x | 7.12x |
| crambin | def2-svp | 2.38x | 1.91x | 2.22x | 2.22x | 7.53x |
| crambin | def2-tzvp | 2.15x | 1.78x | 1.90x | 1.71x | 5.06x |
| crambin | def2-qzvp | 2.20x | 2.03x | 2.06x | 1.70x | 5.74x |
| ubiquitin | def2-svp | 2.52x | 2.39x | 2.50x | 2.25x | 7.61x |
| ubiquitin | def2-tzvp | 2.52x | 2.26x | 2.40x | 1.64x | 4.95x |
| ubiquitin | def2-qzvp | 2.42x | 2.26x | 2.30x | 2.43x | |

The two tables above predate the integrals being formed straight into the values
block and the values blocks no longer being set to zero. Measured again on
fourteen threads with both of those in place, `compute` against the reference:

| molecule | basis | nao | compute | ref compute | speedup |
|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.0005 | 0.0005 | 1.0x |
| tagrisso | def2-tzvp | 1345 | 0.0012 | 0.0031 | 2.6x |
| tagrisso | def2-qzvp | 3099 | 0.0033 | 0.0045 | 1.4x |
| taxol | def2-svp | 1099 | 0.0008 | 0.0011 | 1.4x |
| taxol | def2-tzvp | 2185 | 0.0018 | 0.0016 | 0.9x |
| taxol | def2-qzvp | 4947 | 0.0050 | 0.0065 | 1.3x |
| crambin | def2-svp | 6177 | 0.0064 | 0.0136 | 2.1x |
| crambin | def2-tzvp | 12063 | 0.0090 | 0.0756 | 8.4x |
| crambin | def2-qzvp | 28167 | 0.0216 | 0.3722 | 17.2x |
| ubiquitin | def2-svp | 11577 | 0.0237 | 0.0492 | 2.1x |
| ubiquitin | def2-tzvp | 22442 | 0.0285 | 0.1668 | 5.9x |
| ubiquitin | def2-qzvp | 53197 | 0.0613 | too large | |

The smallest molecules are level with the reference or behind it, as they take
well under a millisecond and their blocks never enter the parallel region.


### What the numbers say

The integrals are twenty to forty two times faster than the reference on a single
thread, which is the comparison of the kernels alone, with the handling of the
matrix excluded from both sides. The advantage is largest where the reference is
weakest, namely many atoms with a small basis, as that is where the screening
discards the most: ubiquitin in def2-svp reaches 41.9x.

The threads are the weak point of the driver, and the cause is not the driver.
The construction of the sparsity pattern, the integrals and the reconstruction of
the dense matrix all saturate at the same 2.2 to 2.5 times, and the last of them
was already parallel before the others were.

**This paragraph originally concluded that the memory bandwidth is the cause.
That conclusion was wrong and is corrected in the section on the phases of
`compute` below.** The integrals scale better than 2.2 times, the figure above
having been measured as `compute` less the sparsity, which silently includes a
serial allocation and zeroing of the values blocks.

The threads therefore let the reference catch up. The advantage of the whole
`compute` falls from between 2.3 and 25.8 times on one thread to between 1.0 and
9.2 times on fourteen, and on the smallest molecules the two are level. Those
cases take under four milliseconds and are below the threshold of the parallel
region of the driver, which is why the driver does not scale on them at all.

The memory is the difference which does not narrow. Ubiquitin in def2-qzvp holds
twenty one gigabytes as a dense matrix and the reference driver cannot compute it
on this machine at all, while the sparse matrix holds 0.872 gigabytes, four
percent of the dense form, and is computed in 0.077 seconds.

The construction of the sparsity pattern is the floor for the large molecules
with a small basis, at 0.048 seconds of the 0.056 seconds of the whole `compute`
for ubiquitin in def2-svp. It threads, as the sort of the atom pairs by distance
was parallel already, but it does not shrink with the basis the way the integrals
do, as it depends on the number of atoms alone.

## The phases of `compute` and what limits their threads

The measurements above take the cost of the integrals as `compute` less the
construction of the sparsity pattern. That is wrong: `compute` also allocates the
values blocks and, at the time, set them to zero, both of which are serial. The
numbers below split `compute` into its four phases directly and supersede the
scaling figures of the previous section.

### Where the time of `compute` goes

Seconds, best of three runs, with the values blocks no longer set to zero.

| molecule | basis | sparsity 1 thr | 14 thr | pair blocks 1 thr | 14 thr | scaling |
|---|---|---|---|---|---|---|
| tagrisso | def2-qzvp | 0.00034 | 0.00040 | 0.00269 | 0.00289 | 0.9x |
| taxol | def2-qzvp | 0.00049 | 0.00043 | 0.00438 | 0.00438 | 1.0x |
| crambin | def2-svp | 0.01085 | 0.00427 | 0.00422 | 0.00198 | 2.1x |
| crambin | def2-tzvp | 0.01121 | 0.00464 | 0.01213 | 0.00444 | 2.7x |
| crambin | def2-qzvp | 0.01186 | 0.00500 | 0.03842 | 0.01159 | 3.3x |
| ubiquitin | def2-svp | 0.05071 | 0.01951 | 0.00850 | 0.00304 | 2.8x |
| ubiquitin | def2-tzvp | 0.04913 | 0.01993 | 0.02424 | 0.00766 | 3.2x |
| ubiquitin | def2-qzvp | 0.05074 | 0.02062 | 0.07841 | 0.02091 | 3.7x |

The allocation of the values blocks is a few tens of microseconds and the
diagonal blocks are tens of nanoseconds, so both are omitted. The two smallest
molecules do not enter the parallel region at all, as their blocks stay below the
threshold of four thousand atom pairs.

The integrals scale between 2.7 and 3.7 times where the molecule is large enough
to enter the parallel region, not the 2.2 times of the previous section. The
construction of the sparsity pattern scales between 2.4 and 2.6 times and is what
the large molecules with a small basis spend their time in: for ubiquitin in
def2-svp it is 0.051 of the 0.059 seconds of the whole `compute`, and the
integrals are almost free beside it.

### What was ruled out, and how

The gap between the threads which are used and the threads which are available
was chased through a series of measurements, each of which excluded a cause.

| cause | measurement | verdict |
|---|---|---|
| memory bandwidth | traffic accounted at 1.86 GB against 0.0298 s, i.e. 62 GB/s of a 296 GB/s ceiling | ruled out, 21 percent used |
| core clocks | pure register arithmetic scales 11.2x on fourteen threads | ruled out |
| page faults | 1 to 256 minor faults per `compute` | ruled out |
| the allocator | 0.0014 s of allocation in a 0.030 s phase | ruled out, five percent |
| the parallel threshold | work in blocks below it is 7.5 percent | ruled out |
| imbalance between combinations | the largest is 1.3 to 3.8 percent of the total | ruled out |
| the working set | batching the atom pairs was slower in both arrangements | ruled out |
| spinning at the barriers | the wait policy changes the wall time by nothing | a symptom, not a cause |

The counters of the CPU showed 2.69e8 cycles per `compute` on one thread against
2.64e9 on ten, i.e. ten times the cycles for twice the speed. That is idle
threads rather than stalled ones, which pointed at the serial phases and led to
the split above.

### What helped and what did not

| change | effect |
|---|---|
| parallel over the combinations of basis functions of a block | the integrals scale 2.7 to 3.7 times |
| the integrals formed straight into the values block | 1.13 to 1.27 times, and it helps a single thread as well |
| the values blocks no longer set to zero | 4 to 10 percent, as every value is written anyway |
| the coordinates and the harmonics made parallel | nothing, measured either way |
| the atom pairs of a block computed in batches | slower, in both arrangements which were tried |

The batches were tried twice. With the batch as the unit of work there are too
few of them, as a block of eighty thousand atom pairs in batches of eight
thousand gives ten tasks for ten threads. With the combinations parallel inside
serial batches the result improves as the batch grows, which is to say as the
batching is switched off. The working set of the harmonics is therefore not what
limits the threads.

## What the kernels are checked against

The integrals have no test in the suite, as no test reaches the driver, so the
evidence that they are right is the comparison below. Every figure is the largest
absolute deviation over the whole matrix unless it says otherwise.

| what | against | over | worst |
|---|---|---|---|
| solid harmonics, l = 1 to 12 | the spherical harmonics of scipy, Racah normalized | 400 random atom pairs, every order and every m | 5.8e-15 relative |
| solid harmonics, l = 1 to 4 | the explicit expressions of Table II of the paper | the same | 1.1e-13 |
| the (s\|l) and (l\|s) kernels, l = 1 to 6 | the analytic form, the overlap of the S functions times the ratio of the exponents raised to l times the harmonic | a custom basis carrying l on one atom and s on the others | 1.7e-16 |
| the (s\|l) and (l\|s) kernels, l = 1 to 6 | the reference overlap driver | the same custom basis | 1.2e-15 |
| every combination up to l = 6 | the reference overlap driver | CO, water and methane in sto-3g, 6-31g, def2-svp, def2-tzvp, def2-qzvp, cc-pvdz, cc-pvtz, cc-pvqz, cc-pv5z and cc-pv6z | 1.8e-15 |
| every combination up to l = 4 | the reference overlap driver | tagrisso, taxol, crambin and ubiquitin in def2-svp, def2-tzvp and def2-qzvp | 1.0e-14 |
| the dense reconstruction | the reference driver through its own dense conversion | the same molecules | 1.0e-14 |

Two of these deserve a note. The reference driver implements angular momenta up
to `I`, i.e. six, so cc-pv6z is the highest basis on which the two can be compared
at all. And the harmonics were taken to twelve because the recursions of the
integrals reach the sum of the two angular momenta, which is twelve when both
sides carry six, so the orders seven to twelve are verified against scipy alone,
no integral of that order existing to compare against.

## The size of the generated code

| | |
|---|---|
| kernels, 49 files | 2609 KB of source, 175 KB of headers |
| driver and dispatcher | 349 KB |
| solid harmonics, 12 orders | 105 KB |
| largest kernel | `SimdOverlapRecIH.cpp`, 267 KB |

The kernels are generated from the recursion descriptors, one file per
combination of angular momenta, in the same shape as the `t2c_overlap` directory
which they parallel and which is 2535 KB across 50 files. The three largest are
`IH`, `HI` and `II`, at 267, 266 and 200 KB, carrying 143, 143 and 91 rows of
angular components.
