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

## What the timings do not cover

The SIMD overlap kernel implements `(s|s)` only, so the sparse matrices timed
above are built with the real screening structure of the basis, allocated and
zeroed, rather than computed. The traversal visits every block, combination of
basis functions, angular component and atom pair exactly as it would for a
computed matrix, and both sides of the comparison are handicapped equally, but
the values are zeros.

The reconstruction was verified against the CMatrix chain to 1.1e-16 on systems
of s functions only, which is as far as the kernel currently reaches. The
angular component and stride handling of the dense mapping runs structurally in
the def2-qzvp cases above, up to g functions, but has never been checked against
reference values for non-zero angular momentum.

## Profile of the (s|s) overlap kernel

Measured with the kernel single threaded, on hydrogen and hydrogen/helium
lattices, which are the largest systems the kernel reaches while only the
`(s|s)` case is implemented.

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
calls to be representative. The decision is postponed until the kernels of higher
angular momenta exist and the choice can be measured on a real molecule.
