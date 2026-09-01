# Benchmarking notes

Notes on the performance of the sparse matrix and SIMD integrals path. Each
entry records what was measured, on what, and what the numbers do and do not
cover.

The tables which describe the three drivers as they stand — the def2, diffuse and
correlation consistent tables of the overlap and the kinetic energy, the parallel
scaling of the kinetic energy, and the fitting set table of the two-center
Coulomb — were all measured again on one build, so they can be read against each
other. The sections which record an intermediate state of the code, the kernel
profile, the sweeps of the blocks and the block floor, the Instruments findings
and the dense reconstruction, keep the numbers of the run which produced them and
were not repeated; they say so where it matters.

## Machine

| | |
|---|---|
| CPU | Apple M4 Max, 14 cores (10 performance) |
| Memory | 36 GB |
| OS | macOS 26.6.2 |
| Build | `make -j 12 release` from `src/` |
| Threads | stated per table, `OMP_NUM_THREADS` set explicitly where it matters |

## Molecules

Geometries are the usual benchmark set.

| molecule | atoms |
|---|---|
| tagrisso | 70 |
| taxol | 110 |
| crambin | 642 |
| ubiquitin | 1231 |

The bases are def2-svp, def2-tzvp and def2-qzvp throughout, their diffuse
counterparts def2-svpd, def2-tzvpd and def2-qzvpd in the section on them, and the
correlation consistent sets to sextuple zeta in the section at the end.

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

That decision was taken, and the level was the combinations of basis functions
of a block. It has since been removed again, along with every other parallel
region of this path, leaving the vectorization alone. See the section on the
single-threaded driver at the end.

## The SIMD overlap driver against the reference overlap driver

**The tables of this section are superseded by the section on the cold and the
warm cost near the end, for two reasons. They were taken before the atom pairs
were sorted by radix, which more than halves the construction of the sparsity
pattern, and they timed the cases above twelve thousand basis functions once and
cold while timing the smaller ones as the best of three and warm, which mixes two
different measurements in one column.**

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

**The scaling columns of this section describe parallel regions which have since
been removed; see the section on the single-threaded driver at the end. The split
of `compute` into its phases, and the causes ruled out below, still stand.**

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

**Three of the rows below concern parallelizations which no longer exist. They
record what the threads bought while they were there, not what the driver does
now.**

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

## The cost when the memory is cold and when it is warm

**The single-thread table of this section still describes the driver. The
fourteen-thread table and the one which follows it do not, the parallel regions
having since been removed; see the section on the single-threaded driver at the
end.**

Every table above this one timed the larger cases once and the smaller ones as
the best of three. That mixes two measurements in one column, as the first call
faults in the pages of a freshly allocated sparse matrix and the later ones do
not. The tables below give every case the same treatment, one cold call followed
by the best of three warm ones, and separate the two rather than choosing between
them. Which is the honest number depends on the caller: a single overlap matrix
pays the cold cost, a loop which computes one repeatedly pays the warm one.

The difference is not small, and it falls on the driver harder than on the
reference, which allocates less: crambin in def2-qzvp on fourteen threads is
0.0385 seconds cold against 0.0141 warm. It also falls harder as the threads
grow, the faults of the first touch being serialized by the kernel, which is why
the cold column scales worse than the warm one.

### Single thread

| molecule | basis | nao | sparse GB | compute cold | compute warm | ref cold | ref warm | cold | warm |
|---|---|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.001 | 0.0019 | 0.0012 | 0.0014 | 0.0012 | 0.71x | 0.99x |
| tagrisso | def2-tzvp | 1345 | 0.004 | 0.0028 | 0.0018 | 0.0046 | 0.0040 | 1.66x | 2.19x |
| tagrisso | def2-qzvp | 3099 | 0.017 | 0.0051 | 0.0037 | 0.0233 | 0.0204 | 4.60x | 5.53x |
| taxol | def2-svp | 1099 | 0.002 | 0.0024 | 0.0014 | 0.0029 | 0.0027 | 1.24x | 1.96x |
| taxol | def2-tzvp | 2185 | 0.009 | 0.0034 | 0.0024 | 0.0106 | 0.0097 | 3.10x | 3.96x |
| taxol | def2-qzvp | 4947 | 0.037 | 0.0078 | 0.0059 | 0.0588 | 0.0524 | 7.58x | 8.93x |
| crambin | def2-svp | 6177 | 0.025 | 0.0093 | 0.0074 | 0.1050 | 0.0977 | 11.26x | 13.14x |
| crambin | def2-tzvp | 12063 | 0.097 | 0.0210 | 0.0150 | 0.3543 | 0.3249 | 16.87x | 21.60x |
| crambin | def2-qzvp | 28167 | 0.409 | 0.0629 | 0.0438 | 1.8982 | 1.7470 | 30.17x | 39.90x |
| ubiquitin | def2-svp | 11577 | 0.052 | 0.0219 | 0.0181 | 0.3850 | 0.3502 | 17.55x | 19.33x |
| ubiquitin | def2-tzvp | 22442 | 0.200 | 0.0444 | 0.0344 | 1.2773 | 1.1472 | 28.79x | 33.34x |
| ubiquitin | def2-qzvp | 53197 | 0.872 | 0.1210 | 0.0877 | too large | | | |

### All fourteen threads

| molecule | basis | compute cold | compute warm | ref cold | ref warm | cold | warm |
|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 0.0016 | 0.0006 | 0.0016 | 0.0016 | 1.06x | 2.48x |
| tagrisso | def2-tzvp | 0.0027 | 0.0013 | 0.0019 | 0.0017 | 0.68x | 1.33x |
| tagrisso | def2-qzvp | 0.0047 | 0.0032 | 0.0065 | 0.0037 | 1.37x | 1.17x |
| taxol | def2-svp | 0.0019 | 0.0009 | 0.0023 | 0.0009 | 1.19x | 1.09x |
| taxol | def2-tzvp | 0.0029 | 0.0019 | 0.0030 | 0.0023 | 1.01x | 1.25x |
| taxol | def2-qzvp | 0.0069 | 0.0048 | 0.0128 | 0.0064 | 1.85x | 1.34x |
| crambin | def2-svp | 0.0057 | 0.0033 | 0.0206 | 0.0127 | 3.59x | 3.83x |
| crambin | def2-tzvp | 0.0129 | 0.0064 | 0.0750 | 0.0440 | 5.79x | 6.92x |
| crambin | def2-qzvp | 0.0385 | 0.0141 | 0.3411 | 0.2099 | 8.86x | 14.85x |
| ubiquitin | def2-svp | 0.0123 | 0.0080 | 0.0697 | 0.0436 | 5.68x | 5.45x |
| ubiquitin | def2-tzvp | 0.0255 | 0.0120 | 0.2386 | 0.1396 | 9.35x | 11.65x |
| ubiquitin | def2-qzvp | 0.0642 | 0.0251 | too large | | | |

### What the threads buy, fourteen against one

| molecule | basis | compute cold | compute warm | reference warm |
|---|---|---|---|---|
| tagrisso | def2-svp | 1.24x | 1.87x | 0.75x |
| tagrisso | def2-tzvp | 1.03x | 1.44x | 2.37x |
| tagrisso | def2-qzvp | 1.07x | 1.16x | 5.48x |
| taxol | def2-svp | 1.23x | 1.63x | 2.93x |
| taxol | def2-tzvp | 1.16x | 1.32x | 4.16x |
| taxol | def2-qzvp | 1.12x | 1.23x | 8.15x |
| crambin | def2-svp | 1.62x | 2.24x | 7.67x |
| crambin | def2-tzvp | 1.62x | 2.36x | 7.38x |
| crambin | def2-qzvp | 1.63x | 3.10x | 8.32x |
| ubiquitin | def2-svp | 1.79x | 2.26x | 8.03x |
| ubiquitin | def2-tzvp | 1.74x | 2.87x | 8.22x |
| ubiquitin | def2-qzvp | 1.89x | 3.50x | |

### What these numbers say

The driver is strongest on a single thread, where it is between eleven and forty
times the reference on the large cases. The reference threads better than the
driver does, seven to eight times against two to three and a half, so its
disadvantage narrows to between four and fifteen times on fourteen threads. The
small molecules are a wash or a loss, tagrisso in def2-svp being level on one
thread and behind on fourteen, as they take a millisecond or two and never enter
the parallel region of the driver at all.

A measurement taken before this one appeared to show `compute` scaling worse than
either of its phases, which would have meant that something inside it grows with
the threads and belongs to neither. It does not. The phases were warm by the time
they were measured within a run while the total was not, so the two were not
comparable. Measured warm throughout, `compute` scales between 2.2 and 3.5 times,
which is what its phases scale to, and the whole of it is accounted for: the
phases and the time the caller sees agree to thirty microseconds.

## The single-threaded driver

Every parallel region of this path was removed: the coordinates, the solid
harmonics, the combinations of basis functions of a block, the description of the
sparsity patterns, the ordering of the atom pairs and the reconstruction of the
dense matrix. The `omp simd` directives of the kernels and of the harmonics are
untouched, so the vectorization is what it was and only the threads are gone.
Every table above which carries a thread count describes the code before that.

Measured with the methodology of the section above, one cold call followed by the
best of three warm ones, each case in its own process. The molecules and the
dimensions of their bases are the same as everywhere else in this file, so the
cases compare one for one.

### Where the time goes

One thread, seconds. `sparsity` is the construction of the sparsity pattern
through the exported `SparseMatrix` constructor, `compute` is the whole driver
call, and the two reference columns are the same call on the reference driver.

| molecule | basis | nao | sparse GB | sparsity cold | sparsity warm | compute cold | compute warm | ref cold | ref warm | cold | warm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.001 | 0.0009 | 0.0008 | 0.0012 | 0.0012 | 0.0021 | 0.0012 | 1.68x | 1.06x |
| tagrisso | def2-tzvp | 1345 | 0.004 | 0.0010 | 0.0009 | 0.0020 | 0.0018 | 0.0052 | 0.0041 | 2.65x | 2.25x |
| tagrisso | def2-qzvp | 3099 | 0.017 | 0.0011 | 0.0010 | 0.0045 | 0.0036 | 0.0226 | 0.0201 | 5.05x | 5.60x |
| taxol | def2-svp | 1099 | 0.002 | 0.0009 | 0.0008 | 0.0013 | 0.0013 | 0.0038 | 0.0028 | 2.82x | 2.15x |
| taxol | def2-tzvp | 2185 | 0.009 | 0.0009 | 0.0009 | 0.0026 | 0.0023 | 0.0113 | 0.0094 | 4.27x | 4.06x |
| taxol | def2-qzvp | 4947 | 0.037 | 0.0011 | 0.0010 | 0.0073 | 0.0055 | 0.0575 | 0.0502 | 7.90x | 9.06x |
| crambin | def2-svp | 6177 | 0.025 | 0.0039 | 0.0036 | 0.0089 | 0.0078 | 0.1074 | 0.0976 | 12.05x | 12.50x |
| crambin | def2-tzvp | 12063 | 0.097 | 0.0041 | 0.0038 | 0.0194 | 0.0157 | 0.3592 | 0.3311 | 18.52x | 21.10x |
| crambin | def2-qzvp | 28167 | 0.409 | 0.0046 | 0.0043 | 0.0564 | 0.0415 | 1.9867 | 1.7735 | 35.23x | 42.77x |
| ubiquitin | def2-svp | 11577 | 0.052 | 0.0109 | 0.0102 | 0.0211 | 0.0186 | 0.3769 | 0.3488 | 17.85x | 18.73x |
| ubiquitin | def2-tzvp | 22442 | 0.200 | 0.0115 | 0.0106 | 0.0417 | 0.0339 | 1.2594 | 1.1611 | 30.17x | 34.20x |
| ubiquitin | def2-qzvp | 53197 | 0.872 | 0.0120 | 0.0111 | 0.1230 | 0.0892 | too large | | | |

### Why there is no fourteen-thread table

**The atom basis pair groups have since been divided into blocks, which gives the
threads something to divide again; see the section on the blocks and the threads
at the end. The single-thread table above still describes the driver, this
subsection no longer does.**

Nothing is left in the driver for the threads to divide, and the numbers say so:
`compute` for crambin in def2-qzvp is 0.0418 warm on fourteen threads against
0.0415 on one, and for ubiquitin in def2-qzvp 0.0891 against 0.0892. What does
change with the threads is the reference driver, which still threads, so the
reference is given at both counts and the driver at one. Warm, seconds.

| molecule | basis | driver 1 thr | driver 14 thr | ref 1 thr | ref 14 thr | against ref, 1 thr | against ref, 14 thr |
|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 0.0012 | 0.0011 | 0.0012 | 0.0005 | 1.06x | 0.42x |
| tagrisso | def2-tzvp | 0.0018 | 0.0018 | 0.0041 | 0.0011 | 2.25x | 0.59x |
| tagrisso | def2-qzvp | 0.0036 | 0.0036 | 0.0201 | 0.0031 | 5.60x | 0.87x |
| taxol | def2-svp | 0.0013 | 0.0013 | 0.0028 | 0.0008 | 2.15x | 0.57x |
| taxol | def2-tzvp | 0.0023 | 0.0023 | 0.0094 | 0.0018 | 4.06x | 0.77x |
| taxol | def2-qzvp | 0.0055 | 0.0054 | 0.0502 | 0.0066 | 9.06x | 1.22x |
| crambin | def2-svp | 0.0078 | 0.0078 | 0.0976 | 0.0129 | 12.50x | 1.66x |
| crambin | def2-tzvp | 0.0157 | 0.0157 | 0.3311 | 0.0417 | 21.10x | 2.65x |
| crambin | def2-qzvp | 0.0415 | 0.0418 | 1.7735 | 0.2044 | 42.77x | 4.89x |
| ubiquitin | def2-svp | 0.0186 | 0.0184 | 0.3488 | 0.0422 | 18.73x | 2.29x |
| ubiquitin | def2-tzvp | 0.0339 | 0.0345 | 1.1611 | 0.1369 | 34.20x | 3.97x |
| ubiquitin | def2-qzvp | 0.0892 | 0.0891 | too large | too large | | |

### What removing the threads cost

The previous commit, `ebd67c3b7`, was built and measured in the same session
rather than compared against the numbers recorded further up, so both sides share
the machine and its state. The ratio is the old time over the new one: above one
the single-threaded version is faster, below one the threaded version was, and
0.22x is the threaded version at four and a half times the speed.

| molecule | basis | 1 thread cold | 1 thread warm | 14 threads cold | 14 threads warm |
|---|---|---|---|---|---|
| tagrisso | def2-svp | 1.03x | 0.99x | 0.93x | 0.89x |
| tagrisso | def2-tzvp | 1.04x | 1.00x | 1.20x | 0.61x |
| tagrisso | def2-qzvp | 1.11x | 0.97x | 0.95x | 0.49x |
| taxol | def2-svp | 1.08x | 1.02x | 0.77x | 0.74x |
| taxol | def2-tzvp | 1.00x | 0.94x | 0.91x | 0.55x |
| taxol | def2-qzvp | 0.94x | 0.97x | 0.79x | 0.35x |
| crambin | def2-svp | 0.99x | 1.01x | 0.53x | 0.44x |
| crambin | def2-tzvp | 1.00x | 0.99x | 0.59x | 0.30x |
| crambin | def2-qzvp | 1.02x | 1.01x | 0.62x | 0.22x |
| ubiquitin | def2-svp | 1.00x | 1.01x | 0.47x | 0.39x |
| ubiquitin | def2-tzvp | 1.00x | 1.02x | 0.55x | 0.30x |
| ubiquitin | def2-qzvp | 1.01x | 1.03x | 0.59x | 0.22x |

### What these numbers say

On a single thread the removal is free. The twelve cases land between 0.94 and
1.11 times, which is the spread of the measurement, and none of them moves in a
way the noise does not explain. That is worth more than it appears: it says the
machinery which is now gone really was gated off when it was given one thread,
and that the sample sort, the predicate which chose it and the thresholds on the
parallel regions cost nothing in the cases where they declined to act.

On fourteen threads the threaded version was 2.3 to 4.6 times faster warm on the
six crambin and ubiquitin cases and 1.1 to 2.9 times on the small ones. That is the price, and
it is paid in full.

Against the reference driver the picture divides by thread count. On one thread
the driver is stronger than anything else recorded in this file, 12.5 to 42.8
times on crambin and ubiquitin, because the reference has no threads either and
the comparison is kernel against kernel. On fourteen threads the reference
catches up and then passes: the driver now loses on all three bases of tagrisso
and on taxol in def2-svp and def2-tzvp, at 0.42 to 0.87 times, where it was ahead
of the reference in every one of those cases before, and its margin on crambin
and ubiquitin falls from between 4.9 and 14.9 times to between 1.7 and 4.9.

Two things do not move with any of this. Ubiquitin in def2-qzvp holds twenty one
gigabytes as a dense matrix, which the reference driver cannot compute on this
machine at all, against 0.872 gigabytes sparse computed in 0.089 seconds. And the
construction of the sparsity pattern is still the floor for a large molecule in a
small basis, 0.0102 of the 0.0186 seconds of ubiquitin in def2-svp.

### A note on the sparsity column

That column is far below the 0.0507 seconds recorded for ubiquitin in def2-svp in
the phases section, 0.0102 warm against it. The gain is not from removing the
threads. It is the radix sort of the bit patterns of the interatomic distances,
which landed between the two measurements, and the phases section predates it.

## The blocks and the threads

The atom basis pair groups are as many as the pairs of the unique atom bases, so
their number is set by the variety of the elements and not by the size of the
molecule: fifteen for crambin, fourteen for ubiquitin, and the largest of them
holds a third of the atom pairs. That is what bounded the earlier attempts at
threading the driver near three times, whatever they parallelised.

The groups are now divided into blocks of a target number of atom pairs, chosen
from the threads and the atom pairs as `npairs / (4 * nthreads)` and no smaller
than 2048, and every stage below works on the blocks: the atom pairs are ordered
one block per thread, the sparsity patterns are described one block per thread,
and the values blocks are computed one block per thread. A single thread leaves
the groups undivided. Crambin in def2-svp gives 15 blocks on one thread, 18 on
two, 26 on four and 65 on fourteen.

### Where the time goes

**Taken at four blocks per thread, which was the target before it was measured.
The table below the sweep supersedes this one.**

Warm, seconds, best of three after one cold call, each case in its own process.
The last two columns are the one thread time over the fourteen thread time.

| molecule | basis | nao | sparsity 1 thr | sparsity 14 thr | compute 1 thr | compute 14 thr | sparsity | compute |
|---|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.0008 | 0.0004 | 0.0011 | 0.0006 | 2.09x | 2.01x |
| tagrisso | def2-tzvp | 1345 | 0.0008 | 0.0004 | 0.0017 | 0.0009 | 2.19x | 2.04x |
| tagrisso | def2-qzvp | 3099 | 0.0010 | 0.0004 | 0.0034 | 0.0011 | 2.30x | 3.11x |
| taxol | def2-svp | 1099 | 0.0008 | 0.0004 | 0.0012 | 0.0007 | 1.83x | 1.79x |
| taxol | def2-tzvp | 2185 | 0.0008 | 0.0004 | 0.0022 | 0.0010 | 1.97x | 2.20x |
| taxol | def2-qzvp | 4947 | 0.0010 | 0.0004 | 0.0051 | 0.0020 | 2.37x | 2.53x |
| crambin | def2-svp | 6177 | 0.0035 | 0.0015 | 0.0074 | 0.0024 | 2.35x | 3.04x |
| crambin | def2-tzvp | 12063 | 0.0037 | 0.0016 | 0.0148 | 0.0040 | 2.33x | 3.69x |
| crambin | def2-qzvp | 28167 | 0.0041 | 0.0017 | 0.0400 | 0.0082 | 2.36x | 4.88x |
| ubiquitin | def2-svp | 11577 | 0.0099 | 0.0036 | 0.0179 | 0.0048 | 2.74x | 3.76x |
| ubiquitin | def2-tzvp | 22442 | 0.0101 | 0.0035 | 0.0331 | 0.0071 | 2.91x | 4.63x |
| ubiquitin | def2-qzvp | 53197 | 0.0107 | 0.0036 | 0.0850 | 0.0154 | 2.96x | 5.52x |

The scaling now rises with the basis, from twice on tagrisso in def2-svp to five
and a half times on ubiquitin in def2-qzvp, which is the ordering the work
predicts: the larger the basis, the more arithmetic there is per atom pair to
divide. The construction of the sparsity pattern scales separately, between 1.8
and 3.0 times, and is what the cheap cases end up waiting on: ubiquitin in
def2-svp spends 0.0036 of its 0.0048 seconds there.

### Against the single threaded driver

| molecule | basis | single threaded | divided, 1 thr | divided, 14 thr | against the single threaded driver |
|---|---|---|---|---|---|
| tagrisso | def2-svp | 0.0012 | 0.0011 | 0.0006 | 2.11x |
| tagrisso | def2-tzvp | 0.0018 | 0.0017 | 0.0009 | 2.12x |
| tagrisso | def2-qzvp | 0.0036 | 0.0034 | 0.0011 | 3.25x |
| taxol | def2-svp | 0.0013 | 0.0012 | 0.0007 | 1.88x |
| taxol | def2-tzvp | 0.0023 | 0.0022 | 0.0010 | 2.33x |
| taxol | def2-qzvp | 0.0055 | 0.0051 | 0.0020 | 2.74x |
| crambin | def2-svp | 0.0078 | 0.0074 | 0.0024 | 3.20x |
| crambin | def2-tzvp | 0.0157 | 0.0148 | 0.0040 | 3.91x |
| crambin | def2-qzvp | 0.0415 | 0.0400 | 0.0082 | 5.06x |
| ubiquitin | def2-svp | 0.0186 | 0.0179 | 0.0048 | 3.91x |
| ubiquitin | def2-tzvp | 0.0339 | 0.0331 | 0.0071 | 4.76x |
| ubiquitin | def2-qzvp | 0.0892 | 0.0850 | 0.0154 | 5.80x |

### The blocks are not what limits this

Per block timings on crambin in def2-qzvp, which divides into 65 blocks on
fourteen threads, with the cost of every block measured separately.

| threads | blocks | summed block work | largest block | work / threads |
|---|---|---|---|---|
| 1 | 15 | 0.0550 | 0.0155 | 0.0550 |
| 2 | 18 | 0.0530 | 0.0104 | 0.0265 |
| 4 | 26 | 0.0648 | 0.0068 | 0.0162 |
| 14 | 65 | 0.1452 | 0.0087 | 0.0104 |

The largest block is 6.1 percent of the work and lies below the work divided by
the threads at every count, so the division itself permits the full fourteen
times and imbalance bounds nothing. The blocks are not the limit.

What the table shows instead is that the same work costs more the more threads
run it: 0.0550 seconds of block work on one thread against 0.1452 on fourteen,
the same blocks inflated 2.6 times. That inflation is the whole of the gap
between the 5 times measured and the 14 times the division allows. One block
spent 76 percent of its time in `make_solid_harmonics` where an identically
shaped sibling spent 4 percent, which points at the per block harmonics, about
1.35 MB at lmax 8, as the contended resource rather than at the kernels.

### A measurement which was wrong

A first sweep of these cases reported the fourteen thread `compute` as slower
than the single thread one on crambin and ubiquitin in def2-qzvp, 0.0447 against
0.0394 and 0.0853 against 0.0846, and a regression was recorded on that basis.
Re-running the same script on the same build gives 0.0082 and 0.0154. The cause
of the first reading was never established; the code was not changed between the
two. It is recorded here because it was reported as a property of the driver and
it was not one.

## How many blocks per thread

The target number of atom pairs of a block is `npairs / (blocks_per_thread *
nthreads)`, no smaller than `min_block_size`, so the two constants decide how many
blocks a molecule is divided into. Both were guessed at four and 2048 when the
division was written. Swept on fourteen threads over crambin and ubiquitin in
def2-svp and def2-qzvp, each point the best of four warm calls, as the geometric
mean of the four cases relative to the best point. Above one is slower.

| blocks per thread | min 512 | min 1024 | min 2048 | min 4096 |
|---|---|---|---|---|
| 2 | 1.01x | 1.00x | 1.00x | 1.01x |
| 4 | 1.04x | 1.04x | 1.05x | 1.03x |
| 8 | 1.32x | 1.31x | 1.26x | 1.15x |
| 16 | 1.57x | 1.52x | 1.40x | 1.22x |
| 32 | 2.03x | 1.70x | 1.47x | 1.21x |
| 64 | 2.40x | 1.92x | 1.46x | 1.20x |

Fewer blocks is better, monotonically, and the effect is large: sixty four blocks
per thread costs between 1.2 and 2.4 times the best. `min_block_size` matters
little at two blocks per thread, where the target is above it anyway, and acts as
a brake at the high end, where it holds the block count down and recovers most of
the loss. Two blocks per thread with a floor of 2048 is the best point and is
what the code now uses, five percent ahead of the four it started with.

The reason is the fixed cost of a block. A block forms its own coordinates and
solid harmonics and bisects the screening once per pair of primitives, none of
which shrinks when the block holds fewer atom pairs. Dividing ubiquitin in
def2-svp into 905 blocks rather than 37 costs 0.0134 seconds against 0.0046, so
the extra 868 blocks cost about ten microseconds of wall time each on fourteen
threads.

This was measured on fourteen threads, where two blocks per thread is 28 blocks
before the per group rounding. The block count follows the threads, so a machine
with 128 cores would get 256 blocks from the same constant, which is the regime
the high end of the sweep emulates: the sweep says that regime is fine as long as
the blocks stay large, and that the floor is what protects it.

### Where the time goes, at two blocks per thread

**Taken before the blocks were ordered by cost and before the check of the values
blocks stopped forming its message on every call. The section on the ordered
blocks at the end supersedes this table.**

Warm, seconds. This supersedes the table of the previous section.

| molecule | basis | nao | sparsity 1 thr | sparsity 14 thr | compute 1 thr | compute 14 thr | sparsity | compute |
|---|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.0008 | 0.0003 | 0.0011 | 0.0006 | 2.32x | 1.88x |
| tagrisso | def2-tzvp | 1345 | 0.0009 | 0.0004 | 0.0017 | 0.0007 | 2.21x | 2.56x |
| tagrisso | def2-qzvp | 3099 | 0.0010 | 0.0004 | 0.0034 | 0.0011 | 2.38x | 2.94x |
| taxol | def2-svp | 1099 | 0.0008 | 0.0004 | 0.0012 | 0.0007 | 1.84x | 1.83x |
| taxol | def2-tzvp | 2185 | 0.0008 | 0.0004 | 0.0022 | 0.0010 | 2.06x | 2.13x |
| taxol | def2-qzvp | 4947 | 0.0010 | 0.0004 | 0.0051 | 0.0017 | 2.22x | 3.05x |
| crambin | def2-svp | 6177 | 0.0035 | 0.0012 | 0.0074 | 0.0021 | 2.78x | 3.48x |
| crambin | def2-tzvp | 12063 | 0.0036 | 0.0013 | 0.0151 | 0.0035 | 2.76x | 4.26x |
| crambin | def2-qzvp | 28167 | 0.0041 | 0.0014 | 0.0404 | 0.0076 | 2.95x | 5.33x |
| ubiquitin | def2-svp | 11577 | 0.0101 | 0.0034 | 0.0181 | 0.0046 | 3.00x | 3.93x |
| ubiquitin | def2-tzvp | 22442 | 0.0103 | 0.0031 | 0.0331 | 0.0074 | 3.27x | 4.48x |
| ubiquitin | def2-qzvp | 53197 | 0.0107 | 0.0033 | 0.0860 | 0.0166 | 3.29x | 5.19x |

### The dense reconstruction

`CSparseMatrix::to_dense` through the `SparseMatrix.to_numpy` binding, warm,
seconds. The zeroing of the dense matrix, the off-diagonal blocks and the
diagonal blocks are each divided among the threads; an element of the dense
matrix belongs to one atom pair and an atom pair to one block, so the blocks
write to disjoint elements and need no synchronization.

| molecule | basis | nao | dense GB | 1 thread | 4 threads | 14 threads | speedup |
|---|---|---|---|---|---|---|---|
| crambin | def2-svp | 6177 | 0.28 | 0.0096 | 0.0033 | 0.0024 | 4.00x |
| crambin | def2-tzvp | 12063 | 1.08 | 0.0419 | 0.0167 | 0.0108 | 3.88x |
| crambin | def2-qzvp | 28167 | 5.91 | 0.4359 | 0.1714 | 0.1411 | 3.09x |
| ubiquitin | def2-svp | 11577 | 1.00 | 0.0214 | 0.0087 | 0.0066 | 3.24x |

The notes above record that the memory bandwidth of this machine saturates at
2.25 times, which would have capped this at about that. It reaches 3.1 to 4.0
times instead. The comparison which makes the point is the largest case: the
whole reconstruction of crambin in def2-qzvp takes 0.141 seconds on fourteen
threads, while numpy takes 0.259 seconds on one thread merely to fill the same
5.91 gigabytes with zeros. The 2.25 figure is a floor for this access pattern and
not a ceiling, first touch spread over the threads doing better than one thread
streaming.

The reconstruction is now much the largest part of the path. Crambin in def2-qzvp
computes in 0.0076 seconds and reconstructs in 0.141, a factor of eighteen.

## The blocks ordered by cost, and what Instruments found

Two changes since the table above. The blocks are now visited from the most
costly to the least, the cost of a block being its atom pairs times its
combinations of basis functions each weighted by the square of the sum of the
angular momenta it carries. And the check guarding the accessors of the values
blocks now forms its message only where it fails, rather than concatenating three
strings on each of the ten thousand calls a `compute` makes.

### Where the time goes

**Taken before the atom pairs of the groups were formed in batches. The section
on the batched groups at the end supersedes this table.**

Warm, seconds, best of three after one cold call, each case in its own process.

| molecule | basis | nao | sparsity 1 thr | sparsity 14 thr | compute 1 thr | compute 14 thr | sparsity | compute |
|---|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 683 | 0.0008 | 0.0004 | 0.0011 | 0.0005 | 1.97x | 2.26x |
| tagrisso | def2-tzvp | 1345 | 0.0008 | 0.0004 | 0.0016 | 0.0008 | 2.12x | 2.10x |
| tagrisso | def2-qzvp | 3099 | 0.0010 | 0.0004 | 0.0030 | 0.0014 | 2.37x | 2.23x |
| taxol | def2-svp | 1099 | 0.0008 | 0.0004 | 0.0012 | 0.0006 | 1.87x | 2.05x |
| taxol | def2-tzvp | 2185 | 0.0008 | 0.0004 | 0.0020 | 0.0008 | 1.94x | 2.47x |
| taxol | def2-qzvp | 4947 | 0.0010 | 0.0003 | 0.0048 | 0.0016 | 3.08x | 3.04x |
| crambin | def2-svp | 6177 | 0.0035 | 0.0013 | 0.0074 | 0.0020 | 2.76x | 3.62x |
| crambin | def2-tzvp | 12063 | 0.0037 | 0.0013 | 0.0145 | 0.0031 | 2.82x | 4.64x |
| crambin | def2-qzvp | 28167 | 0.0041 | 0.0014 | 0.0394 | 0.0069 | 2.99x | 5.74x |
| ubiquitin | def2-svp | 11577 | 0.0099 | 0.0032 | 0.0179 | 0.0044 | 3.15x | 4.09x |
| ubiquitin | def2-tzvp | 22442 | 0.0101 | 0.0031 | 0.0323 | 0.0065 | 3.22x | 4.97x |
| ubiquitin | def2-qzvp | 53197 | 0.0106 | 0.0032 | 0.0840 | 0.0144 | 3.27x | 5.85x |

### Against the reference driver

**Taken before the atom pairs of the groups were formed in batches. The table of
the same name in the section on the batched groups supersedes this one, and in
particular tagrisso in def2-svp is no longer a loss on fourteen threads.**

The same runs, warm, seconds, with the ratio of the reference to the driver.

| molecule | basis | driver 1 thr | reference 1 thr | 1 thr | driver 14 thr | reference 14 thr | 14 thr |
|---|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 0.0011 | 0.0012 | 1.10x | 0.0005 | 0.0004 | 0.92x |
| tagrisso | def2-tzvp | 0.0016 | 0.0040 | 2.49x | 0.0008 | 0.0009 | 1.16x |
| tagrisso | def2-qzvp | 0.0030 | 0.0194 | 6.45x | 0.0014 | 0.0031 | 2.32x |
| taxol | def2-svp | 0.0012 | 0.0027 | 2.21x | 0.0006 | 0.0010 | 1.77x |
| taxol | def2-tzvp | 0.0020 | 0.0091 | 4.48x | 0.0008 | 0.0017 | 2.08x |
| taxol | def2-qzvp | 0.0048 | 0.0481 | 10.12x | 0.0016 | 0.0066 | 4.23x |
| crambin | def2-svp | 0.0074 | 0.0933 | 12.57x | 0.0020 | 0.0124 | 6.03x |
| crambin | def2-tzvp | 0.0145 | 0.3204 | 22.16x | 0.0031 | 0.0447 | 14.35x |
| crambin | def2-qzvp | 0.0394 | 1.7118 | 43.44x | 0.0069 | 0.2042 | 29.76x |
| ubiquitin | def2-svp | 0.0179 | 0.3375 | 18.83x | 0.0044 | 0.0465 | 10.59x |
| ubiquitin | def2-tzvp | 0.0323 | 1.1237 | 34.79x | 0.0065 | 0.1357 | 20.89x |
| ubiquitin | def2-qzvp | 0.0840 | too large |  | 0.0144 | too large |  |

The driver is between 12 and 43 times the reference on a single thread and
between 6 and 30 on fourteen, on the four large cases. Ubiquitin in def2-qzvp has
no reference at all: 53197 basis functions are 21 gigabytes dense, which this
machine cannot hold, against 0.87 gigabytes sparse computed in 0.0144 seconds.
Tagrisso in def2-svp is the one loss on fourteen threads, at 0.92 times. It is
half a millisecond of work and sits below the floor of 2048 atom pairs, so it is
never divided into blocks at all.

The scaling of the driver now rises with the size of the problem throughout,
from 2.1 times on the smallest case to 5.9 times on the largest, which is the
ordering the work predicts and which the earlier attempts at threading inverted.

### What the ordering of the blocks buys

Measured by flipping the ordering within one build, warm, so that nothing but the
order differs.

| molecule | basis | unordered | ordered | gain |
|---|---|---|---|---|
| crambin | def2-qzvp | 0.0079 | 0.0071 | 1.11x |
| ubiquitin | def2-svp | 0.0046 | 0.0043 | 1.07x |
| ubiquitin | def2-tzvp | 0.0075 | 0.0067 | 1.13x |
| ubiquitin | def2-qzvp | 0.0165 | 0.0147 | 1.12x |

The threads draw two or three blocks each, so a costly block drawn last is
finished alone. Per thread instrumentation of the block loop, warm, gives the
idle share of the threads and the ratio of the busiest thread to the mean:

| case | idle, unordered | idle, ordered | imbalance |
|---|---|---|---|
| crambin def2-svp | 39.1% | 25.4% | 1.57 to 1.28 |
| crambin def2-qzvp | 25.6% | 15.7% | 1.34 to 1.18 |
| ubiquitin def2-svp | 28.1% | 23.5% | 1.36 to 1.27 |
| ubiquitin def2-qzvp | 26.7% | 12.0% | 1.36 to 1.13 |

### What Instruments says

Recorded with `xctrace`, ubiquitin in def2-qzvp on fourteen threads.

The System Trace, over eight `compute` calls, gives a mean concurrency of **9.0
of the 14 threads**, so the machine is 64 percent used, with 51 percent of the
window at twelve threads or more and 19 percent below two. The parallel regions
reach the full width of the machine; what is left is the serial remainder around
them, about 3.7 milliseconds of a 15 millisecond call.

The Time Profiler puts `_simd_exp_d2` at the top of the self time at 4.0 percent,
the kernels together at about 12, the allocator at about 5.5 across
`_xzm_malloc_large_huge`, `_xzm_free`, `xzm_segment_group_free_chunk` and
`malloc_type_aligned_alloc`, and the check of the values blocks at 1.0. The last
of those is now fixed and is worth 1.021 times as the geometric mean of five
cases, every one of them improving.

The allocator share is the per block coordinates and solid harmonics, allocated
and freed once per block, and it is the next thing worth attacking. It needs
buffers carried across the blocks rather than a small fix.

### What did not work

| what was tried | result |
|---|---|
| blocks sized so that each carries the same cost, the atom pairs of a group scaled by the angular momenta of its atom bases | 0.98 times, and no exponent from zero to three of the momentum sum reached 1.05 |
| the atom pairs of `to_dense` tiled, so that the writes of a tile fall in fewer rows | a loss at every tile size from 32 to 4096 |
| the values of `CSimdMatrix` taken from a pool held by the thread, so that the coordinates and the solid harmonics of a block reuse those of the block before | 1.017 times, for a custom allocator in a core class and 41 megabytes retained |

The equal cost sizing fails because it gives the most costly groups the smallest
blocks, and those are the groups whose fixed cost per block, chiefly the solid
harmonics formed to the sum of their angular momenta, is the largest. The tiling
of `to_dense` fails because the scatter is already at the limit of the memory
system for random access, 203 gigabytes per second of line traffic against the
352 the machine reaches on a linear write, so there was no locality left to
recover.

The pool was built on the reading that the allocator is 5.5 percent of the
samples, and about 1.05 times was expected of it. It measured 1.017 times, every
one of five cases improving but none by much, and it was removed. The reading was
right and the inference from it was wrong: most of that 5.5 percent is the values
blocks of the matrix, 0.87 gigabytes allocated and freed for every `compute` of
ubiquitin in def2-qzvp through `allocate` and the destructor, which are plain
`new double[]` and never reach `CSimdMatrix`. The pool could only ever reach the
coordinates and the solid harmonics. Peak resident memory went from 1.238 to
1.279 gigabytes with it.

Reusing the values blocks is the part which would pay, and it cannot be done
under the allocator: it needs the driver to fill a matrix it is given rather than
return a new one, so that a caller computing repeatedly keeps one.

## The atom pairs of the groups formed in batches

Instruments had shown a fifth of every `compute` running on one thread, and the
phases of the sparsity were timed to find it. It was not the division of the
groups into blocks, which is 0.6 percent, and not the description of the sparsity
patterns, which is threaded. It was `CMolecularBasis::basis_pair_groups`, which
formed every atom pair of every group with `push_back` in nested loops, before
anything was divided at all: 0.00168 seconds of a 0.0044 second call for
ubiquitin in def2-svp, **37 percent of it**, and 12 percent of the same molecule
in def2-qzvp.

The atom pairs are now formed in batches by the threads. A group is created with
its atom pairs sized but not formed, as their number follows from the atoms of the
atom basis groups it pairs and needs no enumeration, and a pool of chunks
spanning all the groups is then filled by the threads, each chunk finding its atom
pairs from their indices alone. Within one molecular basis an atom carries one
atom basis, so the atom basis groups share no atom: the atom pairs of a symmetric
group are the strict upper triangle of its atoms and those of a pair of groups are
the full rectangle of theirs, and both invert in closed form. The two molecular
bases factory, where an atom does appear on both sides, is left as it was.

### What it bought

Warm, fourteen threads, seconds.

| molecule | basis | sparsity before | sparsity after | compute before | compute after | gain |
|---|---|---|---|---|---|---|
| tagrisso | def2-svp | 0.0004 | 0.0004 | 0.0005 | 0.0005 | 0.95x |
| tagrisso | def2-tzvp | 0.0004 | 0.0003 | 0.0008 | 0.0006 | 1.22x |
| tagrisso | def2-qzvp | 0.0004 | 0.0004 | 0.0014 | 0.0011 | 1.22x |
| taxol | def2-svp | 0.0004 | 0.0005 | 0.0006 | 0.0006 | 0.95x |
| taxol | def2-tzvp | 0.0004 | 0.0004 | 0.0008 | 0.0009 | 0.90x |
| taxol | def2-qzvp | 0.0003 | 0.0004 | 0.0016 | 0.0016 | 0.96x |
| crambin | def2-svp | 0.0013 | 0.0009 | 0.0020 | 0.0016 | 1.29x |
| crambin | def2-tzvp | 0.0013 | 0.0009 | 0.0031 | 0.0029 | 1.09x |
| crambin | def2-qzvp | 0.0014 | 0.0010 | 0.0069 | 0.0065 | 1.06x |
| ubiquitin | def2-svp | 0.0032 | 0.0016 | 0.0044 | 0.0028 | 1.59x |
| ubiquitin | def2-tzvp | 0.0031 | 0.0016 | 0.0065 | 0.0050 | 1.30x |
| ubiquitin | def2-qzvp | 0.0032 | 0.0018 | 0.0144 | 0.0127 | 1.14x |

The sparsity phase roughly halves on the large cases, which is the serial
enumeration going away, and the whole of `compute` gains between 1.06 and 1.59
times. The molecules under a millisecond lose between two and ten percent, which
is the cost of forming the chunks where there is too little work to repay it.

### Where the time goes

Warm, seconds. This supersedes the table of the previous section.

| molecule | basis | nao | sparsity 1 thr | sparsity 14 thr | compute 1 thr | compute 14 thr | sparsity | compute |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-svp | 683 | 0.0008 | 0.0004 | 0.0010 | 0.0005 | 1.81x | 2.21x |
| tagrisso | def2-tzvp | 1345 | 0.0008 | 0.0004 | 0.0016 | 0.0008 | 1.96x | 2.00x |
| tagrisso | def2-qzvp | 3099 | 0.0010 | 0.0005 | 0.0029 | 0.0011 | 2.09x | 2.71x |
| taxol | def2-svp | 1099 | 0.0007 | 0.0005 | 0.0012 | 0.0005 | 1.62x | 2.23x |
| taxol | def2-tzvp | 2185 | 0.0008 | 0.0005 | 0.0020 | 0.0009 | 1.80x | 2.25x |
| taxol | def2-qzvp | 4947 | 0.0010 | 0.0005 | 0.0046 | 0.0016 | 2.02x | 2.77x |
| crambin | def2-svp | 6177 | 0.0032 | 0.0009 | 0.0070 | 0.0016 | 3.57x | 4.48x |
| crambin | def2-tzvp | 12063 | 0.0034 | 0.0010 | 0.0141 | 0.0032 | 3.43x | 4.40x |
| crambin | def2-qzvp | 28167 | 0.0039 | 0.0010 | 0.0386 | 0.0065 | 3.89x | 5.95x |
| ubiquitin | def2-svp | 11577 | 0.0089 | 0.0016 | 0.0167 | 0.0028 | 5.66x | 5.92x |
| ubiquitin | def2-tzvp | 22442 | 0.0092 | 0.0017 | 0.0312 | 0.0050 | 5.53x | 6.19x |
| ubiquitin | def2-qzvp | 53197 | 0.0097 | 0.0018 | 0.0832 | 0.0128 | 5.44x | 6.49x |

The scaling reaches 5.9 to 6.5 times on the three ubiquitin cases and rises with
the size of the problem throughout.

### Against the reference driver

| molecule | basis | driver 1 thr | reference 1 thr | 1 thr | driver 14 thr | reference 14 thr | 14 thr |
| --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-svp | 0.0010 | 0.0011 | 1.09x | 0.0005 | 0.0004 | 0.84x |
| tagrisso | def2-tzvp | 0.0016 | 0.0038 | 2.46x | 0.0008 | 0.0015 | 1.89x |
| tagrisso | def2-qzvp | 0.0029 | 0.0185 | 6.29x | 0.0011 | 0.0031 | 2.86x |
| taxol | def2-svp | 0.0012 | 0.0026 | 2.20x | 0.0005 | 0.0009 | 1.69x |
| taxol | def2-tzvp | 0.0020 | 0.0087 | 4.34x | 0.0009 | 0.0017 | 1.92x |
| taxol | def2-qzvp | 0.0046 | 0.0471 | 10.32x | 0.0016 | 0.0067 | 4.05x |
| crambin | def2-svp | 0.0070 | 0.0921 | 13.23x | 0.0016 | 0.0125 | 8.06x |
| crambin | def2-tzvp | 0.0141 | 0.3154 | 22.44x | 0.0032 | 0.0401 | 12.53x |
| crambin | def2-qzvp | 0.0386 | 1.7087 | 44.29x | 0.0065 | 0.1990 | 30.70x |
| ubiquitin | def2-svp | 0.0167 | 0.3284 | 19.68x | 0.0028 | 0.0414 | 14.70x |
| ubiquitin | def2-tzvp | 0.0312 | 1.1116 | 35.64x | 0.0050 | 0.1344 | 26.68x |
| ubiquitin | def2-qzvp | 0.0832 | 6.2043 | 74.60x | 0.0128 | 0.6904 | 53.91x |

Every case beats the reference on a single thread and all but one of them on
fourteen. Tagrisso in def2-svp is the exception at 0.84 times, which is half a
millisecond of work against four tenths of one, and it has moved either side of
parity between runs: an earlier measurement of the same build put it at 1.06.

The margin runs to 75 times on a single thread and 54 on fourteen, both of them
ubiquitin in def2-qzvp. That case had no reference at all in an earlier version of
this table, on the assumption that its 21 gigabytes as a dense matrix would not
fit. The reference driver does not store the full square, and the case in fact
runs in 13.5 gigabytes of resident memory without swapping, so it is measured
here: 6.2043 seconds on one thread and 0.6904 on fourteen, against 0.0832 and
0.0128 for 0.87 gigabytes sparse.

### What is left

The reconstruction of the dense matrix is now 22 times the integrals it
reconstructs, 0.1406 seconds against 0.0065 for crambin in def2-qzvp, and the
notes above record that it is at the limit of the memory system on both of its
halves. The integrals themselves are no longer where the time of this path goes.

## Why the small molecules are not divided

`min_block_size` is an absolute count of atom pairs, 2048, and tagrisso holds
2415 atom pairs in all with its largest group holding 924. Both are below the
floor, so tagrisso is never divided into blocks whatever the threads, and taxol
divides only its largest group. The two of them scale 2.1 to 3.2 times on
fourteen threads where ubiquitin reaches 6.9.

Expressing the floor in atom pairs is also wrong in principle. What it guards is
the cost a block carries whatever it holds, chiefly the bisection of the screening
over the pairs of primitives of every combination of basis functions, and that
cost is set by the atom bases rather than by a count of atom pairs: the same floor
is far too high for def2-qzvp, where an atom pair is worth many times more, and
about right for def2-svp on a protein, which is the case it was fitted on.

So the floor was made a floor in work instead, `min_block_work / w`, with `w` the
weight of an atom pair of the group, the same expression the ordering of the
blocks uses. Swept on fourteen threads over the molecules where the floor binds,
warm, seconds, with the blocks in brackets.

| min work | tagrisso def2-svp | tagrisso def2-qzvp | taxol def2-svp | taxol def2-qzvp |
|---|---|---|---|---|
| no division | 0.0005 (10) | 0.0010 (10) | 0.0006 (9) | 0.0017 (9) |
| 73728 | 0.0006 (12) | 0.0013 (34) | 0.0007 (16) | 0.0015 (34) |
| 36864 | 0.0007 (15) | 0.0013 (34) | 0.0008 (22) | 0.0015 (34) |
| 18432 | 0.0006 (21) | 0.0012 (34) | 0.0007 (30) | 0.0015 (34) |
| 9216 | 0.0007 (27) | 0.0013 (34) | 0.0010 (32) | 0.0016 (34) |
| 4608 | 0.0009 (31) | 0.0013 (34) | 0.0011 (34) | 0.0016 (34) |
| 2304 | 0.0007 (34) | 0.0012 (34) | 0.0008 (34) | 0.0015 (34) |

Dividing them makes three of the four slower, at every one of the six settings,
so the direction is not noise. Only taxol in def2-qzvp gains, and by about a
tenth. The change was reverted.

The premise was right and the conclusion does not follow from it. The small
molecules are not merely prevented from dividing, they are not worth dividing: a
call of half a millisecond cannot repay the cost a block carries, however the
floor is expressed. The flat floor of 2048 atom pairs happens to produce the
behaviour the sweep finds best for these molecules, for a reason which is not the
one its comment gives.

This is the second block sizing heuristic to lose to the flat rule, after the
blocks sized to carry equal cost recorded above. Both were built on the same
reasoning, that the work per atom pair differs by the atom bases and the blocks
should follow it, and both were beaten by a constant. A third attempt should
weigh that the cost a block carries is large enough that any heuristic which
makes more blocks starts behind.

## The diffuse basis sets

The same four molecules in def2-svpd, def2-tzvpd and def2-qzvpd. Diffuse
functions carry small exponents and so reach far, which is what the screening
lives on, so these are the cases where the bound of the driver discards the least.

### Where the time goes

Warm, seconds.

| molecule | basis | nao | sparse GB | sparsity 14 thr | compute 1 thr | compute 14 thr | scaling | reference 14 thr | against it |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-svpd | 1010 | 0.003 | 0.0004 | 0.0013 | 0.0006 | 2.12x | 0.0006 | 0.94x |
| tagrisso | def2-tzvpd | 1672 | 0.007 | 0.0004 | 0.0020 | 0.0009 | 2.10x | 0.0012 | 1.27x |
| tagrisso | def2-qzvpd | 3426 | 0.024 | 0.0005 | 0.0038 | 0.0013 | 2.88x | 0.0034 | 2.57x |
| taxol | def2-svpd | 1657 | 0.007 | 0.0005 | 0.0017 | 0.0007 | 2.49x | 0.0010 | 1.43x |
| taxol | def2-tzvpd | 2743 | 0.018 | 0.0004 | 0.0029 | 0.0012 | 2.30x | 0.0025 | 1.98x |
| taxol | def2-qzvpd | 5505 | 0.057 | 0.0005 | 0.0064 | 0.0022 | 2.96x | 0.0078 | 3.57x |
| crambin | def2-svpd | 9294 | 0.106 | 0.0009 | 0.0139 | 0.0026 | 5.24x | 0.0234 | 8.84x |
| crambin | def2-tzvpd | 15180 | 0.262 | 0.0009 | 0.0272 | 0.0053 | 5.12x | 0.0611 | 11.48x |
| crambin | def2-qzvpd | 31284 | 0.791 | 0.0011 | 0.0674 | 0.0110 | 6.15x | 0.2393 | 21.85x |
| ubiquitin | def2-svpd | 17433 | 0.240 | 0.0017 | 0.0330 | 0.0052 | 6.38x | 0.0787 | 15.20x |
| ubiquitin | def2-tzvpd | 28298 | 0.593 | 0.0017 | 0.0629 | 0.0092 | 6.85x | 0.1889 | 20.56x |
| ubiquitin | def2-qzvpd | 59053 | 1.821 | 0.0018 | 0.1541 | 0.0229 | 6.74x | 0.8266 | 36.17x |

### What the diffuse functions cost

The basis grows by half and the sparse matrix by four and a half.

| molecule | basis | nao | nao diffuse | sparse GB | sparse GB diffuse | growth | compute 14 thr | diffuse | growth |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-svp | 683 | 1010 | 0.001 | 0.003 | 2.9x | 0.0005 | 0.0006 | 1.3x |
| tagrisso | def2-tzvp | 1345 | 1672 | 0.004 | 0.007 | 1.9x | 0.0008 | 0.0009 | 1.2x |
| tagrisso | def2-qzvp | 3099 | 3426 | 0.017 | 0.024 | 1.5x | 0.0011 | 0.0013 | 1.2x |
| taxol | def2-svp | 1099 | 1657 | 0.002 | 0.007 | 3.2x | 0.0005 | 0.0007 | 1.3x |
| taxol | def2-tzvp | 2185 | 2743 | 0.009 | 0.018 | 2.1x | 0.0009 | 0.0012 | 1.4x |
| taxol | def2-qzvp | 4947 | 5505 | 0.037 | 0.057 | 1.5x | 0.0016 | 0.0022 | 1.3x |
| crambin | def2-svp | 6177 | 9294 | 0.025 | 0.106 | 4.2x | 0.0016 | 0.0026 | 1.7x |
| crambin | def2-tzvp | 12063 | 15180 | 0.097 | 0.262 | 2.7x | 0.0032 | 0.0053 | 1.7x |
| crambin | def2-qzvp | 28167 | 31284 | 0.409 | 0.791 | 1.9x | 0.0065 | 0.0110 | 1.7x |
| ubiquitin | def2-svp | 11577 | 17433 | 0.052 | 0.240 | 4.6x | 0.0028 | 0.0052 | 1.8x |
| ubiquitin | def2-tzvp | 22442 | 28298 | 0.200 | 0.593 | 3.0x | 0.0050 | 0.0092 | 1.8x |
| ubiquitin | def2-qzvp | 53197 | 59053 | 0.872 | 1.821 | 2.1x | 0.0128 | 0.0229 | 1.8x |

Ubiquitin in def2-svp against def2-svpd is the clearest of them: 11577 basis
functions become 17433, a factor of 1.5, while the sparse matrix goes from 0.052
to 0.240 gigabytes, a factor of 4.6, and `compute` from 0.0028 to 0.0052 seconds.
The extra cost is not the basis, it is the atom pairs which the screening no
longer discards and the combinations of basis functions which survive on the atom
pairs it keeps.

### What these numbers say

The scaling holds and improves a little, 6.74 times on ubiquitin in def2-qzvpd
against 6.49 for the same molecule in def2-qzvp, and 6.85 on ubiquitin in
def2-tzvpd is the best of these twelve. More surviving atom pairs mean more work
in a block and a smaller share for the cost a block carries whatever it holds.

The advantage over the reference narrows, as it should. Crambin in def2-qzvp is
30.7 times on fourteen threads and in def2-qzvpd 21.9. The driver wins on the
screening, and diffuse functions are precisely what defeats screening, so the
case where it discards the least is the case where the reference closes the most
ground. Tagrisso in def2-svpd is the one loss, 0.94 times on six tenths of a
millisecond.

The construction of the sparsity pattern does not move at all: every ubiquitin
case is 0.0017 or 0.0018 seconds on fourteen threads whatever the basis, as it
follows the atoms of the molecule and not the functions on them.

Ubiquitin in def2-qzvpd is the largest case of this section, 59053 basis
functions, 1.82 gigabytes sparse against 26 dense, computed in 0.0229 seconds.
The reference driver holds it after all, in 16.7 gigabytes resident and without
swapping, and takes 0.8266 seconds over it.

## The correlation consistent basis sets

The four molecules in cc-pVXZ and aug-cc-pVXZ from double to sextuple zeta, forty
cases. These reach angular momentum six, so the recursions of the integrals reach
twelve, and they are the only cases in these notes which exercise the solid
harmonics above order eight.

### A bug they found

The first attempt at this table failed on thirteen of the forty cases with a
memory error, which AddressSanitizer placed exactly. Every one of the forty nine
kernels sized the buffer its primitives accumulate into as `dimensions.back()`,
on the assumption that the last pair of primitives is the one reaching furthest.
The primitives are sorted by descending exponent, so that holds if the reach
followed the decay alone, but the bound of a pair of primitives carries their
prefactor as well, and for exponents which are close the prefactor decides. On
crambin in aug-cc-pvdz, an eight primitive contraction against one diffuse
primitive gives `back` of 48787 against a true largest of 49175, and the loop
then writes some three kilobytes past its buffer.

The values were never wrong: the atom pairs between the two are ones whose
integrals fall below the screening threshold, and the reference driver agrees to
1.0e-14 with the bug in place as it does without it. The def2 sets never trip it.
It reproduces on one thread, so it was never a race; the threads only decided
whether the stray write landed somewhere the allocator noticed.

The kernels now take the largest of the dimensions rather than the last, and
AddressSanitizer is clean on the cases which reported the overflow.

### Where the time goes

Warm, seconds, best of three after one cold call, each case in its own process.
The reference driver is run only where the dense matrix fits in memory.

| molecule | basis | nao | dense GB | sparse GB | compute 1 thr | compute 14 thr | scaling | reference 14 thr | against it |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | cc-pvdz | 683 | 0.0 | 0.00 | 0.0013 | 0.0005 | 2.36x | 0.0009 | 1.56x |
| tagrisso | cc-pvtz | 1572 | 0.0 | 0.00 | 0.0018 | 0.0008 | 2.39x | 0.0015 | 1.99x |
| tagrisso | cc-pvqz | 3025 | 0.1 | 0.02 | 0.0029 | 0.0011 | 2.76x | 0.0033 | 3.13x |
| tagrisso | cc-pv5z | 5182 | 0.2 | 0.04 | 0.0057 | 0.0018 | 3.07x | 0.0087 | 4.71x |
| tagrisso | cc-pv6z | 8183 | 0.5 | 0.10 | 0.0116 | 0.0036 | 3.19x | 0.0235 | 6.43x |
| tagrisso | aug-cc-pvdz | 1148 | 0.0 | 0.00 | 0.0018 | 0.0007 | 2.59x | 0.0011 | 1.56x |
| tagrisso | aug-cc-pvtz | 2461 | 0.0 | 0.02 | 0.0030 | 0.0011 | 2.82x | 0.0022 | 2.07x |
| tagrisso | aug-cc-pvqz | 4478 | 0.1 | 0.05 | 0.0057 | 0.0019 | 3.06x | 0.0058 | 3.10x |
| tagrisso | aug-cc-pv5z | 7339 | 0.4 | 0.12 | 0.0124 | 0.0041 | 3.03x | 0.0179 | 4.36x |
| tagrisso | aug-cc-pv6z | 11184 | 0.9 | 0.24 | 0.0261 | 0.0085 | 3.08x | 0.0445 | 5.24x |
| taxol | cc-pvdz | 1099 | 0.0 | 0.00 | 0.0015 | 0.0007 | 2.26x | 0.0013 | 1.87x |
| taxol | cc-pvtz | 2516 | 0.0 | 0.01 | 0.0024 | 0.0009 | 2.59x | 0.0027 | 2.81x |
| taxol | cc-pvqz | 4825 | 0.2 | 0.03 | 0.0043 | 0.0015 | 2.80x | 0.0064 | 4.16x |
| taxol | cc-pv5z | 8246 | 0.5 | 0.09 | 0.0097 | 0.0032 | 3.01x | 0.0192 | 5.92x |
| taxol | cc-pv6z | 12999 | 1.3 | 0.21 | 0.0213 | 0.0070 | 3.03x | 0.0554 | 7.90x |
| taxol | aug-cc-pvdz | 1844 | 0.0 | 0.01 | 0.0025 | 0.0010 | 2.59x | 0.0015 | 1.56x |
| taxol | aug-cc-pvtz | 3933 | 0.1 | 0.04 | 0.0049 | 0.0017 | 2.92x | 0.0041 | 2.46x |
| taxol | aug-cc-pvqz | 7134 | 0.4 | 0.11 | 0.0107 | 0.0035 | 3.07x | 0.0127 | 3.65x |
| taxol | aug-cc-pv5z | 11667 | 1.0 | 0.26 | 0.0238 | 0.0072 | 3.32x | 0.0387 | 5.39x |
| taxol | aug-cc-pv6z | 17752 | 2.3 | 0.54 | 0.0527 | 0.0147 | 3.58x | 0.1077 | 7.31x |
| crambin | cc-pvdz | 6177 | 0.3 | 0.03 | 0.0089 | 0.0018 | 4.83x | 0.0189 | 10.23x |
| crambin | cc-pvtz | 14244 | 1.5 | 0.11 | 0.0172 | 0.0031 | 5.50x | 0.0602 | 19.19x |
| crambin | cc-pvqz | 27459 | 5.6 | 0.35 | 0.0357 | 0.0058 | 6.20x | 0.1964 | 34.11x |
| crambin | cc-pv5z | 47106 | 16.5 | 0.95 | 0.0828 | 0.0125 | 6.60x | 0.6280 | 50.07x |
| crambin | cc-pv6z | 74469 | 41.3 | 2.21 | 0.1870 | 0.0273 | 6.85x | 1.8725 | 68.54x |
| crambin | aug-cc-pvdz | 10380 | 0.8 | 0.17 | 0.0261 | 0.0040 | 6.48x | 0.0372 | 9.25x |
| crambin | aug-cc-pvtz | 22311 | 3.7 | 0.60 | 0.0567 | 0.0085 | 6.69x | 0.1201 | 14.17x |
| crambin | aug-cc-pvqz | 40674 | 12.3 | 1.61 | 0.1284 | 0.0181 | 7.09x | 0.4051 | 22.37x |
| crambin | aug-cc-pv5z | 66753 | 33.2 | 3.73 | 0.2906 | 0.0415 | 7.01x | 1.2618 | 30.42x |
| crambin | aug-cc-pv6z | 101832 | 77.3 | 7.63 | 0.6109 | 0.0858 | 7.12x | too large |  |
| ubiquitin | cc-pvdz | 11577 | 1.0 | 0.06 | 0.0196 | 0.0032 | 6.20x | 0.0634 | 20.09x |
| ubiquitin | cc-pvtz | 26870 | 5.4 | 0.24 | 0.0360 | 0.0055 | 6.59x | 0.1946 | 35.67x |
| ubiquitin | cc-pvqz | 51984 | 20.1 | 0.74 | 0.0739 | 0.0107 | 6.92x | 0.6708 | 62.77x |
| ubiquitin | cc-pv5z | 89381 | 59.5 | 2.03 | 0.1739 | 0.0259 | 6.72x | too large |  |
| ubiquitin | cc-pv6z | 141523 | 149.2 | 4.74 | 0.3931 | 0.0619 | 6.35x | too large |  |
| ubiquitin | aug-cc-pvdz | 19511 | 2.8 | 0.42 | 0.0649 | 0.0091 | 7.16x | 0.1228 | 13.54x |
| ubiquitin | aug-cc-pvtz | 42163 | 13.2 | 1.48 | 0.1419 | 0.0192 | 7.41x | 0.4242 | 22.15x |
| ubiquitin | aug-cc-pvqz | 77098 | 44.3 | 3.95 | 0.3105 | 0.0461 | 6.74x | too large |  |
| ubiquitin | aug-cc-pv5z | 126778 | 119.8 | 9.09 | 0.6977 | 0.1079 | 6.47x | too large |  |
| ubiquitin | aug-cc-pv6z | 193665 | 279.4 | 18.53 | 1.6873 | 0.2354 | 7.17x | too large |  |

### How much of this is noise

The whole table was measured twice, on the same build, to see what a single run
is worth.

| | median | worst |
|---|---|---|
| compute, fourteen threads | 1.7% | 9.5% |
| compute, one thread | 2.8% | 7.6% |

A single run therefore says nothing below about five percent, which is worth
remembering when reading any one number in these notes. The differences which
this session acted on were measured by flipping a switch within one build, and
one which was not, a supposed regression of the parallel block loop, turned out
to be this noise and nothing else.

### What these numbers say

The scaling holds from 2.3 to 3.6 times on the two small molecules and from 4.8
to 7.4 on crambin and ubiquitin, the best of them being ubiquitin in aug-cc-pvtz
at 7.4. The small molecules gain over their def2 counterparts, where they reached
2.0 to 2.8, because the higher angular momenta give a block more work to carry
against the cost it pays whatever it holds.

Six of the forty are beyond what the reference driver can hold. The dense GB
column overstates what it actually needs, its resident memory coming out at
about two thirds of that, so several cases which an earlier version of this table
gave up on are measured here rather than assumed. The smallest of the six,
ubiquitin in aug-cc-pvqz, was attempted and then rejected: it reached 27
gigabytes resident and spent 60 seconds of system time against 73 of user time,
which is the memory system and not the driver. The other five are ninety thousand
basis functions and above and were not attempted.

Where it can be compared the reference is between 1.6 and 68.5 times slower, the
widest margin being crambin in cc-pv6z, 1.8725 seconds against 0.0273. That is
the largest reference measured anywhere in these notes, 24.6 gigabytes resident
and no swapping, and it is the only accepted one whose system time reaches a sixth of
its user time.

The largest case of these notes is ubiquitin in aug-cc-pv6z, 193665 basis
functions, 279 gigabytes as a dense matrix against 18.5 sparse, computed in 0.24
seconds on fourteen threads and 1.69 on one.

## The dense reconstruction, and where its time actually goes

The reconstruction became the largest part of this path once the integrals were
threaded, twenty times `compute` on crambin in def2-qzvp, so it was instrumented
in place: a timer around each of its parts and the minor faults of the process
counted across the zeroing.

| part | crambin def2-qzvp, fourteen threads |
|---|---|
| zeroing the dense matrix | 0.0861 s, and 387393 minor faults |
| the off-diagonal blocks | 0.0527 s, 109700914 values written |
| the diagonal blocks | 0.0002 s |

The two halves are limited by different things and neither is the code.

The scatter writes 109700914 values of eight bytes each, and each of them pulls a
cache line of sixty four, so it moves 7.0 gigabytes in 0.0527 seconds, 133
gigabytes per second. A linear write on this machine reaches 352. That is the
limit of the memory system for writes scattered over six gigabytes, not a limit
of the loop, and the attempt to improve its locality by tiling, recorded above,
lost at every tile size.

The zeroing faults in every page of the dense matrix: 6.35 gigabytes over pages
of sixteen kilobytes is 387891 pages against the 387393 faults counted. The same
zeroing of a warm buffer takes 0.0168 seconds. It is therefore not bandwidth but
the first touch of newly mapped memory.

### The zeroing is not redundant work, it is the allocation

The obvious idea, to allocate the array zeroed and skip the zeroing, was measured
and dropped. A fresh allocation with a parallel fill takes 0.1099 seconds and a
fresh zeroed allocation touched once per page takes 0.1084: writing the zeros
costs 1.5 milliseconds once the fault is paid, one percent of the whole
reconstruction. Every page is touched by the scatter in any case, 10 to 18
percent of the elements being non-zero but spread so that no page of any of the
cases is empty, so there is nothing for a lazily zeroed page to save.

### What does help is not allocating

The faults are paid on every call only because the binding hands numpy a new
array each time.

| case | dense | first call | later calls | faults per call |
|---|---|---|---|---|
| crambin def2-tzvp | 1.08 GB | 0.0287 | 0.0112 | 71064 then 0 |
| crambin def2-qzvp | 5.91 GB | 0.1414 | 0.1362, 0.1306 | 387393 every time |

Below about a gigabyte the allocator keeps the block when it is freed and the
second and later reconstructions cost nothing in faults. Above it the block is
returned to the system and every call faults the whole matrix in again.

`fill_numpy` reconstructs into an array the caller keeps, so that a caller
converting one matrix after another pays the faults once.

| case | dense | to_numpy | fill_numpy, kept buffer | gain | faults |
|---|---|---|---|---|---|
| crambin def2-tzvp | 1.08 GB | 0.0116 | 0.0114 | 1.01x | 0 to 1 |
| ubiquitin def2-svp | 1.00 GB | 0.0069 | 0.0069 | 0.99x | 0 to 0 |
| crambin def2-qzvp | 5.91 GB | 0.1189 | 0.0535 | 2.22x | 387396 to 0 |

The results are bit-identical to `to_numpy`. The gain appears only where there
were faults to remove, which is the point: the two smaller cases are unchanged
because the allocator was already keeping their memory.

What is left of the reconstruction on the large case is 0.0535 seconds, which is
the scatter and almost nothing else, and the scatter is at the limit quoted
above. A single reconstruction of a matrix which has never been built cannot be
made faster than the cost of faulting in its memory.

## The kinetic energy driver

The kinetic energy integrals follow the overlap: the same sparse matrix, the same
division of the atom basis pair groups into blocks, the same ordering of the
blocks by cost and the same parallel loop over them. Only the kernels differ, and
they carry more arithmetic for the same atom pair.

The reference driver of VeloxChem implements the kinetic energy to angular
momentum four alone, so it is compared against only where the basis stays within
that. The column distinguishes the two reasons it can be missing rather than
being left blank: `not valid` where the basis reaches angular momentum five or
six, and `too large` where the reference could not be held in memory.

### The def2 sets

Warm, seconds, best of three after one cold call, each case in its own process.

| molecule | basis | nao | lmax | sparse GB | compute 1 thr | compute 14 thr | scaling | reference 14 thr | against it |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-svp | 683 | 2 | 0.00 | 0.0011 | 0.0005 | 2.20x | 0.0005 | 1.07x |
| tagrisso | def2-tzvp | 1345 | 3 | 0.00 | 0.0017 | 0.0008 | 2.03x | 0.0011 | 1.36x |
| tagrisso | def2-qzvp | 3099 | 4 | 0.02 | 0.0036 | 0.0013 | 2.85x | 0.0037 | 2.93x |
| taxol | def2-svp | 1099 | 2 | 0.00 | 0.0013 | 0.0006 | 2.22x | 0.0008 | 1.39x |
| taxol | def2-tzvp | 2185 | 3 | 0.01 | 0.0022 | 0.0010 | 2.23x | 0.0022 | 2.17x |
| taxol | def2-qzvp | 4947 | 4 | 0.04 | 0.0060 | 0.0022 | 2.77x | 0.0087 | 4.05x |
| crambin | def2-svp | 6177 | 2 | 0.03 | 0.0075 | 0.0016 | 4.63x | 0.0144 | 8.85x |
| crambin | def2-tzvp | 12063 | 3 | 0.10 | 0.0167 | 0.0032 | 5.26x | 0.0493 | 15.56x |
| crambin | def2-qzvp | 28167 | 4 | 0.44 | 0.0522 | 0.0080 | 6.50x | 0.2625 | 32.71x |
| ubiquitin | def2-svp | 11577 | 2 | 0.06 | 0.0178 | 0.0030 | 6.03x | 0.0488 | 16.51x |
| ubiquitin | def2-tzvp | 22442 | 3 | 0.22 | 0.0361 | 0.0057 | 6.38x | 0.1658 | 29.30x |
| ubiquitin | def2-qzvp | 53197 | 4 | 0.95 | 0.1116 | 0.0155 | 7.22x | 0.9114 | 58.91x |

### The def2 sets with diffuse functions

| molecule | basis | nao | lmax | sparse GB | compute 1 thr | compute 14 thr | scaling | reference 14 thr | against it |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-svpd | 1010 | 2 | 0.00 | 0.0014 | 0.0006 | 2.38x | 0.0007 | 1.14x |
| tagrisso | def2-tzvpd | 1672 | 3 | 0.01 | 0.0021 | 0.0009 | 2.48x | 0.0013 | 1.50x |
| tagrisso | def2-qzvpd | 3426 | 4 | 0.03 | 0.0046 | 0.0016 | 2.89x | 0.0043 | 2.72x |
| taxol | def2-svpd | 1657 | 2 | 0.01 | 0.0018 | 0.0007 | 2.62x | 0.0013 | 1.80x |
| taxol | def2-tzvpd | 2743 | 3 | 0.02 | 0.0033 | 0.0014 | 2.42x | 0.0028 | 2.05x |
| taxol | def2-qzvpd | 5505 | 4 | 0.06 | 0.0083 | 0.0027 | 3.06x | 0.0106 | 3.94x |
| crambin | def2-svpd | 9294 | 2 | 0.11 | 0.0159 | 0.0029 | 5.40x | 0.0269 | 9.14x |
| crambin | def2-tzvpd | 15180 | 3 | 0.28 | 0.0338 | 0.0058 | 5.86x | 0.0733 | 12.72x |
| crambin | def2-qzvpd | 31284 | 4 | 0.84 | 0.0903 | 0.0134 | 6.72x | 0.3115 | 23.17x |
| ubiquitin | def2-svpd | 17433 | 2 | 0.26 | 0.0374 | 0.0059 | 6.35x | 0.0902 | 15.28x |
| ubiquitin | def2-tzvpd | 28298 | 3 | 0.63 | 0.0779 | 0.0107 | 7.31x | 0.2313 | 21.69x |
| ubiquitin | def2-qzvpd | 59053 | 4 | 1.94 | 0.2069 | 0.0266 | 7.78x | 1.0796 | 40.61x |

### The correlation consistent sets

| molecule | basis | nao | lmax | sparse GB | compute 1 thr | compute 14 thr | scaling | reference 14 thr | against it |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | cc-pvdz | 683 | 2 | 0.00 | 0.0013 | 0.0006 | 2.38x | 0.0010 | 1.75x |
| tagrisso | cc-pvtz | 1572 | 3 | 0.01 | 0.0020 | 0.0008 | 2.38x | 0.0017 | 2.07x |
| tagrisso | cc-pvqz | 3025 | 4 | 0.02 | 0.0036 | 0.0013 | 2.84x | 0.0043 | 3.41x |
| tagrisso | cc-pv5z | 5182 | 5 | 0.04 | 0.0077 | 0.0025 | 3.04x | not valid |  |
| tagrisso | cc-pv6z | 8183 | 6 | 0.10 | 0.0180 | 0.0057 | 3.15x | not valid |  |
| tagrisso | aug-cc-pvdz | 1148 | 2 | 0.00 | 0.0019 | 0.0007 | 2.53x | 0.0013 | 1.68x |
| tagrisso | aug-cc-pvtz | 2461 | 3 | 0.02 | 0.0035 | 0.0012 | 2.88x | 0.0026 | 2.15x |
| tagrisso | aug-cc-pvqz | 4478 | 4 | 0.05 | 0.0077 | 0.0025 | 3.05x | 0.0080 | 3.19x |
| tagrisso | aug-cc-pv5z | 7339 | 5 | 0.12 | 0.0188 | 0.0058 | 3.22x | not valid |  |
| tagrisso | aug-cc-pv6z | 11184 | 6 | 0.25 | 0.0438 | 0.0141 | 3.11x | not valid |  |
| taxol | cc-pvdz | 1099 | 2 | 0.00 | 0.0016 | 0.0007 | 2.29x | 0.0013 | 1.83x |
| taxol | cc-pvtz | 2516 | 3 | 0.01 | 0.0027 | 0.0011 | 2.56x | 0.0037 | 3.49x |
| taxol | cc-pvqz | 4825 | 4 | 0.04 | 0.0056 | 0.0021 | 2.73x | 0.0107 | 5.18x |
| taxol | cc-pv5z | 8246 | 5 | 0.10 | 0.0138 | 0.0047 | 2.97x | not valid |  |
| taxol | cc-pv6z | 12999 | 6 | 0.22 | 0.0341 | 0.0107 | 3.19x | not valid |  |
| taxol | aug-cc-pvdz | 1844 | 2 | 0.01 | 0.0027 | 0.0011 | 2.58x | 0.0019 | 1.80x |
| taxol | aug-cc-pvtz | 3933 | 3 | 0.04 | 0.0058 | 0.0020 | 2.87x | 0.0052 | 2.56x |
| taxol | aug-cc-pvqz | 7134 | 4 | 0.11 | 0.0148 | 0.0046 | 3.22x | 0.0179 | 3.87x |
| taxol | aug-cc-pv5z | 11667 | 5 | 0.27 | 0.0370 | 0.0113 | 3.27x | not valid |  |
| taxol | aug-cc-pv6z | 17752 | 6 | 0.57 | 0.0905 | 0.0265 | 3.41x | not valid |  |
| crambin | cc-pvdz | 6177 | 2 | 0.03 | 0.0099 | 0.0021 | 4.76x | 0.0226 | 10.88x |
| crambin | cc-pvtz | 14244 | 3 | 0.13 | 0.0212 | 0.0035 | 6.00x | 0.0739 | 20.89x |
| crambin | cc-pvqz | 27459 | 4 | 0.38 | 0.0471 | 0.0071 | 6.64x | 0.2593 | 36.57x |
| crambin | cc-pv5z | 47106 | 5 | 1.03 | 0.1255 | 0.0167 | 7.53x | not valid |  |
| crambin | cc-pv6z | 74469 | 6 | 2.38 | 0.3090 | 0.0397 | 7.79x | not valid |  |
| crambin | aug-cc-pvdz | 10380 | 2 | 0.17 | 0.0297 | 0.0049 | 6.05x | 0.0433 | 8.83x |
| crambin | aug-cc-pvtz | 22311 | 3 | 0.62 | 0.0721 | 0.0101 | 7.16x | 0.1493 | 14.83x |
| crambin | aug-cc-pvqz | 40674 | 4 | 1.69 | 0.1838 | 0.0252 | 7.29x | 0.5456 | 21.64x |
| crambin | aug-cc-pv5z | 66753 | 5 | 3.93 | 0.4625 | 0.0569 | 8.13x | not valid |  |
| crambin | aug-cc-pv6z | 101832 | 6 | 8.05 | 1.0900 | 0.1231 | 8.85x | not valid |  |
| ubiquitin | cc-pvdz | 11577 | 2 | 0.06 | 0.0212 | 0.0034 | 6.15x | 0.0752 | 21.83x |
| ubiquitin | cc-pvtz | 26870 | 3 | 0.26 | 0.0419 | 0.0063 | 6.67x | 0.2393 | 38.06x |
| ubiquitin | cc-pvqz | 51984 | 4 | 0.80 | 0.0985 | 0.0130 | 7.56x | 0.8936 | 68.63x |
| ubiquitin | cc-pv5z | 89381 | 5 | 2.21 | 0.2606 | 0.0386 | 6.75x | not valid |  |
| ubiquitin | cc-pv6z | 141523 | 6 | 5.15 | 0.6434 | 0.0837 | 7.69x | not valid |  |
| ubiquitin | aug-cc-pvdz | 19511 | 2 | 0.44 | 0.0707 | 0.0091 | 7.80x | 0.1436 | 15.85x |
| ubiquitin | aug-cc-pvtz | 42163 | 3 | 1.56 | 0.1762 | 0.0231 | 7.63x | 0.5334 | 23.10x |
| ubiquitin | aug-cc-pvqz | 77098 | 4 | 4.16 | 0.4427 | 0.0583 | 7.59x | too large |  |
| ubiquitin | aug-cc-pv5z | 126778 | 5 | 9.61 | 1.0944 | 0.1370 | 7.99x | not valid |  |
| ubiquitin | aug-cc-pv6z | 193665 | 6 | 19.61 | 2.7294 | 0.3085 | 8.85x | not valid |  |

### Parallel scaling

The speedup of `compute` against one thread, warm.

| molecule | basis | nao | 1 thr | 2 thr | 4 thr | 6 thr | 8 thr | 10 thr | 12 thr | 14 thr |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-qzvp | 3099 | 1.00x | 1.72x | 2.48x | 2.81x | 2.86x | 2.96x | 2.93x | 2.48x |
| taxol | def2-qzvp | 4947 | 1.00x | 1.85x | 2.89x | 3.03x | 3.02x | 3.02x | 2.92x | 2.94x |
| crambin | def2-svp | 6177 | 1.00x | 1.76x | 2.91x | 3.79x | 4.25x | 4.53x | 4.80x | 4.81x |
| crambin | def2-qzvp | 28167 | 1.00x | 1.83x | 3.06x | 4.12x | 4.95x | 5.42x | 6.30x | 6.61x |
| ubiquitin | def2-svp | 11577 | 1.00x | 1.88x | 3.06x | 4.07x | 5.23x | 6.01x | 6.18x | 6.16x |
| ubiquitin | def2-tzvp | 22442 | 1.00x | 1.81x | 3.07x | 4.27x | 4.93x | 5.14x | 6.41x | 6.33x |
| ubiquitin | def2-qzvp | 53197 | 1.00x | 1.84x | 3.17x | 4.42x | 5.56x | 6.16x | 7.09x | 7.15x |
| crambin | cc-pv5z | 47106 | 1.00x | 1.85x | 3.39x | 4.62x | 5.01x | 5.63x | 6.83x | 7.49x |
| ubiquitin | aug-cc-pvtz | 42163 | 1.00x | 1.88x | 3.23x | 4.58x | 5.64x | 6.00x | 7.61x | 7.94x |

### What these numbers say

The scaling reaches **8.85 times** on crambin in aug-cc-pv6z and the same 8.85 on
ubiquitin in it, which is the best measured on either driver. It is better than
the overlap, 7.12 times on the same case, because an atom pair of the kinetic
energy carries more arithmetic and the cost a block pays whatever it holds is a
smaller share of it.

The scaling is near linear to two threads, 1.7 to 1.9 times, and then divides by
the size of the molecule: tagrisso reaches 2.9 times and stops, holding too few
atom pairs to be divided into blocks at all, while ubiquitin climbs to fourteen.
The knee visible at eight to ten threads on several rows is the machine and not
the code, its ten performance cores being filled before the four efficiency cores
are reached.

Against the reference where it is valid, the driver is between 8.9 and 58.9 times
faster on the four large def2 cases and between 1.1 and 4.1 on the two small ones,
and it wins every case of every table. Tagrisso in def2-svp, which was the one
loss of the previous run at 0.93 times, is 1.07 times here; it is half a
millisecond of work and it sits on the parity line either way.

The largest case is ubiquitin in aug-cc-pv6z, 193665 basis functions and 19.6
gigabytes sparse, whose kinetic energy matrix is formed in 0.31 seconds on
fourteen threads.

## The two-center Coulomb driver

The two-center Coulomb integrals are not screenable. The operator decays as the
inverse of the interatomic distance, so no atom pair and no pair of primitives
falls below a threshold at any separation a molecule reaches, and the matrix is
dense for every molecule. Three things follow, and they are what this driver does
differently from the overlap and the kinetic energy.

The matrix is a `CPackedMatrix` rather than a `CSparseMatrix`: the lower triangle
of a symmetric matrix, in one allocation, with no sparsity pattern to describe.
For crambin in the jkfit set that is 3.52 gigabytes where the full square would
be 7.05, and for ubiquitin in it 12.09 against 24.18, which is the difference
between fitting in this machine and not.

The blocks divide the work alone and no longer divide the storage. Each of them
computes one combination of basis functions at a time into a buffer of its own
and adds it to the matrix, which the blocks share and write disjoint elements of.

The exponential of a pair of primitives is replaced by the Boys function, which
`simdfunc::compute_boys_function` evaluates for all pairs of primitives of a
kernel in one call. That routine is scalar and branchy, unlike the vector
exponential the other two drivers lean on.

### The def2 universal fitting sets

Warm, seconds, best of three after one cold call, each case in its own process.
Both drivers are timed the same way. The reference reaches angular momentum six
and both fitting sets stop at four, so it is valid everywhere here, ubiquitin
included: its dense matrix for jkfit peaks at 15.9 gigabytes resident, not the
24.18 the square would suggest, and it runs without swapping.

| molecule | basis | nao | lmax | packed GB | compute 1 thr | compute 14 thr | scaling | reference 14 thr | against it |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | jfit | 2176 | 4 | 0.02 | 0.0063 | 0.0013 | 4.83x | 0.0024 | 1.83x |
| tagrisso | jkfit | 3387 | 4 | 0.04 | 0.0137 | 0.0025 | 5.53x | 0.0040 | 1.61x |
| taxol | jfit | 3528 | 4 | 0.05 | 0.0168 | 0.0027 | 6.15x | 0.0048 | 1.76x |
| taxol | jkfit | 5489 | 4 | 0.11 | 0.0363 | 0.0060 | 6.06x | 0.0088 | 1.47x |
| crambin | jfit | 19500 | 4 | 1.42 | 0.5864 | 0.0717 | 8.18x | 0.1234 | 1.72x |
| crambin | jkfit | 30751 | 4 | 3.52 | 1.4821 | 0.1572 | 9.43x | 0.2543 | 1.62x |
| ubiquitin | jfit | 36419 | 4 | 4.94 | 2.2069 | 0.3166 | 6.97x | 0.4098 | 1.29x |
| ubiquitin | jkfit | 56971 | 4 | 12.09 | 5.5086 | 0.7383 | 7.46x | 0.8483 | 1.15x |

The values were checked against the reference on tagrisso and taxol in both sets
while the timings were taken: the largest deviation is 6.1e-13 on elements
reaching 92.7, which is 6.6e-15 relative.

### The allocation was a third of ubiquitin

The first version of this driver scaled to 5.8 times on crambin in jfit but only
3.7 on ubiquitin, and lost more ground to the reference the larger the molecule
got. Neither was the integrals. `CPackedMatrix` allocated its values with
`_values.resize(n, 0.0)`, a single serial zero fill of the whole triangle, and
`compute` paid it on every call: a tenth of the call up to crambin and **a third
of it on ubiquitin**, at 23 gigabytes per second, which is one core writing at the
bandwidth of the memory.

The values are now allocated without being value initialized, through a
default-init allocator, and zeroed by the threads with the chunked `zero()` the
class already had. The fill cannot be dropped altogether: the off-diagonal blocks
write every element they own, but the diagonal atom blocks do not, as the
one-center integral is diagonal in the angular components and the elements of two
differing momenta or components are left at the zero the matrix was constructed
with.

Seconds, fourteen threads, the construction of the packed matrix alone.

| molecule | basis | packed GB | serial fill | threaded fill | gain |
|---|---|---|---|---|---|
| taxol | jkfit | 0.11 | 0.0009 | 0.0003 | 2.95x |
| crambin | jfit | 1.42 | 0.0108 | 0.0038 | 2.82x |
| crambin | jkfit | 3.52 | 0.0269 | 0.0094 | 2.87x |
| ubiquitin | jfit | 4.94 | 0.2149 | 0.0734 | 2.93x |
| ubiquitin | jkfit | 12.09 | 0.5210 | 0.1862 | 2.80x |

The fill is the one part of this driver which is purely bandwidth bound, and it
shows: it saturates at 1.48 times on crambin and 2.84 on ubiquitin whatever the
threads, where the integrals reach seven to ten. That contrast is what settles
whether anything else here is memory bound, and nothing else is.

### The blocks were one to two orders of magnitude too large

`block_size` was `max(npairs / (blocks_per_thread * nthreads), min_block_size)`
with `blocks_per_thread` at the two inherited from the sparse matrix. On fourteen
threads that is 7348 atom pairs a block for crambin and 27038 for ubiquitin. The
right size is between one and three hundred.

A block of the screened path carries the bisection of the screening and is worth
making large. A block here is cheap to start, so a large molecule is better
divided into many small blocks: the dynamic loop then has enough of them to even
out the ones which differ in cost, and each of them keeps its Boys function and
its buffer in the cache. The size is now capped by `max_block_size`, 256, which
is a ceiling on the block rather than a count of blocks and so does not move with
the threads.

Milliseconds, fourteen threads, best of three, against the block size the setting
produces.

| block size | crambin jfit | crambin jkfit | ubiquitin jkfit |
|---|---|---|---|
| 7348 / 27038, the old value | 96.15 | 271.30 | 1252.54 |
| 918 / 3379 | 75.62 | 194.64 | 1003.15 |
| 306 / 1126 | 72.85 | 158.97 | 875.34 |
| 153 / 563 | 74.02 | 158.04 | 797.35 |
| 128 / 422 | 74.55 | 161.57 | 776.70 |
| 128 / 211 | 75.66 | 160.25 | 730.84 |
| 128 / 128 | 74.50 | 160.77 | 718.92 |

The time is flat below about three hundred and rises steeply above it, so the
exact value matters little and the ceiling is taken inside the flat range rather
than at the floor, to keep a molecule between the two from being divided more
finely than it repays. The sweep ran from one to five hundred and twelve blocks
per thread, rebuilding for each, and the table is the resulting block size because
that, not the count, is what the timings follow.

The parallel scaling is what moved. Speedup of `compute` against one thread.

| case | 1 thr | 2 thr | 4 thr | 6 thr | 8 thr | 10 thr | 12 thr | 14 thr |
|---|---|---|---|---|---|---|---|---|
| crambin jfit, before | 1.00x | 1.69x | 3.00x | 3.17x | 4.72x | 4.82x | 5.72x | 5.96x |
| crambin jfit, after | 1.00x | 1.97x | 3.60x | 5.16x | 6.65x | 7.93x | 8.11x | 8.23x |
| crambin jkfit, before | 1.00x | 1.57x | 2.87x | 3.19x | 4.54x | 4.60x | 5.18x | 5.39x |
| crambin jkfit, after | 1.00x | 2.23x | 4.15x | 5.92x | 7.49x | 8.88x | 9.17x | 9.83x |
| ubiquitin jkfit, before | 1.00x | 1.55x | 2.61x | 3.05x | 3.41x | 4.05x | 4.28x | 4.31x |
| ubiquitin jkfit, after | 1.00x | 2.04x | 3.76x | 5.03x | 6.24x | 7.24x | 7.37x | 7.39x |

Crambin in jkfit at **9.83 times** is the best scaling in these notes, ahead of
the 8.85 of the kinetic energy. The plateau from four to six threads which the
old setting showed on every one of the three cases, and which the kinetic energy
did not show on the same machine, was the block decomposition and is gone.

The single threaded column is untouched by any of this: `make_block_size` returns
zero below two threads and the groups are then left undivided, so no setting of
either bound can reach it.

### Where the time goes, and what the profile got wrong

Time Profiler on crambin in jkfit, restricted to the samples inside the driver.

| phase | one thread | fourteen threads |
|---|---|---|
| the kernels, about twenty of them | 47.7% | 61.0% |
| `compute_boys_function` itself | 19.7% | 16.3% |
| `exp` from libm, called by it | 19.4% | 10.1% |
| memmove and bzero | 5.3% | 3.1% |
| the scatter into the packed matrix | 1.2% | 0.6% |
| the coordinates and the solid harmonics | 0.2% | 0.1% |

The kernel share has no hotspot, the largest single kernel being 3.5 per cent and
the rest tailing off evenly, so anything done there is a change to the generator
and not to one file.

**The sampling overstates `exp` by about a factor of two, and acting on it would
have been a mistake.** Building with the call replaced by a cheap expression
measures its cost directly, and it is 11 per cent of the driver on one thread and
**2 per cent on fourteen**: crambin in jkfit goes from 1.4688 to 1.3076 seconds on
one thread and from 0.2665 to 0.2600 on fourteen, and ubiquitin in jkfit shows no
gain at fourteen at all. The Boys function is the second largest symbol in the
profile and very nearly free at the thread count the driver is used with.

Two things were tried against it and rejected on measurement:

The exponential is not evaluated where it underflows, `exp(-745.2)` being already
zero in double precision, which is 27 per cent of the arguments on crambin and 35
on ubiquitin. It is worth 0.18 nanoseconds of a 4.67 nanosecond point, because
libm already returns quickly for a large negative argument, so about 0.7 per cent
on crambin. It was implemented, measured and reverted rather than kept for the
extra branch.

`CSimdVariableMatrix` is sharply sensitive to the stride between the blocks a row
is spread over: at order eight a micro-benchmark runs at 26 million points a
second with 4096 columns and 81 million with 4032, a threefold swing from cache
set conflicts between the order and one streams. The driver does not reach it. The
byte stride is `nprims * pitch_of(npairs) * 8` and lands on a multiple of 4096 for
0.2 per cent of the kernel calls of crambin and 4 per cent of ubiquitin, the large
blocks giving `pitch_of(7348) = 2^6 * 919` and `pitch_of(27038) = 2^8 * 845`. At
the real block sizes the Boys function runs at 7.10 to 7.69 nanoseconds a point
against 7.37 for a small in cache case, and 9.59 to 9.91 for the power of two
geometries. No padding was added.

The cost of the Boys function itself is 2.65 nanoseconds fixed and 0.61 for each
order, in cache, so the recursion is about two thirds of it at order eight and the
exponential and the table the rest.

### The block floor of the screened path is wrong here

The driver started with the `min_block_size` of the sparse matrix, two thousand
and forty eight. That floor exists because a block of the screened path carries
the bisection of the screening as its fixed cost, and a block too small pays it
for too little work. There is no bisection here, and a Coulomb block does much
more arithmetic per atom pair, so the fixed cost is repaid by far fewer of them.

The floor was costing about half the throughput on anything below a few hundred
atoms. Tagrisso holds 2415 atom pairs in ten groups and its largest group holds
924, so at a floor of two thousand and forty eight nothing was ever divided and
the driver fell back to parallelism over the groups, whose ceiling is the share
of that largest group. It is 38 per cent, which bounds the speedup at 2.6 times,
and 2.59 is what it measured.

Milliseconds, best of three, on fourteen threads.

| case | 2048 | 1024 | 512 | 256 | 128 | 64 | 32 | 16 |
|---|---|---|---|---|---|---|---|---|
| tagrisso jfit | 2.7 | 2.8 | 2.7 | 1.9 | 1.4 | 1.5 | 1.6 | 1.5 |
| tagrisso jkfit | 6.2 | 5.9 | 5.9 | 3.6 | 2.8 | 2.9 | 2.9 | 3.0 |
| taxol jfit | 6.6 | 6.7 | 4.0 | 3.3 | 3.3 | 3.2 | 3.2 | 3.2 |
| taxol jkfit | 15.5 | 15.9 | 8.7 | 6.8 | 6.9 | 6.9 | 6.5 | 7.2 |

Everything from two hundred and fifty six down is within a few per cent of
everything else, and the whole of the gain is in leaving two thousand behind. The
number of blocks tells the same story: tagrisso goes from ten blocks to
twenty five at a floor of one hundred and twenty eight, and no further, because
below eighty six the size computed from the threads takes over and the floor
stops binding.

### Sixty four against one hundred and twenty eight

The two candidates needed more than one run each to separate. Medians of seven
process repeats, each the best of three, in milliseconds. The fragments are the
first thirty and forty five atoms of taxol, sizes at which the floor still binds.

| case | 128 | 64 | 64 against 128 |
|---|---|---|---|
| 30 atoms jkfit | 2.10 | 1.60 | 1.31x |
| 45 atoms jkfit | 2.70 | 3.00 | 0.90x |
| tagrisso jfit | 1.40 | 1.50 | 0.93x |
| tagrisso jkfit | 2.90 | 3.00 | 0.97x |
| taxol jfit | 3.20 | 3.30 | 0.97x |
| taxol jkfit | 6.80 | 6.80 | 1.00x |

There is a real crossover near forty atoms rather than a single best value. Sixty
four is a third faster below it and three to ten per cent slower above it, and
the driver uses one hundred and twenty eight: the molecules sixty four wins are
already under two milliseconds, and the ones it loses are the ones whose cost is
worth anything. A first pass with three repeats had chosen sixty four, which four
of these six cases do not support.

The floor binds for the small molecules alone. Crambin asked for blocks of 7348
atom pairs of its own accord at the time, so no floor below that could reach it,
and measuring it confirmed as much: 0.1096 seconds at two thousand and forty eight
against 0.1079 at sixty four, thirty seven blocks either way.

Both of these subsections describe the driver before the ceiling of the section
above existed, and the sentence they end on is the reason it does. A large
molecule was asking for blocks two orders of magnitude larger than it wanted, and
the floor could say nothing about that because a floor only ever raises. The
measurements themselves stand: the floor is what decides the small molecules and
the ceiling what decides the large ones, and between them lies a band of a factor
of two. On fourteen threads tagrisso is held at the floor of 128, taxol asks for
214 of its own accord and gets it, and crambin and ubiquitin ask for 7348 and
27038 and are held at the ceiling of 256.

### A trap in timing the reference

The reference driver returns a dense matrix, which for crambin in the jkfit set
is 7.05 gigabytes. A single call therefore pays for faulting in every page of it,
and timing one call rather than the best of several overstates it by a factor of
three: 0.0111 seconds against 0.0024 on tagrisso in jfit, 0.5869 against 0.2635
on crambin in jkfit. An early comparison was against cold calls and reported the
driver as 1.2 to 1.9 times faster than the reference. It is not. The table above
times both warm. The four numbers in this paragraph are from the run which found
the trap and were not measured again.

### What these numbers say

The driver beats the reference everywhere, from **1.15 times** on ubiquitin in
jkfit to **1.83** on tagrisso in jfit. It did not at first: crambin in jkfit was
12 per cent slower and ubiquitin in jkfit 48 per cent, and the loss grew
monotonically with the size of the matrix. Two changes closed it, and the largest
case more than halved, 1.6366 seconds to 0.7383.

| molecule | basis | first version | now | gain | against the reference, then and now |
|---|---|---|---|---|---|
| tagrisso | jfit | 0.0014 | 0.0013 | 1.04x | 1.77x → 1.83x |
| tagrisso | jkfit | 0.0028 | 0.0025 | 1.13x | 1.43x → 1.61x |
| taxol | jfit | 0.0031 | 0.0027 | 1.14x | 1.54x → 1.76x |
| taxol | jkfit | 0.0067 | 0.0060 | 1.12x | 1.31x → 1.47x |
| crambin | jfit | 0.0999 | 0.0717 | 1.39x | 1.24x → 1.72x |
| crambin | jkfit | 0.2900 | 0.1572 | 1.85x | 0.88x → 1.62x |
| ubiquitin | jfit | 0.5995 | 0.3166 | 1.89x | 0.68x → 1.29x |
| ubiquitin | jkfit | 1.6366 | 0.7383 | 2.22x | 0.52x → 1.15x |

**Neither change was in the integrals.** One was the zero fill of the packed
matrix, which is memory and not arithmetic, and the other the size of a block,
which is scheduling. The kernels are 61 per cent of the driver on fourteen
threads and were not touched. The line of attack this session first proposed, the
exponential of the Boys function, was measured at 2 per cent of the driver at
that thread count and abandoned. On crambin in jkfit the two changes which were
made save 0.1328 seconds where removing the exponential outright saves 0.0065,
twenty times as much, and on ubiquitin the exponential is worth nothing
measurable at all.

The margin over the reference is still nothing like the 30 to 75 times the
overlap and the kinetic energy reach, and that is expected rather than a defect.
Almost all of their win comes from screening the work away, and there is nothing
to screen here. What is left is the efficiency of the kernels alone, and on that
footing 1.15 to 1.83 times is the honest measure of them.

Angular momentum still moves the driver against the reference, though it no longer
decides a win from a loss: jkfit is worse than jfit for the same molecule at every
size, and jkfit is the set carrying more of the high momenta. The kernels of two
non-zero angular momenta are the closed form expansions and grow quickly with the
momenta, (h|J|i) alone holding 2909 terms, where the reference walks a recursion.
Size moves it too, and what is left of that after the threaded fill has not been
run down.

What is left to try is the kernels themselves, which means the generator: they are
61 per cent of the time, spread evenly over about twenty of them with no hotspot,
and they are not memory bound, since the integrals scale to seven and ten times
where the zero fill of the same data saturates at 1.5 and 2.8. The ceiling and the
floor on the block size meet within a factor of two of each other and the timings
are flat between them, so there is nothing more to win there.
