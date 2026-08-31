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
| Threads | stated per table, `OMP_NUM_THREADS` set explicitly where it matters |

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
