# Benchmarking notes

Notes on the performance of the sparse matrix and SIMD integrals path. Each
entry records what was measured, on what, and what the numbers do and do not
cover.

The tables which describe the three drivers as they stand — the def2 tables of the
overlap and the kinetic energy and the fitting set table of the two-center Coulomb
— were measured again on the build which carries the opt-in reuse of the freed
blocks of values, so they can be read against each other. The diffuse and
correlation consistent tables of the overlap and the kinetic energy and the
parallel scaling of the kinetic energy were measured on the build before it; the
reuse is off for these three drivers, so those numbers still describe them, but
they were taken on a different day and the machine drifts by a few per cent
between runs. The sections which record an intermediate state of the code, the
kernel profile, the sweeps of the blocks and the block floor, the Instruments
findings and the dense reconstruction, keep the numbers of the run which produced
them and were not repeated; they say so where it matters.

Every overlap and kinetic energy table before the two sections on the generated
kernels describes the hand written kernels, which were replaced by a generated set,
and was measured under the block size constants which preceded the fit recorded
there. They are kept as the record of what was measured at the time and do not
describe the drivers as they stand. The sections named `The generated overlap
kernels`, `The generated kinetic energy kernels` and `The generated two-center
Coulomb kernels` do.

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
| c60 | 60 |
| tagrisso | 70 |
| taxol | 110 |
| Cu(PPh3)4 cation | 137 |
| paracetamol cluster | 320 |
| crambin | 642 |
| ubiquitin | 1231 |

The bases are def2-svp, def2-tzvp and def2-qzvp throughout, their diffuse
counterparts def2-svpd, def2-tzvpd and def2-qzvpd in the section on them, and the
correlation consistent sets to sextuple zeta in the section at the end. The section
on the generated overlap kernels adds c60 and the paracetamol cluster, the ten def2
sets, and the correlation consistent sets with and without diffuse functions.

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
| tagrisso | def2-svp | 683 | 0.0008 | 0.0003 | 0.0011 | 0.0005 | 2.54x | 2.12x |
| tagrisso | def2-tzvp | 1345 | 0.0009 | 0.0003 | 0.0017 | 0.0007 | 2.73x | 2.30x |
| tagrisso | def2-qzvp | 3099 | 0.0011 | 0.0003 | 0.0033 | 0.0011 | 3.20x | 2.96x |
| taxol | def2-svp | 1099 | 0.0008 | 0.0004 | 0.0013 | 0.0005 | 2.29x | 2.40x |
| taxol | def2-tzvp | 2185 | 0.0009 | 0.0003 | 0.0022 | 0.0009 | 2.60x | 2.44x |
| taxol | def2-qzvp | 4947 | 0.0011 | 0.0004 | 0.0052 | 0.0017 | 2.85x | 3.04x |
| crambin | def2-svp | 6177 | 0.0035 | 0.0009 | 0.0077 | 0.0016 | 3.93x | 4.93x |
| crambin | def2-tzvp | 12063 | 0.0037 | 0.0009 | 0.0152 | 0.0030 | 4.14x | 5.05x |
| crambin | def2-qzvp | 28167 | 0.0043 | 0.0011 | 0.0412 | 0.0065 | 3.93x | 6.34x |
| ubiquitin | def2-svp | 11577 | 0.0097 | 0.0016 | 0.0176 | 0.0027 | 5.88x | 6.49x |
| ubiquitin | def2-tzvp | 22442 | 0.0100 | 0.0017 | 0.0334 | 0.0053 | 5.84x | 6.34x |
| ubiquitin | def2-qzvp | 53197 | 0.0104 | 0.0018 | 0.0875 | 0.0131 | 5.78x | 6.68x |

The scaling reaches 6.3 to 6.7 times on the three ubiquitin cases and rises with
the size of the problem throughout.

### Against the reference driver

| molecule | basis | driver 1 thr | reference 1 thr | 1 thr | driver 14 thr | reference 14 thr | 14 thr |
| --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-svp | 0.0011 | 0.0012 | 1.09x | 0.0005 | 0.0004 | 0.72x |
| tagrisso | def2-tzvp | 0.0017 | 0.0042 | 2.46x | 0.0007 | 0.0012 | 1.62x |
| tagrisso | def2-qzvp | 0.0033 | 0.0203 | 6.12x | 0.0011 | 0.0061 | 5.42x |
| taxol | def2-svp | 0.0013 | 0.0029 | 2.18x | 0.0005 | 0.0006 | 1.16x |
| taxol | def2-tzvp | 0.0022 | 0.0096 | 4.31x | 0.0009 | 0.0017 | 1.84x |
| taxol | def2-qzvp | 0.0052 | 0.0513 | 9.88x | 0.0017 | 0.0066 | 3.86x |
| crambin | def2-svp | 0.0077 | 0.0990 | 12.85x | 0.0016 | 0.0125 | 8.00x |
| crambin | def2-tzvp | 0.0152 | 0.3380 | 22.19x | 0.0030 | 0.0422 | 14.00x |
| crambin | def2-qzvp | 0.0412 | 1.8104 | 43.93x | 0.0065 | 0.2066 | 31.78x |
| ubiquitin | def2-svp | 0.0176 | 0.3543 | 20.10x | 0.0027 | 0.0426 | 15.68x |
| ubiquitin | def2-tzvp | 0.0334 | 1.1749 | 35.16x | 0.0053 | 0.1363 | 25.86x |
| ubiquitin | def2-qzvp | 0.0875 | 6.4360 | 73.56x | 0.0131 | 0.6944 | 53.02x |

Every case beats the reference on a single thread and all but one of them on
fourteen. Tagrisso in def2-svp is the exception at 0.72 times, which is half a
millisecond of work against four tenths of one, and it has moved either side of
parity between runs: earlier measurements put it at 0.84 and at 1.06.

The margin runs to 74 times on a single thread and 53 on fourteen, both of them
ubiquitin in def2-qzvp. That case had no reference at all in an earlier version of
this table, on the assumption that its 21 gigabytes as a dense matrix would not
fit. The reference driver does not store the full square, and the case in fact
runs in 13.5 gigabytes of resident memory without swapping, so it is measured
here: 6.4360 seconds on one thread and 0.6944 on fourteen, against 0.0875 and
0.0131 for 0.87 gigabytes sparse.

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
| tagrisso | def2-svp | 683 | 2 | 0.00 | 0.0012 | 0.0005 | 2.24x | 0.0005 | 0.91x |
| tagrisso | def2-tzvp | 1345 | 3 | 0.00 | 0.0018 | 0.0008 | 2.25x | 0.0010 | 1.28x |
| tagrisso | def2-qzvp | 3099 | 4 | 0.02 | 0.0040 | 0.0013 | 2.96x | 0.0041 | 3.07x |
| taxol | def2-svp | 1099 | 2 | 0.00 | 0.0014 | 0.0006 | 2.20x | 0.0008 | 1.26x |
| taxol | def2-tzvp | 2185 | 3 | 0.01 | 0.0025 | 0.0010 | 2.39x | 0.0022 | 2.15x |
| taxol | def2-qzvp | 4947 | 4 | 0.04 | 0.0065 | 0.0022 | 2.92x | 0.0087 | 3.89x |
| crambin | def2-svp | 6177 | 2 | 0.03 | 0.0082 | 0.0016 | 5.04x | 0.0154 | 9.42x |
| crambin | def2-tzvp | 12063 | 3 | 0.10 | 0.0179 | 0.0034 | 5.27x | 0.0503 | 14.79x |
| crambin | def2-qzvp | 28167 | 4 | 0.44 | 0.0549 | 0.0082 | 6.68x | 0.2666 | 32.45x |
| ubiquitin | def2-svp | 11577 | 2 | 0.06 | 0.0191 | 0.0029 | 6.65x | 0.0489 | 16.98x |
| ubiquitin | def2-tzvp | 22442 | 3 | 0.22 | 0.0381 | 0.0058 | 6.59x | 0.1676 | 29.00x |
| ubiquitin | def2-qzvp | 53197 | 4 | 0.95 | 0.1172 | 0.0150 | 7.79x | 0.9166 | 60.95x |

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

Against the reference where it is valid, the driver is between 9.4 and 61.0 times
faster on the six crambin and ubiquitin def2 cases and between 1.26 and 3.89 on
the five larger tagrisso and taxol ones. The single loss of the def2 table is tagrisso
in def2-svp at 0.91 times, which has been measured at 0.93, at 1.07 and at 0.91 on
three runs of builds which do not differ in this driver: it is half a millisecond
of work against four tenths of one and it sits on the parity line.

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
| tagrisso | jfit | 2176 | 4 | 0.02 | 0.0069 | 0.0014 | 4.79x | 0.0024 | 1.70x |
| tagrisso | jkfit | 3387 | 4 | 0.04 | 0.0147 | 0.0026 | 5.71x | 0.0042 | 1.62x |
| taxol | jfit | 3528 | 4 | 0.05 | 0.0179 | 0.0029 | 6.21x | 0.0049 | 1.70x |
| taxol | jkfit | 5489 | 4 | 0.11 | 0.0391 | 0.0061 | 6.43x | 0.0088 | 1.45x |
| crambin | jfit | 19500 | 4 | 1.42 | 0.6163 | 0.0729 | 8.45x | 0.1355 | 1.86x |
| crambin | jkfit | 30751 | 4 | 3.52 | 1.5387 | 0.1603 | 9.60x | 0.2702 | 1.69x |
| ubiquitin | jfit | 36419 | 4 | 4.94 | 2.2433 | 0.3120 | 7.19x | 0.4160 | 1.33x |
| ubiquitin | jkfit | 56971 | 4 | 12.09 | 5.5422 | 0.7321 | 7.57x | 0.8610 | 1.18x |

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

The driver beats the reference everywhere, from **1.18 times** on ubiquitin in
jkfit to **1.86** on crambin in jfit. It did not at first: crambin in jkfit was
12 per cent slower and ubiquitin in jkfit 48 per cent, and the loss grew
monotonically with the size of the matrix. Two changes closed it, and the largest
case more than halved, 1.6366 seconds to 0.7321.

The table below compares the first version against the run which was current when
the two changes were made; its `now` column is that run and not the fitting set
table above, which was measured later.

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

## The blocks of values a thread frees, and who should reuse them

`CSimdMatrix` allocates the values of every matrix through `::operator new[]` with
a 128 byte alignment, and the drivers form those matrices inside a parallel
region. When several threads ask for the same large aligned size at the same
moment, the allocator of the system serializes them. The cost does not divide by
the threads; it grows with them.

Measured on the coordinates and the solid harmonics of one atom on c side, the
shapes the three-center driver forms, against a loop which allocates and frees the
same seven buffers and computes nothing:

| ncols | threads | coordinates and harmonics | allocation alone |
| --- | --- | --- | --- |
| 1830 | 1 | 15.758 ms | 0.348 ms |
| 1830 | 14 | 3.977 ms | **4.229 ms** |
| 10133 | 1 | 56.509 ms | 0.199 ms |
| 10133 | 14 | 8.377 ms | **3.838 ms** |

The allocation alone takes 4.2 milliseconds of wall time on fourteen threads
against 0.35 on one, which is 170 times the processor time for the same work. It
is a lock, and on the small case it costs as much as the arithmetic it serves.

### A cache of the freed blocks, and why it cannot be on for everyone

Each thread keeps the blocks it frees, keyed by size, and takes one back when a
matrix of that size is formed. Two blocks per size, sixteen sizes, sixteen
megabytes per thread, nothing below four kilobytes, and the size which has gone
unused for the longest is evicted when either bound is reached.

Turned on for every matrix this is not an improvement. The two-center drivers form
a few matrices per block of atom pairs, where the allocations are already rare
against the work of the block, and they pay for the memory the cache holds back
from the allocator:

| case | cache off | always on | opt in |
| --- | --- | --- | --- |
| overlap crambin def2-qzvp | 0.00641 | 0.00673 (+5.0%) | 0.00655 (+2.2%) |
| overlap ubiquitin def2-qzvp | 0.01261 | 0.01335 (+5.9%) | 0.01259 (−0.2%) |
| kinetic ubiquitin def2-qzvp | 0.01498 | 0.01618 (+8.0%) | 0.01506 (+0.5%) |
| coulomb crambin jkfit | 0.16227 | 0.15783 (−2.7%) | 0.15874 (−2.2%) |
| coulomb ubiquitin jfit | 0.31609 | 0.29996 (−5.1%) | 0.30737 (−2.8%) |

The two Coulomb rows are the control. The reuse is off for that driver in both the
`always on` and the `opt in` column, so the 2 to 5 per cent they move is the noise
of the measurement and nothing else. Read against them, the overlap and the
kinetic energy lose 5 to 8 per cent when the cache is on for everyone and lose
nothing when it is not.

### No floor and no budget separates the two

The obvious repair is to cache only the sizes which benefit. There are none. The
budget scales the harm and the supposed gain together, and the floor which frees
the overlap is the floor which destroys what the three-center driver gains:

| floor | overlap ubiquitin def2-qzvp | three-center pattern, ncols 1830 |
| --- | --- | --- |
| off | 0.01271 | 4.760 ms |
| 4 KB | +5.3% | **1.368 ms** |
| 64 KB | +6.5% | 1.424 ms |
| 512 KB | 0.0% | 3.915 ms |

The three-center driver's per-atom buffers are 44 to 200 kilobytes at a typical
block, which is exactly the range whose caching costs the overlap 5 to 6 per cent.
The two populations are the same sizes. What separates them is not the size of the
block but how often it is asked for, and only the caller knows that.

### The guard

`CSimdMatrix::CBlockReuse` turns the reuse on for the calling thread while it is
alive. It is off otherwise, so a driver which does not construct it allocates and
frees exactly as it did before the cache existed, and the two-center drivers are
untouched by construction rather than by measurement. The three-center driver
constructs one per block of atom pairs, inside the parallel region and before the
matrices it governs, so it spans the loop over the atoms on c side and frees what
the thread holds when the block ends.

On the pattern it serves, with the guard placed as the driver places it, fourteen
threads: **3.774 ms to 2.099 ms, 1.80 times**. An earlier figure of 2.9 times for
the same case came from a cache which also persisted through the warm up of the
measurement and is an upper bound rather than what the driver sees.

## The generated overlap kernels

The overlap kernels were replaced by a generated set covering every combination
of angular momenta to l = 6 in both orders, built from a vertical recurrence, a
contraction and a horizontal transfer rather than from the closed form of the
earlier hand written kernels. Everything below was measured on that build, on an
otherwise idle machine, with `OMP_NUM_THREADS=14` unless a table says otherwise.
Timings are the best of two to five runs, so the scheduler and the first touch of
the pages do not count.

Three columns appear throughout. **ref** is `OverlapDriver`, which computes every
atom pair and carries no threshold, so it is timed once per molecule and basis and
the same number serves both threshold rows. **auto** is the SIMD driver at the
block size the heuristic chose before this section was written,
`blocks_per_thread = 2` and `min_block_size = 2048`. **fit** is the same driver at the block size the
constants which now ship choose, `blocks_per_thread = 4` and
`min_block_size = 256`. The benchmark predates that change and reached those sizes
through the driver's `block_size` argument, so the **fit** column is what the
driver does today and the **auto** column is what it did before.

`dense too big` marks the cases where the reference cannot be run at all: it
returns a dense matrix, and ubiquitin in aug-cc-pV6Z would be 300 GB. Those rows
carry the two SIMD timings and no ratio rather than an invented one.

### tagrisso

| basis | threshold | nao | ref ms | auto ms | fit ms | x auto | x fit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 683 | 0.85 | 0.60 | 0.78 | 1.4 | 1.1 |
| def2-svp | 1e-12 | 683 | 0.85 | 0.64 | 0.77 | 1.3 | 1.1 |
| def2-svpd | 1e-14 | 1010 | 1.76 | 0.83 | 0.93 | 2.1 | 1.9 |
| def2-svpd | 1e-12 | 1010 | 1.76 | 0.79 | 0.92 | 2.2 | 1.9 |
| def2-tzvp | 1e-14 | 1345 | 1.72 | 1.06 | 1.08 | 1.6 | 1.6 |
| def2-tzvp | 1e-12 | 1345 | 1.72 | 1.02 | 1.05 | 1.7 | 1.6 |
| def2-tzvpp | 1e-14 | 1609 | 1.18 | 0.94 | 0.95 | 1.3 | 1.2 |
| def2-tzvpp | 1e-12 | 1609 | 1.18 | 1.03 | 1.10 | 1.1 | 1.1 |
| def2-tzvpd | 1e-14 | 1672 | 1.77 | 1.37 | 1.13 | 1.3 | 1.6 |
| def2-tzvpd | 1e-12 | 1672 | 1.77 | 1.32 | 1.27 | 1.3 | 1.4 |
| def2-tzvppd | 1e-14 | 1936 | 2.43 | 1.26 | 1.17 | 1.9 | 2.1 |
| def2-tzvppd | 1e-12 | 1936 | 2.43 | 1.19 | 1.33 | 2.0 | 1.8 |
| def2-qzvp | 1e-14 | 3099 | 2.98 | 2.14 | 1.87 | 1.4 | 1.6 |
| def2-qzvp | 1e-12 | 3099 | 2.98 | 2.08 | 2.09 | 1.4 | 1.4 |
| def2-qzvpp | 1e-14 | 3099 | 4.92 | 2.12 | 1.87 | 2.3 | 2.6 |
| def2-qzvpp | 1e-12 | 3099 | 4.92 | 2.09 | 2.20 | 2.4 | 2.2 |
| def2-qzvpd | 1e-14 | 3426 | 3.51 | 2.68 | 2.17 | 1.3 | 1.6 |
| def2-qzvpd | 1e-12 | 3426 | 3.51 | 2.68 | 2.31 | 1.3 | 1.5 |
| def2-qzvppd | 1e-14 | 3426 | 3.56 | 2.72 | 2.32 | 1.3 | 1.5 |
| def2-qzvppd | 1e-12 | 3426 | 3.56 | 2.67 | 2.55 | 1.3 | 1.4 |
| cc-pvdz | 1e-14 | 683 | 0.86 | 0.61 | 0.66 | 1.4 | 1.3 |
| cc-pvdz | 1e-12 | 683 | 0.86 | 0.69 | 0.77 | 1.2 | 1.1 |
| cc-pvtz | 1e-14 | 1572 | 1.33 | 0.97 | 0.93 | 1.4 | 1.4 |
| cc-pvtz | 1e-12 | 1572 | 1.33 | 1.05 | 1.06 | 1.3 | 1.3 |
| cc-pvqz | 1e-14 | 3025 | 4.06 | 2.05 | 1.72 | 2.0 | 2.4 |
| cc-pvqz | 1e-12 | 3025 | 4.06 | 1.94 | 1.97 | 2.1 | 2.1 |
| cc-pv5z | 1e-14 | 5182 | 8.85 | 5.57 | 4.13 | 1.6 | 2.1 |
| cc-pv5z | 1e-12 | 5182 | 8.85 | 5.15 | 3.91 | 1.7 | 2.3 |
| cc-pv6z | 1e-14 | 8183 | 24.68 | 16.57 | 11.25 | 1.5 | 2.2 |
| cc-pv6z | 1e-12 | 8183 | 24.68 | 15.09 | 10.63 | 1.6 | 2.3 |
| aug-cc-pvdz | 1e-14 | 1148 | 1.04 | 0.92 | 0.86 | 1.1 | 1.2 |
| aug-cc-pvdz | 1e-12 | 1148 | 1.04 | 1.01 | 1.04 | 1.0 | 1.0 |
| aug-cc-pvtz | 1e-14 | 2461 | 2.00 | 1.87 | 1.58 | 1.1 | 1.3 |
| aug-cc-pvtz | 1e-12 | 2461 | 2.00 | 1.95 | 1.98 | 1.0 | 1.0 |
| aug-cc-pvqz | 1e-14 | 4478 | 5.74 | 4.92 | 3.74 | 1.2 | 1.5 |
| aug-cc-pvqz | 1e-12 | 4478 | 5.74 | 4.78 | 3.62 | 1.2 | 1.6 |
| aug-cc-pv5z | 1e-14 | 7339 | 16.45 | 15.28 | 10.01 | 1.1 | 1.6 |
| aug-cc-pv5z | 1e-12 | 7339 | 16.45 | 13.73 | 9.54 | 1.2 | 1.7 |
| aug-cc-pv6z | 1e-14 | 11184 | 44.94 | 45.36 | 28.02 | 1.0 | 1.6 |
| aug-cc-pv6z | 1e-12 | 11184 | 44.94 | 41.62 | 26.36 | 1.1 | 1.7 |

### c60

| basis | threshold | nao | ref ms | auto ms | fit ms | x auto | x fit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 840 | 0.44 | 0.90 | 0.56 | 0.5 | 0.8 |
| def2-svp | 1e-12 | 840 | 0.44 | 0.84 | 0.63 | 0.5 | 0.7 |
| def2-svpd | 1e-14 | 1200 | 0.62 | 1.42 | 0.85 | 0.4 | 0.7 |
| def2-svpd | 1e-12 | 1200 | 0.62 | 1.42 | 0.82 | 0.4 | 0.8 |
| def2-tzvp | 1e-14 | 1860 | 1.34 | 2.73 | 1.27 | 0.5 | 1.0 |
| def2-tzvp | 1e-12 | 1860 | 1.34 | 2.63 | 0.97 | 0.5 | 1.4 |
| def2-tzvpp | 1e-14 | 1860 | 1.59 | 2.75 | 1.26 | 0.6 | 1.3 |
| def2-tzvpp | 1e-12 | 1860 | 1.59 | 2.64 | 1.23 | 0.6 | 1.3 |
| def2-tzvpd | 1e-14 | 2220 | 2.11 | 3.89 | 1.24 | 0.5 | 1.7 |
| def2-tzvpd | 1e-12 | 2220 | 2.11 | 3.80 | 1.60 | 0.6 | 1.3 |
| def2-tzvppd | 1e-14 | 2220 | 2.02 | 3.89 | 1.70 | 0.5 | 1.2 |
| def2-tzvppd | 1e-12 | 2220 | 2.02 | 3.82 | 1.62 | 0.5 | 1.2 |
| def2-qzvp | 1e-14 | 3420 | 4.01 | 9.06 | 2.35 | 0.4 | 1.7 |
| def2-qzvp | 1e-12 | 3420 | 4.01 | 8.31 | 2.35 | 0.5 | 1.7 |
| def2-qzvpp | 1e-14 | 3420 | 3.68 | 8.97 | 2.43 | 0.4 | 1.5 |
| def2-qzvpp | 1e-12 | 3420 | 3.68 | 8.41 | 2.62 | 0.4 | 1.4 |
| def2-qzvpd | 1e-14 | 3780 | 4.56 | 11.07 | 2.94 | 0.4 | 1.6 |
| def2-qzvpd | 1e-12 | 3780 | 4.56 | 10.57 | 3.06 | 0.4 | 1.5 |
| def2-qzvppd | 1e-14 | 3780 | 3.91 | 11.21 | 2.87 | 0.3 | 1.4 |
| def2-qzvppd | 1e-12 | 3780 | 3.91 | 10.52 | 2.65 | 0.4 | 1.5 |
| cc-pvdz | 1e-14 | 840 | 1.33 | 1.23 | 0.83 | 1.1 | 1.6 |
| cc-pvdz | 1e-12 | 840 | 1.33 | 1.17 | 0.82 | 1.1 | 1.6 |
| cc-pvtz | 1e-14 | 1800 | 1.53 | 2.91 | 1.35 | 0.5 | 1.1 |
| cc-pvtz | 1e-12 | 1800 | 1.53 | 2.76 | 1.32 | 0.6 | 1.2 |
| cc-pvqz | 1e-14 | 3300 | 3.66 | 8.67 | 2.47 | 0.4 | 1.5 |
| cc-pvqz | 1e-12 | 3300 | 3.66 | 8.12 | 2.57 | 0.4 | 1.4 |
| cc-pv5z | 1e-14 | 5460 | 10.13 | 29.19 | 6.29 | 0.3 | 1.6 |
| cc-pv5z | 1e-12 | 5460 | 10.13 | 26.82 | 6.43 | 0.4 | 1.6 |
| cc-pv6z | 1e-14 | 8400 | 25.78 | 91.24 | 20.09 | 0.3 | 1.3 |
| cc-pv6z | 1e-12 | 8400 | 25.78 | 84.64 | 19.01 | 0.3 | 1.4 |
| aug-cc-pvdz | 1e-14 | 1380 | 1.22 | 2.35 | 1.17 | 0.5 | 1.0 |
| aug-cc-pvdz | 1e-12 | 1380 | 1.22 | 2.30 | 1.15 | 0.5 | 1.1 |
| aug-cc-pvtz | 1e-14 | 2760 | 2.78 | 6.83 | 1.91 | 0.4 | 1.5 |
| aug-cc-pvtz | 1e-12 | 2760 | 2.78 | 6.68 | 2.45 | 0.4 | 1.1 |
| aug-cc-pvqz | 1e-14 | 4800 | 7.20 | 23.31 | 4.97 | 0.3 | 1.4 |
| aug-cc-pvqz | 1e-12 | 4800 | 7.20 | 22.16 | 4.74 | 0.3 | 1.5 |
| aug-cc-pv5z | 1e-14 | 7620 | 19.58 | 78.10 | 15.08 | 0.3 | 1.3 |
| aug-cc-pv5z | 1e-12 | 7620 | 19.58 | 74.72 | 16.22 | 0.3 | 1.2 |
| aug-cc-pv6z | 1e-14 | 11340 | 50.11 | 245.77 | 47.81 | 0.2 | 1.0 |
| aug-cc-pv6z | 1e-12 | 11340 | 50.11 | 229.13 | 45.27 | 0.2 | 1.1 |

### taxol

| basis | threshold | nao | ref ms | auto ms | fit ms | x auto | x fit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 1099 | 0.91 | 0.64 | 0.78 | 1.4 | 1.2 |
| def2-svp | 1e-12 | 1099 | 0.91 | 0.61 | 0.77 | 1.5 | 1.2 |
| def2-svpd | 1e-14 | 1657 | 1.07 | 0.95 | 0.99 | 1.1 | 1.1 |
| def2-svpd | 1e-12 | 1657 | 1.07 | 0.96 | 0.98 | 1.1 | 1.1 |
| def2-tzvp | 1e-14 | 2185 | 2.15 | 1.34 | 1.14 | 1.6 | 1.9 |
| def2-tzvp | 1e-12 | 2185 | 2.15 | 1.29 | 1.12 | 1.7 | 1.9 |
| def2-tzvpp | 1e-14 | 2577 | 1.90 | 1.38 | 1.28 | 1.4 | 1.5 |
| def2-tzvpp | 1e-12 | 2577 | 1.90 | 1.34 | 1.25 | 1.4 | 1.5 |
| def2-tzvpd | 1e-14 | 2743 | 2.08 | 2.02 | 1.57 | 1.0 | 1.3 |
| def2-tzvpd | 1e-12 | 2743 | 2.08 | 1.95 | 1.49 | 1.1 | 1.4 |
| def2-tzvppd | 1e-14 | 3135 | 2.57 | 2.06 | 1.73 | 1.2 | 1.5 |
| def2-tzvppd | 1e-12 | 3135 | 2.57 | 1.97 | 1.74 | 1.3 | 1.5 |
| def2-qzvp | 1e-14 | 4947 | 6.50 | 3.68 | 2.82 | 1.8 | 2.3 |
| def2-qzvp | 1e-12 | 4947 | 6.50 | 3.49 | 2.70 | 1.9 | 2.4 |
| def2-qzvpp | 1e-14 | 4947 | 7.81 | 3.72 | 2.83 | 2.1 | 2.8 |
| def2-qzvpp | 1e-12 | 4947 | 7.81 | 3.35 | 2.63 | 2.3 | 3.0 |
| def2-qzvpd | 1e-14 | 5505 | 7.74 | 4.89 | 3.64 | 1.6 | 2.1 |
| def2-qzvpd | 1e-12 | 5505 | 7.74 | 4.77 | 3.49 | 1.6 | 2.2 |
| def2-qzvppd | 1e-14 | 5505 | 9.19 | 4.83 | 3.57 | 1.9 | 2.6 |
| def2-qzvppd | 1e-12 | 5505 | 9.19 | 4.74 | 3.41 | 1.9 | 2.7 |
| cc-pvdz | 1e-14 | 1099 | 2.62 | 0.78 | 0.86 | 3.4 | 3.1 |
| cc-pvdz | 1e-12 | 1099 | 2.62 | 0.75 | 0.84 | 3.5 | 3.1 |
| cc-pvtz | 1e-14 | 2516 | 3.45 | 1.42 | 1.29 | 2.4 | 2.7 |
| cc-pvtz | 1e-12 | 2516 | 3.45 | 1.37 | 1.26 | 2.5 | 2.7 |
| cc-pvqz | 1e-14 | 4825 | 9.59 | 3.53 | 2.64 | 2.7 | 3.6 |
| cc-pvqz | 1e-12 | 4825 | 9.59 | 3.26 | 2.52 | 2.9 | 3.8 |
| cc-pv5z | 1e-14 | 8246 | 23.13 | 10.50 | 6.61 | 2.2 | 3.5 |
| cc-pv5z | 1e-12 | 8246 | 23.13 | 9.82 | 6.16 | 2.4 | 3.8 |
| cc-pv6z | 1e-14 | 12999 | 56.78 | 32.49 | 16.69 | 1.7 | 3.4 |
| cc-pv6z | 1e-12 | 12999 | 56.78 | 30.31 | 15.51 | 1.9 | 3.7 |
| aug-cc-pvdz | 1e-14 | 1844 | 4.20 | 1.39 | 1.21 | 3.0 | 3.5 |
| aug-cc-pvdz | 1e-12 | 1844 | 4.20 | 1.31 | 1.18 | 3.2 | 3.5 |
| aug-cc-pvtz | 1e-14 | 3933 | 3.96 | 3.38 | 2.30 | 1.2 | 1.7 |
| aug-cc-pvtz | 1e-12 | 3933 | 3.96 | 3.27 | 2.42 | 1.2 | 1.6 |
| aug-cc-pvqz | 1e-14 | 7134 | 13.08 | 9.63 | 5.89 | 1.4 | 2.2 |
| aug-cc-pvqz | 1e-12 | 7134 | 13.08 | 9.50 | 5.93 | 1.4 | 2.2 |
| aug-cc-pv5z | 1e-14 | 11667 | 39.21 | 31.56 | 16.15 | 1.2 | 2.4 |
| aug-cc-pv5z | 1e-12 | 11667 | 39.21 | 29.36 | 15.34 | 1.3 | 2.6 |
| aug-cc-pv6z | 1e-14 | 17752 | 110.75 | 95.77 | 44.79 | 1.2 | 2.5 |
| aug-cc-pv6z | 1e-12 | 17752 | 110.75 | 90.64 | 43.71 | 1.2 | 2.5 |

### paracetamol_cluster

| basis | threshold | nao | ref ms | auto ms | fit ms | x auto | x fit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 3184 | 3.29 | 1.18 | 1.53 | 2.8 | 2.1 |
| def2-svp | 1e-12 | 3184 | 3.29 | 1.10 | 1.48 | 3.0 | 2.2 |
| def2-svpd | 1e-14 | 4768 | 6.23 | 2.06 | 2.42 | 3.0 | 2.6 |
| def2-svpd | 1e-12 | 4768 | 6.23 | 1.92 | 2.29 | 3.2 | 2.7 |
| def2-tzvp | 1e-14 | 6320 | 11.58 | 2.35 | 2.68 | 4.9 | 4.3 |
| def2-tzvp | 1e-12 | 6320 | 11.58 | 2.21 | 2.56 | 5.2 | 4.5 |
| def2-tzvpp | 1e-14 | 7472 | 14.19 | 2.55 | 3.08 | 5.6 | 4.6 |
| def2-tzvpp | 1e-12 | 7472 | 14.19 | 2.43 | 2.98 | 5.8 | 4.8 |
| def2-tzvpd | 1e-14 | 7904 | 16.13 | 4.33 | 4.31 | 3.7 | 3.7 |
| def2-tzvpd | 1e-12 | 7904 | 16.13 | 3.88 | 4.17 | 4.2 | 3.9 |
| def2-tzvppd | 1e-14 | 9056 | 22.71 | 4.61 | 5.05 | 4.9 | 4.5 |
| def2-tzvppd | 1e-12 | 9056 | 22.71 | 4.23 | 4.73 | 5.4 | 4.8 |
| def2-qzvp | 1e-14 | 14352 | 51.56 | 6.87 | 7.57 | 7.5 | 6.8 |
| def2-qzvp | 1e-12 | 14352 | 51.56 | 6.30 | 6.98 | 8.2 | 7.4 |
| def2-qzvpp | 1e-14 | 14352 | 52.83 | 6.95 | 7.54 | 7.6 | 7.0 |
| def2-qzvpp | 1e-12 | 14352 | 52.83 | 6.34 | 7.06 | 8.3 | 7.5 |
| def2-qzvpd | 1e-14 | 15936 | 68.45 | 10.62 | 11.15 | 6.4 | 6.1 |
| def2-qzvpd | 1e-12 | 15936 | 68.45 | 9.89 | 10.50 | 6.9 | 6.5 |
| def2-qzvppd | 1e-14 | 15936 | 64.11 | 10.71 | 11.34 | 6.0 | 5.7 |
| def2-qzvppd | 1e-12 | 15936 | 64.11 | 9.35 | 10.27 | 6.9 | 6.2 |
| cc-pvdz | 1e-14 | 3184 | 5.86 | 1.33 | 1.68 | 4.4 | 3.5 |
| cc-pvdz | 1e-12 | 3184 | 5.86 | 1.27 | 1.66 | 4.6 | 3.5 |
| cc-pvtz | 1e-14 | 7296 | 14.88 | 2.61 | 3.16 | 5.7 | 4.7 |
| cc-pvtz | 1e-12 | 7296 | 14.88 | 2.37 | 2.94 | 6.3 | 5.1 |
| cc-pvqz | 1e-14 | 14000 | 51.57 | 6.47 | 7.03 | 8.0 | 7.3 |
| cc-pvqz | 1e-12 | 14000 | 51.57 | 5.93 | 6.30 | 8.7 | 8.2 |
| cc-pv5z | 1e-14 | 23936 | 177.39 | 18.88 | 18.50 | 9.4 | 9.6 |
| cc-pv5z | 1e-12 | 23936 | 177.39 | 16.83 | 16.65 | 10.5 | 10.7 |
| cc-pv6z | 1e-14 | 37744 | -- | 59.41 | 52.89 | -- | -- |
| cc-pv6z | 1e-12 | 37744 | -- | 50.50 | 48.18 | -- | -- |
| aug-cc-pvdz | 1e-14 | 5344 | 10.04 | 2.91 | 3.30 | 3.4 | 3.0 |
| aug-cc-pvdz | 1e-12 | 5344 | 10.04 | 2.77 | 3.14 | 3.6 | 3.2 |
| aug-cc-pvtz | 1e-14 | 11408 | 36.85 | 7.61 | 7.88 | 4.8 | 4.7 |
| aug-cc-pvtz | 1e-12 | 11408 | 36.85 | 6.83 | 7.49 | 5.4 | 4.9 |
| aug-cc-pvqz | 1e-14 | 20704 | 111.02 | 23.63 | 21.49 | 4.7 | 5.2 |
| aug-cc-pvqz | 1e-12 | 20704 | 111.02 | 20.44 | 20.68 | 5.4 | 5.4 |
| aug-cc-pv5z | 1e-14 | 33872 | -- | 75.23 | 62.09 | -- | -- |
| aug-cc-pv5z | 1e-12 | 33872 | -- | 66.22 | 56.98 | -- | -- |
| aug-cc-pv6z | 1e-14 | 51552 | -- | 252.50 | 197.79 | -- | -- |
| aug-cc-pv6z | 1e-12 | 51552 | -- | 227.32 | 176.71 | -- | -- |

### crambin

| basis | threshold | nao | ref ms | auto ms | fit ms | x auto | x fit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 6177 | 12.71 | 1.84 | 2.11 | 6.9 | 6.0 |
| def2-svp | 1e-12 | 6177 | 12.71 | 1.76 | 2.06 | 7.2 | 6.2 |
| def2-svpd | 1e-14 | 9294 | 23.34 | 3.84 | 3.93 | 6.1 | 5.9 |
| def2-svpd | 1e-12 | 9294 | 23.34 | 3.53 | 3.60 | 6.6 | 6.5 |
| def2-tzvp | 1e-14 | 12063 | 42.07 | 4.04 | 4.32 | 10.4 | 9.7 |
| def2-tzvp | 1e-12 | 12063 | 42.07 | 3.72 | 4.00 | 11.3 | 10.5 |
| def2-tzvpp | 1e-14 | 14613 | 58.91 | 4.40 | 4.85 | 13.4 | 12.1 |
| def2-tzvpp | 1e-12 | 14613 | 58.91 | 4.01 | 4.50 | 14.7 | 13.1 |
| def2-tzvpd | 1e-14 | 15180 | 65.75 | 8.84 | 7.96 | 7.4 | 8.3 |
| def2-tzvpd | 1e-12 | 15180 | 65.75 | 6.98 | 7.69 | 9.4 | 8.6 |
| def2-tzvppd | 1e-14 | 17730 | 82.14 | 8.99 | 9.46 | 9.1 | 8.7 |
| def2-tzvppd | 1e-12 | 17730 | 82.14 | 8.01 | 8.20 | 10.3 | 10.0 |
| def2-qzvp | 1e-14 | 28167 | 218.21 | 12.76 | 13.54 | 17.1 | 16.1 |
| def2-qzvp | 1e-12 | 28167 | 218.21 | 12.01 | 11.71 | 18.2 | 18.6 |
| def2-qzvpp | 1e-14 | 28167 | 212.17 | 13.58 | 12.88 | 15.6 | 16.5 |
| def2-qzvpp | 1e-12 | 28167 | 212.17 | 11.14 | 11.57 | 19.0 | 18.3 |
| def2-qzvpd | 1e-14 | 31284 | -- | 23.00 | 22.02 | -- | -- |
| def2-qzvpd | 1e-12 | 31284 | -- | 20.61 | 19.31 | -- | -- |
| def2-qzvppd | 1e-14 | 31284 | -- | 22.58 | 23.19 | -- | -- |
| def2-qzvppd | 1e-12 | 31284 | -- | 20.46 | 20.45 | -- | -- |
| cc-pvdz | 1e-14 | 6177 | 20.01 | 2.19 | 2.41 | 9.2 | 8.3 |
| cc-pvdz | 1e-12 | 6177 | 20.01 | 2.08 | 2.31 | 9.6 | 8.7 |
| cc-pvtz | 1e-14 | 14244 | 60.89 | 4.46 | 4.83 | 13.6 | 12.6 |
| cc-pvtz | 1e-12 | 14244 | 60.89 | 3.91 | 4.39 | 15.6 | 13.9 |
| cc-pvqz | 1e-14 | 27459 | 219.01 | 11.48 | 11.69 | 19.1 | 18.7 |
| cc-pvqz | 1e-12 | 27459 | 219.01 | 10.16 | 10.63 | 21.6 | 20.6 |
| cc-pv5z | 1e-14 | 47106 | -- | 40.16 | 36.47 | -- | -- |
| cc-pv5z | 1e-12 | 47106 | -- | 33.59 | 32.42 | -- | -- |
| cc-pv6z | 1e-14 | 74469 | -- | 138.99 | 112.60 | -- | -- |
| cc-pv6z | 1e-12 | 74469 | -- | 113.51 | 101.92 | -- | -- |
| aug-cc-pvdz | 1e-14 | 10380 | 37.65 | 6.23 | 6.16 | 6.0 | 6.1 |
| aug-cc-pvdz | 1e-12 | 10380 | 37.65 | 5.25 | 5.48 | 7.2 | 6.9 |
| aug-cc-pvtz | 1e-14 | 22311 | 130.14 | 16.19 | 16.32 | 8.0 | 8.0 |
| aug-cc-pvtz | 1e-12 | 22311 | 130.14 | 14.35 | 14.64 | 9.1 | 8.9 |
| aug-cc-pvqz | 1e-14 | 40674 | -- | 58.02 | 49.83 | -- | -- |
| aug-cc-pvqz | 1e-12 | 40674 | -- | 52.04 | 45.13 | -- | -- |
| aug-cc-pv5z | 1e-14 | 66753 | -- | 214.02 | 173.04 | -- | -- |
| aug-cc-pv5z | 1e-12 | 66753 | -- | 178.37 | 152.82 | -- | -- |
| aug-cc-pv6z | 1e-14 | 101832 | -- | 664.99 | 561.23 | -- | -- |
| aug-cc-pv6z | 1e-12 | 101832 | -- | 563.69 | 488.05 | -- | -- |

### ubiquitin

| basis | threshold | nao | ref ms | auto ms | fit ms | x auto | x fit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 11577 | 42.61 | 3.61 | 3.49 | 11.8 | 12.2 |
| def2-svp | 1e-12 | 11577 | 42.61 | 3.25 | 3.46 | 13.1 | 12.3 |
| def2-svpd | 1e-14 | 17433 | 129.16 | 7.92 | 7.81 | 16.3 | 16.5 |
| def2-svpd | 1e-12 | 17433 | 129.16 | 7.06 | 6.89 | 18.3 | 18.7 |
| def2-tzvp | 1e-14 | 22442 | 145.02 | 8.37 | 7.70 | 17.3 | 18.8 |
| def2-tzvp | 1e-12 | 22442 | 145.02 | 7.50 | 6.86 | 19.3 | 21.1 |
| def2-tzvpp | 1e-14 | 27479 | 239.93 | 9.10 | 8.50 | 26.4 | 28.2 |
| def2-tzvpp | 1e-12 | 27479 | 239.93 | 7.80 | 7.66 | 30.8 | 31.3 |
| def2-tzvpd | 1e-14 | 28298 | 195.33 | 19.99 | 16.92 | 9.8 | 11.5 |
| def2-tzvpd | 1e-12 | 28298 | 195.33 | 17.60 | 14.71 | 11.1 | 13.3 |
| def2-tzvppd | 1e-14 | 33335 | -- | 22.01 | 19.65 | -- | -- |
| def2-tzvppd | 1e-12 | 33335 | -- | 19.00 | 16.78 | -- | -- |
| def2-qzvp | 1e-14 | 53197 | -- | 28.37 | 24.94 | -- | -- |
| def2-qzvp | 1e-12 | 53197 | -- | 25.11 | 22.29 | -- | -- |
| def2-qzvpp | 1e-14 | 53197 | -- | 28.22 | 25.05 | -- | -- |
| def2-qzvpp | 1e-12 | 53197 | -- | 24.57 | 21.87 | -- | -- |
| def2-qzvpd | 1e-14 | 59053 | -- | 64.42 | 52.00 | -- | -- |
| def2-qzvpd | 1e-12 | 59053 | -- | 52.09 | 43.93 | -- | -- |
| def2-qzvppd | 1e-14 | 59053 | -- | 65.25 | 52.66 | -- | -- |
| def2-qzvppd | 1e-12 | 59053 | -- | 54.70 | 45.81 | -- | -- |
| cc-pvdz | 1e-14 | 11577 | 64.63 | 3.94 | 4.41 | 16.4 | 14.6 |
| cc-pvdz | 1e-12 | 11577 | 64.63 | 3.55 | 3.81 | 18.2 | 17.0 |
| cc-pvtz | 1e-14 | 26870 | 217.75 | 8.45 | 8.33 | 25.8 | 26.1 |
| cc-pvtz | 1e-12 | 26870 | 217.75 | 7.39 | 7.55 | 29.5 | 28.9 |
| cc-pvqz | 1e-14 | 51984 | -- | 25.39 | 23.40 | -- | -- |
| cc-pvqz | 1e-12 | 51984 | -- | 22.63 | 19.37 | -- | -- |
| cc-pv5z | 1e-14 | 89381 | -- | 105.61 | 78.96 | -- | -- |
| cc-pv5z | 1e-12 | 89381 | -- | 82.92 | 67.42 | -- | -- |
| cc-pv6z | 1e-14 | 141523 | -- | 357.68 | 268.35 | -- | -- |
| cc-pv6z | 1e-12 | 141523 | -- | 289.64 | 224.50 | -- | -- |
| aug-cc-pvdz | 1e-14 | 19511 | 138.27 | 13.54 | 12.80 | 10.2 | 10.8 |
| aug-cc-pvdz | 1e-12 | 19511 | 138.27 | 12.31 | 11.74 | 11.2 | 11.8 |
| aug-cc-pvtz | 1e-14 | 42163 | -- | 44.13 | 40.11 | -- | -- |
| aug-cc-pvtz | 1e-12 | 42163 | -- | 40.76 | 34.70 | -- | -- |
| aug-cc-pvqz | 1e-14 | 77098 | -- | 192.13 | 152.97 | -- | -- |
| aug-cc-pvqz | 1e-12 | 77098 | -- | 163.96 | 126.19 | -- | -- |
| aug-cc-pv5z | 1e-14 | 126778 | -- | 620.91 | 681.71 | -- | -- |
| aug-cc-pv5z | 1e-12 | 126778 | -- | 537.40 | 444.73 | -- | -- |
| aug-cc-pv6z | 1e-14 | 193665 | -- | 2614.54 | 1826.12 | -- | -- |
| aug-cc-pv6z | 1e-12 | 193665 | -- | 1556.12 | 1577.73 | -- | -- |

### The advantage against the reference

Geometric mean and range of the ratio, over the rows which have a reference.

| molecule | rows | x ref at auto | x ref at fit |
| --- | --- | --- | --- |
| tagrisso | 40 | 1.4 (1.0 to 2.4) | 1.6 (1.0 to 2.6) |
| c60 | 40 | 0.4 (0.2 to 1.1) | 1.3 (0.7 to 1.7) |
| taxol | 40 | 1.7 (1.0 to 3.5) | 2.2 (1.1 to 3.8) |
| paracetamol_cluster | 34 | 5.4 (2.8 to 10.5) | 4.8 (2.1 to 10.7) |
| crambin | 26 | 10.8 (6.0 to 21.6) | 10.2 (5.9 to 20.6) |
| ubiquitin | 16 | 16.7 (9.8 to 30.8) | 17.2 (10.8 to 31.3) |
| **all** | **196** | **2.3** (0.2 to 30.8) | **3.1** (0.7 to 31.3) |

The advantage grows with the molecule, which is what a screening win looks like:
the fraction of atom pairs which survive falls as the molecule grows, and the
reference computes all of them. It falls with the diffuse sets for the same
reason in reverse.

### The block size the heuristic picked was wrong at both ends

A sweep of the `block_size` argument over 720 configurations, six molecules and
fifteen bases at a threshold of 1.0e-14. The number of blocks is what the driver
actually formed. `auto` is the heuristic's own choice.

c60 / def2-qzvpd, 3780 nao

| block size | auto | 256 | 512 | 1024 | 2048 | 4096 | 8192 | 16384 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blocks | 1 | 7 | 4 | 2 | 1 | 1 | 1 | 1 |
| ms | 11.90 | 2.71 | 3.73 | 6.30 | 11.12 | 11.50 | 11.56 | 11.11 |

c60 / cc-pv6z, 8400 nao

| block size | auto | 256 | 512 | 1024 | 2048 | 4096 | 8192 | 16384 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blocks | 1 | 7 | 4 | 2 | 1 | 1 | 1 | 1 |
| ms | 90.26 | 19.56 | 30.64 | 49.40 | 90.46 | 91.05 | 89.81 | 90.52 |

tagrisso / cc-pv6z, 8183 nao

| block size | auto | 256 | 512 | 1024 | 2048 | 4096 | 8192 | 16384 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blocks | 10 | 16 | 12 | 10 | 10 | 10 | 10 | 10 |
| ms | 16.07 | 11.21 | 15.77 | 15.96 | 15.87 | 15.86 | 15.88 | 15.66 |

ubiquitin / def2-qzvpd, 59053 nao

| block size | auto | 256 | 512 | 1024 | 2048 | 4096 | 8192 | 16384 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blocks | 37 | 2964 | 1486 | 747 | 379 | 195 | 102 | 56 |
| ms | 64.39 | 146.33 | 103.59 | 79.13 | 65.74 | 57.92 | 53.59 | 53.22 |

ubiquitin / cc-pv6z, 141523 nao

| block size | auto | 256 | 512 | 1024 | 2048 | 4096 | 8192 | 16384 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| blocks | 37 | 2964 | 1486 | 747 | 379 | 195 | 102 | 56 |
| ms | 328.05 | 559.02 | 414.35 | 318.13 | 263.62 | 241.50 | 252.81 | 270.21 |

c60 has 1770 atom pairs, which is below `min_block_size = 2048`, so the heuristic
gives it **one block** and the driver runs serially. Forcing a smaller block is
worth up to 4.6 times. ubiquitin fails the other way: the heuristic aims at
`blocks_per_thread * nthreads` blocks, 28, and lands on 37, where the optimum is
56 to 195. Its 256 column is 2.3 times worse than auto, so the per block fixed
cost is real and this is a genuine optimum rather than a monotone preference for
more blocks.

### The thread scan confirms it is serial, not merely slow

Best of six runs, `OMP_NUM_THREADS` set per column, threshold 1.0e-14.

| molecule | basis | block | 1 | 2 | 4 | 8 | 14 | blocks |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| c60 | def2-qzvpd | auto | 11.055 | 11.147 | 11.033 | 11.176 | 11.328 | 1 |
| c60 | def2-qzvpd | 256 | 13.132 | 7.639 | 4.529 | 2.999 | 3.12 | 7 |
| tagrisso | cc-pv6z | auto | 51.493 | 27.373 | 20.539 | 16.432 | 16.663 | 10 |
| tagrisso | cc-pv6z | 256 | 55.58 | 29.311 | 19.375 | 12.457 | 11.637 | 16 |
| crambin | def2-tzvp | auto | 23.688 | 13.176 | 8.142 | 5.338 | 4.088 | 37 |
| ubiquitin | def2-qzvp | auto | 192.416 | 107.434 | 64.429 | 36.698 | 28.421 | 37 |
| ubiquitin | def2-qzvp | 8192 | 229.257 | 116.684 | 66.069 | 36.426 | 28.165 | 102 |

The first row is flat from one thread to fourteen. That is the whole of the c60
defect: not a driver which threads badly, but one which does not thread at all.
The same case at 256 costs 19 per cent on one thread, which is the fixed cost of
seven blocks against one, and returns 3.7 times on eight.

### Fitting the two constants

`CAtomBasisPairGroup::make_block_size` computes
`npairs / (blocks_per_thread * nthreads)` and clamps it up to `min_block_size`.
Both constants were fitted by emulating the formula through the driver's
`block_size` argument and measuring, over ninety molecule and basis cases at
fourteen threads. The objective is the geometric mean of each case's slowdown
against the best any candidate reaches for that case.

| geometric mean | worst case | total ms | blocks_per_thread | min_block_size |
| --- | --- | --- | --- | --- |
| 1.081 | 1.33 | 1131.2 | 4 | 256 |
| 1.081 | 1.41 | 1284.1 | 2 | 256 |
| 1.107 | 1.41 | 1282.6 | 2 | 128 |
| 1.112 | 1.77 | 1156.6 | 4 | 512 |
| 1.113 | 1.77 | 1309.6 | 2 | 512 |
| 1.135 | 1.73 | 1136.7 | 4 | 128 |
| 1.331 | 5.55 | 1463.8 | 2 | 2048 | (previous, rank 29 of 56) |

`(4, 256)` ties `(2, 256)` on the mean and wins on both the worst case and the
total, and is the pair which now ships. It lands on the measured optimum for the cases
which were broken: c60 gets `max(1770/56, 256) = 256`, which is the best column of
its sweep, and ubiquitin gets `757065/56 = 13519`, which is 56 blocks and also its
best.

| molecule | basis | now ms | new ms | gain | size now | size new |
| --- | --- | --- | --- | --- | --- | --- |
| tagrisso | def2-svp | 0.60 | 0.72 | 0.83 | 2048 | 256 |
| tagrisso | def2-qzvpd | 2.80 | 2.21 | 1.26 | 2048 | 256 |
| tagrisso | cc-pv6z | 16.25 | 11.19 | 1.45 | 2048 | 256 |
| c60 | def2-svp | 0.89 | 0.67 | 1.34 | 2048 | 256 |
| c60 | def2-qzvpd | 11.12 | 2.75 | 4.05 | 2048 | 256 |
| c60 | cc-pv6z | 92.69 | 21.11 | 4.39 | 2048 | 256 |
| taxol | def2-svp | 0.73 | 0.96 | 0.76 | 2048 | 256 |
| taxol | def2-qzvpd | 5.12 | 3.71 | 1.38 | 2048 | 256 |
| taxol | cc-pv6z | 32.53 | 17.29 | 1.88 | 2048 | 256 |
| paracetamol_cluster | def2-svp | 1.55 | 1.53 | 1.01 | 2048 | 911 |
| paracetamol_cluster | def2-qzvpd | 10.23 | 11.15 | 0.92 | 2048 | 911 |
| paracetamol_cluster | cc-pv6z | 56.04 | 53.20 | 1.05 | 2048 | 911 |
| crambin | def2-svp | 1.87 | 2.14 | 0.87 | 7348 | 3674 |
| crambin | def2-qzvpd | 22.91 | 22.47 | 1.02 | 7348 | 3674 |
| crambin | cc-pv6z | 128.90 | 114.47 | 1.13 | 7348 | 3674 |
| ubiquitin | def2-svp | 3.35 | 3.61 | 0.93 | 27038 | 13519 |
| ubiquitin | def2-qzvpd | 66.74 | 53.10 | 1.26 | 27038 | 13519 |
| ubiquitin | cc-pv6z | 325.55 | 251.22 | 1.30 | 27038 | 13519 |
| **all ninety cases** | | **1463.8** | **1131.2** | **1.29** | | |

It is not free. The worst regression is 1.32 times, on taxol in def2-svp, and the
cheap small cases lose a little throughout: tagrisso def2-svp 0.83, taxol def2-svp
0.76, crambin def2-svp 0.87. Those are all sub-millisecond runs where the extra
per block fixed cost is not repaid. The trade is tenths of a millisecond on the
cheapest cases against 4.4 times on c60 and 1.3 times on the largest.

Two limits on the fit. It was measured at fourteen threads on this machine; the
formula divides by the thread count so `blocks_per_thread` should carry, while
`min_block_size` is an absolute floor and is the more machine specific of the two.
And it was fitted on the overlap alone, while the same constants serve every
driver built on `sparsity::make_pattern`, whose per block costs differ. The
earlier note in `SparsityPattern.hpp`, that two blocks per thread beat four, was
measured against the hand written kernels and no longer holds for this one.

### What these numbers say

The driver beats the reference by 3.0 times in the geometric mean over the rows
which have one, from 1.0 to 31 times, and the advantage tracks molecule size
rather than basis size. Under the previous constants thirty six of c60's forty
rows were slower than the reference; under the ones which now ship, three are, all
within 20 per cent.

Spot checks after the constants were changed, against the numbers above: c60 in
aug-cc-pV6Z 245.8 ms to 52.8 ms, c60 in def2-qzvpd 11.1 ms to 2.96 ms, ubiquitin in
cc-pV6Z 357.7 ms to 242.4 ms, and the block counts move from 1 to 7 for c60 and
from 37 to 64 for ubiquitin.

Loosening the threshold from 1.0e-14 to 1.0e-12 is worth 5 to 20 per cent almost
everywhere. It is a real effect and a second order one next to the block size.

Two rows should not be read. ubiquitin in aug-cc-pV5Z at 1.0e-14 shows the fitted
size slower than auto, 682 ms against 621, and aug-cc-pV6Z at 1.0e-12 shows them
level, both against the trend of every neighbouring row. At 127 to 194 thousand
basis functions those allocate ten gigabytes or more per matrix and are bound by
the memory rather than the arithmetic, and the best of three runs is not enough to
separate them.

## The generated kinetic energy kernels

The kinetic energy driver was added on the same shape as the overlap driver, with
the pair of molecular bases dropped: the operator is symmetric in the two sides,
so it takes one basis and describes a symmetric quantity. Its kernels are
generated to l = 6 in both orders, building both sides with vertical recurrences
and reaching the spherical components through a transform rather than a transfer,
and they consume the overlap intermediates on the way.

Measured on an otherwise idle machine at `OMP_NUM_THREADS=14`, best of two to five
runs, with the block size constants fitted on the overlap, `blocks_per_thread = 4`
and `min_block_size = 256`. **ref** is `KineticEnergyDriver`, which computes every
atom pair and carries no threshold, so it is timed once per molecule and basis and
the same number serves both threshold rows.

**The bases stop at g.** `KineticEnergyFunc.hpp` in the reference dispatches only
up to angular momentum four and returns zeros above it, silently, although the
`PrimRec` kernels for h and i are present in `t2c_kinetic_energy`. Every basis
below therefore reaches g at most, asserted per case, so the reference is
computing rather than returning zeros. The correlation consistent sets at quintuple
and sextuple zeta are left out for that reason and not because the driver cannot
do them.

`max |d|` is the largest absolute difference over the whole matrix, filled where
both the sparse and the dense matrix fit under eight thousand basis functions.
`dense too big` marks the cases where the reference cannot be run at all.

### tagrisso

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 683 | 1.70 | 0.79 | 2.2 | 8.94e-15 |
| def2-svp | 1e-12 | 2 | 683 | 1.70 | 0.80 | 2.1 | 8.94e-15 |
| def2-svpd | 1e-14 | 2 | 1010 | 1.89 | 0.87 | 2.2 | 8.94e-15 |
| def2-svpd | 1e-12 | 2 | 1010 | 1.89 | 1.01 | 1.9 | 8.94e-15 |
| def2-tzvp | 1e-14 | 3 | 1345 | 1.74 | 1.30 | 1.3 | 9.10e-15 |
| def2-tzvp | 1e-12 | 3 | 1345 | 1.74 | 1.23 | 1.4 | 9.10e-15 |
| def2-tzvpp | 1e-14 | 3 | 1609 | 1.78 | 1.16 | 1.5 | 9.10e-15 |
| def2-tzvpp | 1e-12 | 3 | 1609 | 1.78 | 1.14 | 1.6 | 9.10e-15 |
| def2-tzvpd | 1e-14 | 3 | 1672 | 1.43 | 1.37 | 1.0 | 9.10e-15 |
| def2-tzvpd | 1e-12 | 3 | 1672 | 1.43 | 1.40 | 1.0 | 9.10e-15 |
| def2-tzvppd | 1e-14 | 3 | 1936 | 2.07 | 1.48 | 1.4 | 9.10e-15 |
| def2-tzvppd | 1e-12 | 3 | 1936 | 2.07 | 1.71 | 1.2 | 9.10e-15 |
| def2-qzvp | 1e-14 | 4 | 3099 | 4.06 | 2.51 | 1.6 | 9.16e-15 |
| def2-qzvp | 1e-12 | 4 | 3099 | 4.06 | 2.51 | 1.6 | 9.16e-15 |
| def2-qzvpp | 1e-14 | 4 | 3099 | 3.82 | 2.58 | 1.5 | 9.16e-15 |
| def2-qzvpp | 1e-12 | 4 | 3099 | 3.82 | 2.50 | 1.5 | 9.16e-15 |
| def2-qzvpd | 1e-14 | 4 | 3426 | 4.29 | 3.12 | 1.4 | 9.16e-15 |
| def2-qzvpd | 1e-12 | 4 | 3426 | 4.29 | 3.04 | 1.4 | 9.16e-15 |
| def2-qzvppd | 1e-14 | 4 | 3426 | 4.30 | 3.08 | 1.4 | 9.16e-15 |
| def2-qzvppd | 1e-12 | 4 | 3426 | 4.30 | 3.06 | 1.4 | 9.16e-15 |
| cc-pvdz | 1e-14 | 2 | 683 | 1.03 | 0.72 | 1.4 | 8.23e-15 |
| cc-pvdz | 1e-12 | 2 | 683 | 1.03 | 0.75 | 1.4 | 8.23e-15 |
| cc-pvtz | 1e-14 | 3 | 1572 | 1.58 | 1.17 | 1.4 | 9.15e-15 |
| cc-pvtz | 1e-12 | 3 | 1572 | 1.58 | 1.11 | 1.4 | 9.15e-15 |
| cc-pvqz | 1e-14 | 4 | 3025 | 4.02 | 2.49 | 1.6 | 9.15e-15 |
| cc-pvqz | 1e-12 | 4 | 3025 | 4.02 | 2.38 | 1.7 | 9.15e-15 |
| aug-cc-pvdz | 1e-14 | 2 | 1148 | 1.20 | 1.00 | 1.2 | 8.23e-15 |
| aug-cc-pvdz | 1e-12 | 2 | 1148 | 1.20 | 0.98 | 1.2 | 8.23e-15 |
| aug-cc-pvtz | 1e-14 | 3 | 2461 | 2.50 | 2.15 | 1.2 | 9.15e-15 |
| aug-cc-pvtz | 1e-12 | 3 | 2461 | 2.50 | 2.10 | 1.2 | 9.15e-15 |
| aug-cc-pvqz | 1e-14 | 4 | 4478 | 8.56 | 5.48 | 1.6 | 9.15e-15 |
| aug-cc-pvqz | 1e-12 | 4 | 4478 | 8.56 | 5.19 | 1.6 | 9.15e-15 |

### c60

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 840 | 1.12 | 0.61 | 1.8 | 1.95e-15 |
| def2-svp | 1e-12 | 2 | 840 | 1.12 | 0.75 | 1.5 | 1.95e-15 |
| def2-svpd | 1e-14 | 2 | 1200 | 1.63 | 1.04 | 1.6 | 1.95e-15 |
| def2-svpd | 1e-12 | 2 | 1200 | 1.63 | 0.81 | 2.0 | 1.95e-15 |
| def2-tzvp | 1e-14 | 3 | 1860 | 2.13 | 1.34 | 1.6 | 8.67e-15 |
| def2-tzvp | 1e-12 | 3 | 1860 | 2.13 | 1.30 | 1.6 | 8.67e-15 |
| def2-tzvpp | 1e-14 | 3 | 1860 | 2.11 | 1.28 | 1.6 | 8.67e-15 |
| def2-tzvpp | 1e-12 | 3 | 1860 | 2.11 | 1.28 | 1.6 | 8.67e-15 |
| def2-tzvpd | 1e-14 | 3 | 2220 | 1.96 | 1.71 | 1.1 | 8.67e-15 |
| def2-tzvpd | 1e-12 | 3 | 2220 | 1.96 | 1.65 | 1.2 | 8.67e-15 |
| def2-tzvppd | 1e-14 | 3 | 2220 | 1.87 | 1.68 | 1.1 | 8.67e-15 |
| def2-tzvppd | 1e-12 | 3 | 2220 | 1.87 | 1.63 | 1.1 | 8.67e-15 |
| def2-qzvp | 1e-14 | 4 | 3420 | 4.94 | 3.56 | 1.4 | 9.17e-15 |
| def2-qzvp | 1e-12 | 4 | 3420 | 4.94 | 3.22 | 1.5 | 9.17e-15 |
| def2-qzvpp | 1e-14 | 4 | 3420 | 6.10 | 3.31 | 1.8 | 9.17e-15 |
| def2-qzvpp | 1e-12 | 4 | 3420 | 6.10 | 3.47 | 1.8 | 9.17e-15 |
| def2-qzvpd | 1e-14 | 4 | 3780 | 6.42 | 4.33 | 1.5 | 9.17e-15 |
| def2-qzvpd | 1e-12 | 4 | 3780 | 6.42 | 4.16 | 1.5 | 9.17e-15 |
| def2-qzvppd | 1e-14 | 4 | 3780 | 6.83 | 3.95 | 1.7 | 9.17e-15 |
| def2-qzvppd | 1e-12 | 4 | 3780 | 6.83 | 4.18 | 1.6 | 9.17e-15 |
| cc-pvdz | 1e-14 | 2 | 840 | 1.54 | 0.74 | 2.1 | 1.78e-15 |
| cc-pvdz | 1e-12 | 2 | 840 | 1.54 | 0.70 | 2.2 | 1.78e-15 |
| cc-pvtz | 1e-14 | 3 | 1800 | 2.08 | 1.41 | 1.5 | 9.16e-15 |
| cc-pvtz | 1e-12 | 3 | 1800 | 2.08 | 1.33 | 1.6 | 9.16e-15 |
| cc-pvqz | 1e-14 | 4 | 3300 | 5.31 | 3.23 | 1.6 | 9.15e-15 |
| cc-pvqz | 1e-12 | 4 | 3300 | 5.31 | 3.13 | 1.7 | 9.15e-15 |
| aug-cc-pvdz | 1e-14 | 2 | 1380 | 1.51 | 1.11 | 1.4 | 1.78e-15 |
| aug-cc-pvdz | 1e-12 | 2 | 1380 | 1.51 | 1.07 | 1.4 | 1.78e-15 |
| aug-cc-pvtz | 1e-14 | 3 | 2760 | 3.30 | 2.63 | 1.3 | 9.16e-15 |
| aug-cc-pvtz | 1e-12 | 3 | 2760 | 3.30 | 2.55 | 1.3 | 9.16e-15 |
| aug-cc-pvqz | 1e-14 | 4 | 4800 | 11.56 | 7.30 | 1.6 | 9.15e-15 |
| aug-cc-pvqz | 1e-12 | 4 | 4800 | 11.56 | 6.90 | 1.7 | 9.15e-15 |

### taxol

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 1099 | 1.12 | 0.83 | 1.4 | 9.13e-15 |
| def2-svp | 1e-12 | 2 | 1099 | 1.12 | 0.83 | 1.4 | 9.13e-15 |
| def2-svpd | 1e-14 | 2 | 1657 | 1.95 | 1.13 | 1.7 | 9.13e-15 |
| def2-svpd | 1e-12 | 2 | 1657 | 1.95 | 1.43 | 1.4 | 9.13e-15 |
| def2-tzvp | 1e-14 | 3 | 2185 | 2.29 | 1.39 | 1.7 | 9.14e-15 |
| def2-tzvp | 1e-12 | 3 | 2185 | 2.29 | 1.36 | 1.7 | 9.14e-15 |
| def2-tzvpp | 1e-14 | 3 | 2577 | 2.97 | 1.58 | 1.9 | 9.14e-15 |
| def2-tzvpp | 1e-12 | 3 | 2577 | 2.97 | 1.52 | 1.9 | 9.14e-15 |
| def2-tzvpd | 1e-14 | 3 | 2743 | 3.01 | 1.92 | 1.6 | 9.14e-15 |
| def2-tzvpd | 1e-12 | 3 | 2743 | 3.01 | 1.86 | 1.6 | 9.14e-15 |
| def2-tzvppd | 1e-14 | 3 | 3135 | 3.14 | 2.16 | 1.5 | 9.14e-15 |
| def2-tzvppd | 1e-12 | 3 | 3135 | 3.14 | 2.14 | 1.5 | 9.14e-15 |
| def2-qzvp | 1e-14 | 4 | 4947 | 9.64 | 3.82 | 2.5 | 9.12e-15 |
| def2-qzvp | 1e-12 | 4 | 4947 | 9.64 | 3.61 | 2.7 | 9.12e-15 |
| def2-qzvpp | 1e-14 | 4 | 4947 | 9.64 | 3.73 | 2.6 | 9.12e-15 |
| def2-qzvpp | 1e-12 | 4 | 4947 | 9.64 | 3.67 | 2.6 | 9.12e-15 |
| def2-qzvpd | 1e-14 | 4 | 5505 | 10.06 | 5.00 | 2.0 | 9.12e-15 |
| def2-qzvpd | 1e-12 | 4 | 5505 | 10.06 | 4.81 | 2.1 | 9.12e-15 |
| def2-qzvppd | 1e-14 | 4 | 5505 | 11.87 | 4.81 | 2.5 | 9.12e-15 |
| def2-qzvppd | 1e-12 | 4 | 5505 | 11.87 | 4.66 | 2.5 | 9.12e-15 |
| cc-pvdz | 1e-14 | 2 | 1099 | 2.49 | 0.90 | 2.8 | 8.96e-15 |
| cc-pvdz | 1e-12 | 2 | 1099 | 2.49 | 0.92 | 2.7 | 8.96e-15 |
| cc-pvtz | 1e-14 | 3 | 2516 | 2.96 | 1.61 | 1.8 | 9.04e-15 |
| cc-pvtz | 1e-12 | 3 | 2516 | 2.96 | 1.54 | 1.9 | 9.04e-15 |
| cc-pvqz | 1e-14 | 4 | 4825 | 8.79 | 3.63 | 2.4 | 9.17e-15 |
| cc-pvqz | 1e-12 | 4 | 4825 | 8.79 | 3.47 | 2.5 | 9.17e-15 |
| aug-cc-pvdz | 1e-14 | 2 | 1844 | 1.80 | 1.35 | 1.3 | 8.96e-15 |
| aug-cc-pvdz | 1e-12 | 2 | 1844 | 1.80 | 1.36 | 1.3 | 8.96e-15 |
| aug-cc-pvtz | 1e-14 | 3 | 3933 | 5.20 | 3.21 | 1.6 | 9.04e-15 |
| aug-cc-pvtz | 1e-12 | 3 | 3933 | 5.20 | 3.14 | 1.7 | 9.04e-15 |
| aug-cc-pvqz | 1e-14 | 4 | 7134 | 17.44 | 8.36 | 2.1 | 9.17e-15 |
| aug-cc-pvqz | 1e-12 | 4 | 7134 | 17.44 | 8.14 | 2.1 | 9.17e-15 |

### paracetamol_cluster

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 3184 | 3.96 | 1.76 | 2.3 | 9.13e-15 |
| def2-svp | 1e-12 | 2 | 3184 | 3.96 | 1.68 | 2.4 | 9.13e-15 |
| def2-svpd | 1e-14 | 2 | 4768 | 7.23 | 3.11 | 2.3 | 9.13e-15 |
| def2-svpd | 1e-12 | 2 | 4768 | 7.23 | 2.98 | 2.4 | 9.13e-15 |
| def2-tzvp | 1e-14 | 3 | 6320 | 13.55 | 3.64 | 3.7 | 9.15e-15 |
| def2-tzvp | 1e-12 | 3 | 6320 | 13.55 | 3.43 | 4.0 | 9.15e-15 |
| def2-tzvpp | 1e-14 | 3 | 7472 | 17.79 | 4.34 | 4.1 | 9.15e-15 |
| def2-tzvpp | 1e-12 | 3 | 7472 | 17.79 | 3.96 | 4.5 | 9.15e-15 |
| def2-tzvpd | 1e-14 | 3 | 7904 | 20.96 | 6.32 | 3.3 | 9.15e-15 |
| def2-tzvpd | 1e-12 | 3 | 7904 | 20.96 | 5.91 | 3.5 | 9.15e-15 |
| def2-tzvppd | 1e-14 | 3 | 9056 | 24.72 | 6.85 | 3.6 | -- |
| def2-tzvppd | 1e-12 | 3 | 9056 | 24.72 | 6.50 | 3.8 | -- |
| def2-qzvp | 1e-14 | 4 | 14352 | 69.56 | 11.15 | 6.2 | -- |
| def2-qzvp | 1e-12 | 4 | 14352 | 69.56 | 10.30 | 6.8 | -- |
| def2-qzvpp | 1e-14 | 4 | 14352 | 69.00 | 11.08 | 6.2 | -- |
| def2-qzvpp | 1e-12 | 4 | 14352 | 69.00 | 10.24 | 6.7 | -- |
| def2-qzvpd | 1e-14 | 4 | 15936 | 86.02 | 17.47 | 4.9 | -- |
| def2-qzvpd | 1e-12 | 4 | 15936 | 86.02 | 15.75 | 5.5 | -- |
| def2-qzvppd | 1e-14 | 4 | 15936 | 85.74 | 17.31 | 5.0 | -- |
| def2-qzvppd | 1e-12 | 4 | 15936 | 85.74 | 15.48 | 5.5 | -- |
| cc-pvdz | 1e-14 | 2 | 3184 | 6.43 | 1.92 | 3.3 | 9.11e-15 |
| cc-pvdz | 1e-12 | 2 | 3184 | 6.43 | 1.84 | 3.5 | 9.11e-15 |
| cc-pvtz | 1e-14 | 3 | 7296 | 19.01 | 4.21 | 4.5 | 9.15e-15 |
| cc-pvtz | 1e-12 | 3 | 7296 | 19.01 | 3.91 | 4.9 | 9.15e-15 |
| cc-pvqz | 1e-14 | 4 | 14000 | 67.72 | 10.39 | 6.5 | -- |
| cc-pvqz | 1e-12 | 4 | 14000 | 67.72 | 9.38 | 7.2 | -- |
| aug-cc-pvdz | 1e-14 | 2 | 5344 | 11.53 | 4.11 | 2.8 | 9.11e-15 |
| aug-cc-pvdz | 1e-12 | 2 | 5344 | 11.53 | 3.92 | 2.9 | 9.11e-15 |
| aug-cc-pvtz | 1e-14 | 3 | 11408 | 41.96 | 11.27 | 3.7 | -- |
| aug-cc-pvtz | 1e-12 | 3 | 11408 | 41.96 | 10.52 | 4.0 | -- |
| aug-cc-pvqz | 1e-14 | 4 | 20704 | 157.24 | 34.51 | 4.6 | -- |
| aug-cc-pvqz | 1e-12 | 4 | 20704 | 157.24 | 31.00 | 5.1 | -- |

### crambin

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 6177 | 14.91 | 2.53 | 5.9 | 9.17e-15 |
| def2-svp | 1e-12 | 2 | 6177 | 14.91 | 2.43 | 6.1 | 9.17e-15 |
| def2-svpd | 1e-14 | 2 | 9294 | 27.24 | 5.35 | 5.1 | -- |
| def2-svpd | 1e-12 | 2 | 9294 | 27.24 | 5.14 | 5.3 | -- |
| def2-tzvp | 1e-14 | 3 | 12063 | 50.94 | 6.25 | 8.2 | -- |
| def2-tzvp | 1e-12 | 3 | 12063 | 50.94 | 5.54 | 9.2 | -- |
| def2-tzvpp | 1e-14 | 3 | 14613 | 69.78 | 6.79 | 10.3 | -- |
| def2-tzvpp | 1e-12 | 3 | 14613 | 69.78 | 6.24 | 11.2 | -- |
| def2-tzvpd | 1e-14 | 3 | 15180 | 74.67 | 12.62 | 5.9 | -- |
| def2-tzvpd | 1e-12 | 3 | 15180 | 74.67 | 11.44 | 6.5 | -- |
| def2-tzvppd | 1e-14 | 3 | 17730 | 92.74 | 14.33 | 6.5 | -- |
| def2-tzvppd | 1e-12 | 3 | 17730 | 92.74 | 12.29 | 7.5 | -- |
| def2-qzvp | 1e-14 | 4 | 28167 | 274.07 | 20.77 | 13.2 | -- |
| def2-qzvp | 1e-12 | 4 | 28167 | 274.07 | 19.15 | 14.3 | -- |
| def2-qzvpp | 1e-14 | 4 | 28167 | 276.43 | 21.59 | 12.8 | -- |
| def2-qzvpp | 1e-12 | 4 | 28167 | 276.43 | 18.72 | 14.8 | -- |
| def2-qzvpd | 1e-14 | 4 | 31284 | -- | 39.92 | -- | -- |
| def2-qzvpd | 1e-12 | 4 | 31284 | -- | 34.03 | -- | -- |
| def2-qzvppd | 1e-14 | 4 | 31284 | -- | 38.93 | -- | -- |
| def2-qzvppd | 1e-12 | 4 | 31284 | -- | 34.41 | -- | -- |
| cc-pvdz | 1e-14 | 2 | 6177 | 23.11 | 3.04 | 7.6 | 1.42e-14 |
| cc-pvdz | 1e-12 | 2 | 6177 | 23.11 | 2.83 | 8.2 | 1.42e-14 |
| cc-pvtz | 1e-14 | 3 | 14244 | 81.09 | 6.70 | 12.1 | -- |
| cc-pvtz | 1e-12 | 3 | 14244 | 81.09 | 6.10 | 13.3 | -- |
| cc-pvqz | 1e-14 | 4 | 27459 | 269.24 | 19.45 | 13.8 | -- |
| cc-pvqz | 1e-12 | 4 | 27459 | 269.24 | 17.20 | 15.7 | -- |
| aug-cc-pvdz | 1e-14 | 2 | 10380 | 43.94 | 8.48 | 5.2 | -- |
| aug-cc-pvdz | 1e-12 | 2 | 10380 | 43.94 | 7.72 | 5.7 | -- |
| aug-cc-pvtz | 1e-14 | 3 | 22311 | 156.07 | 25.99 | 6.0 | -- |
| aug-cc-pvtz | 1e-12 | 3 | 22311 | 156.07 | 23.19 | 6.7 | -- |
| aug-cc-pvqz | 1e-14 | 4 | 40674 | -- | 98.76 | -- | -- |
| aug-cc-pvqz | 1e-12 | 4 | 40674 | -- | 86.79 | -- | -- |

### ubiquitin

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 11577 | 55.21 | 4.40 | 12.5 | -- |
| def2-svp | 1e-12 | 2 | 11577 | 55.21 | 4.04 | 13.7 | -- |
| def2-svpd | 1e-14 | 2 | 17433 | 94.84 | 12.14 | 7.8 | -- |
| def2-svpd | 1e-12 | 2 | 17433 | 94.84 | 10.03 | 9.5 | -- |
| def2-tzvp | 1e-14 | 3 | 22442 | 176.50 | 11.13 | 15.9 | -- |
| def2-tzvp | 1e-12 | 3 | 22442 | 176.50 | 9.93 | 17.8 | -- |
| def2-tzvpp | 1e-14 | 3 | 27479 | 259.35 | 13.56 | 19.1 | -- |
| def2-tzvpp | 1e-12 | 3 | 27479 | 259.35 | 11.47 | 22.6 | -- |
| def2-tzvpd | 1e-14 | 3 | 28298 | 274.51 | 31.25 | 8.8 | -- |
| def2-tzvpd | 1e-12 | 3 | 28298 | 274.51 | 26.43 | 10.4 | -- |
| def2-tzvppd | 1e-14 | 3 | 33335 | -- | 36.56 | -- | -- |
| def2-tzvppd | 1e-12 | 3 | 33335 | -- | 30.75 | -- | -- |
| def2-qzvp | 1e-14 | 4 | 53197 | -- | 50.25 | -- | -- |
| def2-qzvp | 1e-12 | 4 | 53197 | -- | 40.43 | -- | -- |
| def2-qzvpp | 1e-14 | 4 | 53197 | -- | 48.06 | -- | -- |
| def2-qzvpp | 1e-12 | 4 | 53197 | -- | 40.90 | -- | -- |
| def2-qzvpd | 1e-14 | 4 | 59053 | -- | 112.50 | -- | -- |
| def2-qzvpd | 1e-12 | 4 | 59053 | -- | 97.71 | -- | -- |
| def2-qzvppd | 1e-14 | 4 | 59053 | -- | 110.37 | -- | -- |
| def2-qzvppd | 1e-12 | 4 | 59053 | -- | 94.84 | -- | -- |
| cc-pvdz | 1e-14 | 2 | 11577 | 77.27 | 5.12 | 15.1 | -- |
| cc-pvdz | 1e-12 | 2 | 11577 | 77.27 | 4.79 | 16.1 | -- |
| cc-pvtz | 1e-14 | 3 | 26870 | 266.12 | 13.03 | 20.4 | -- |
| cc-pvtz | 1e-12 | 3 | 26870 | 266.12 | 11.00 | 24.2 | -- |
| cc-pvqz | 1e-14 | 4 | 51984 | -- | 42.57 | -- | -- |
| cc-pvqz | 1e-12 | 4 | 51984 | -- | 36.29 | -- | -- |
| aug-cc-pvdz | 1e-14 | 2 | 19511 | 150.63 | 19.88 | 7.6 | -- |
| aug-cc-pvdz | 1e-12 | 2 | 19511 | 150.63 | 17.46 | 8.6 | -- |
| aug-cc-pvtz | 1e-14 | 3 | 42163 | -- | 82.23 | -- | -- |
| aug-cc-pvtz | 1e-12 | 3 | 42163 | -- | 69.78 | -- | -- |
| aug-cc-pvqz | 1e-14 | 4 | 77098 | -- | 321.00 | -- | -- |
| aug-cc-pvqz | 1e-12 | 4 | 77098 | -- | 271.11 | -- | -- |

### The advantage against the reference

| molecule | rows | x ref |
| --- | --- | --- |
| tagrisso | 32 | 1.5 (1.0 to 2.2) |
| c60 | 32 | 1.5 (1.1 to 2.2) |
| taxol | 32 | 1.9 (1.3 to 2.8) |
| paracetamol_cluster | 32 | 4.1 (2.3 to 7.2) |
| crambin | 26 | 8.5 (5.1 to 15.7) |
| ubiquitin | 16 | 13.4 (7.6 to 24.2) |
| **all** | **170** | **3.0** (1.0 to 24.2) |

### What these numbers say

The advantage tracks the size of the molecule, as it does for the overlap, and no
case is slower than the reference. The smallest ratio is 1.0 and it belongs to
tagrisso in def2-tzvpd, which is a 1.4 millisecond run.

It is consistently a smaller win than the overlap on the same case: ubiquitin in
cc-pVTZ is 24.2 here against 29.5 for the overlap, crambin in def2-qzvp 14.3
against 18.2. The kinetic kernels carry the overlap intermediates as well as their
own recurrence, so there is more arithmetic behind each surviving atom pair, while
the reference's work grows in much the same way for both operators.

Loosening the threshold from 1.0e-14 to 1.0e-12 is worth 5 to 20 per cent, the
same as for the overlap.

The c60 rows are worth reading against the overlap section. Under the block size
constants which preceded this measurement c60 was a single block and the overlap
lost to the reference on thirty six of its forty rows. Here it wins on all
thirty two, between 1.1 and 2.2 times, which is the refitted floor rather than
anything about the kinetic kernels.

### What the kernels are checked against

The same center integrals go through a closed form rather than a kernel:
(2 l + 3) a b / (a + b) times the overlap of the two primitives, which follows
from applying the operator to a solid harmonic Gaussian. It reproduces the
reference exactly for l = 0 to 4, to 1.8e-15.

The combinations above g cannot be checked against the reference at all. They were
checked instead against the operator identity

    T = (2 l_b + 3) b S + 2 b^2 dS/db

which holds because the operator on a solid harmonic Gaussian is
`(2l + 3) b - 2 b^2 r^2` times the function itself, and `r^2 exp(-b r^2)` is the
derivative of the exponential in the exponent. The derivative was taken by central
difference on the overlap driver, which does reach l = 6. The thirteen
combinations involving h or i agree with it to between 2e-10 and 6e-8 relative,
against a control of 5e-10 on the (s|s) combination, which is itself exact against
the reference. That is the accuracy of the finite difference and not of the
kernels.

## The generated two-center Coulomb kernels

The two-center Coulomb driver takes one molecular basis and returns a
`CPackedMatrix`, the lower triangle of a symmetric matrix held in full. There is
no sparsity pattern and no threshold: the operator decays as the inverse of the
interatomic distance, so no atom pair of any molecule falls below a threshold and
the matrix is dense however large the molecule is. Only the geometry half of
`sparsity::make_blocks` is formed, to divide the atom pairs for the threads.

Measured on an otherwise idle machine at `OMP_NUM_THREADS=14`, best of two to five
runs. **ref** is `TwoCenterElectronRepulsionDriver`. `max |d|` is the largest
absolute difference over the whole matrix.

The bases are the def2 universal fitting sets, which is what this operator is
used with. Both reach angular momentum four.

| molecule | basis | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- |
| tagrisso | jfit | 2176 | 2.59 | 3.32 | 0.8 | 3.48e-13 |
| tagrisso | jkfit | 3387 | 5.17 | 5.82 | 0.9 | 5.12e-13 |
| c60 | jfit | 2940 | 4.00 | 4.13 | 1.0 | 1.28e-13 |
| c60 | jkfit | 4500 | 6.68 | 7.79 | 0.9 | 5.51e-13 |
| taxol | jfit | 3528 | 4.59 | 5.46 | 0.8 | 4.12e-13 |
| taxol | jkfit | 5489 | 12.31 | 9.66 | 1.3 | 5.97e-13 |
| paracetamol_cluster | jfit | 10208 | 40.50 | 34.11 | 1.2 | 4.05e-13 |
| paracetamol_cluster | jkfit | 15888 | 78.19 | 74.08 | 1.1 | 4.97e-13 |
| crambin | jfit | 19500 | 128.64 | 123.19 | 1.0 | 3.55e-13 |
| crambin | jkfit | 30751 | 591.28 | 283.95 | 2.1 | -- |
| ubiquitin | jfit | 36419 | 627.88 | 492.91 | 1.3 | -- |
| ubiquitin | jkfit | 56971 | -- | 1076.53 | -- | -- |

Two cases carry a caveat. crambin in jkfit and ubiquitin in jfit gave the
reference one or two runs only, as each allocates seven to ten gigabytes of dense
matrix, so those two ratios are the softest here. ubiquitin in jkfit has no
reference at all: twenty four gigabytes dense against twelve packed.

### What these numbers say

This is a different picture from the overlap and the kinetic energy, and the
reason is structural. Those two win by not computing the atom pairs the reference
computes. Here both compute every pair, so the comparison is arithmetic against
arithmetic and the ratio sits near one. Five of the ten measured cases are at or
below the reference, and only the two largest pull ahead, crambin in jkfit at 2.1
and ubiquitin in jfit at 1.3.

The agreement is 1e-13 to 6e-13 absolute against integrals of order one hundred,
which is about 1e-15 relative and is where a scheme built on the Boys function
should sit.

### The block size is not what limits this driver

A sweep of the `block_size` argument over the twelve cases. `auto` is the
heuristic's own choice, which is 256 for every molecule.

| molecule | basis | auto | 32 | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096 | 8192 | 16384 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | jfit | 3.51 | 3.29 | 2.84 | 2.62 | 3.24 | 5.14 | 5.00 | 4.86 | 5.16 | 5.11 | 5.11 |
| tagrisso | jkfit | 5.72 | 5.26 | 4.75 | 4.56 | 5.67 | 9.28 | 9.24 | 9.12 | 9.10 | 9.16 | 9.39 |
| c60 | jfit | 4.04 | 4.52 | 3.88 | 3.68 | 4.08 | 6.23 | 10.96 | 20.19 | 20.09 | 20.03 | 20.04 |
| c60 | jkfit | 7.75 | 8.11 | 7.22 | 7.17 | 7.69 | 12.13 | 21.24 | 40.28 | 40.25 | 39.73 | 39.82 |
| taxol | jfit | 5.58 | 6.75 | 5.39 | 4.89 | 5.32 | 7.05 | 12.27 | 12.15 | 12.11 | 12.15 | 12.13 |
| taxol | jkfit | 9.96 | 11.79 | 9.76 | 9.30 | 10.12 | 13.79 | 24.24 | 24.15 | 23.99 | 24.20 | 23.88 |
| paracetamol_cluster | jfit | 33.70 | 54.84 | 42.05 | 37.35 | 34.60 | 35.73 | 41.71 | 61.17 | 107.11 | 117.57 | 121.07 |
| paracetamol_cluster | jkfit | 72.77 | 101.87 | 81.51 | 75.46 | 70.92 | 77.93 | 93.90 | 152.59 | 227.37 | 266.92 | 249.09 |
| crambin | jfit | 122.01 | 199.34 | 153.65 | 135.05 | 121.25 | 119.26 | 130.10 | 145.52 | 161.60 | 189.31 | 240.59 |
| crambin | jkfit | 258.27 | 356.48 | 287.42 | 272.12 | 253.90 | 261.99 | 289.60 | 331.52 | 379.80 | 442.55 | 579.93 |
| ubiquitin | jfit | 512.76 | 767.66 | 597.25 | 528.02 | 488.24 | 481.71 | 526.43 | 575.96 | 609.02 | 730.12 | 816.56 |
| ubiquitin | jkfit | 1071.08 | 1361.77 | 1139.90 | 1091.31 | 1068.72 | 1116.13 | 1285.46 | 1341.58 | 1494.20 | 1775.46 | 1954.27 |

Fitting the ceiling against this, with the geometric mean of each case's slowdown
against the best any ceiling reaches for it:

| max_block_size | geometric mean | worst case | total ms | sizes it picks |
| --- | --- | --- | --- | --- |
| 128 | 1.038 | 1.13 | 2171.5 | 128 |
| 256 (current) | 1.069 | 1.24 | 2073.7 | 256 |
| 512 | 1.084 | 1.24 | 2128.9 | 256, 512 |
| 64 | 1.116 | 1.29 | 2335.6 | 64 |
| 32 | 1.355 | 1.67 | 2881.7 | 32 |

Every ceiling from 64 to 512 lands within four to eight per cent of the best any
of them reaches, where the wrong floor cost the overlap driver a factor of four.
A ceiling of 128 wins the mean and 256 wins the total time, and they disagree
because the effect splits by molecule: 128 buys 7 to 24 per cent on the three
small molecules and costs 6 to 10 per cent on the three large ones, which carry
the time. The value stays at 256 for that reason. Ceilings above 512 were not
fitted, as the curves are already rising there.

One thing the fit turned up which is not about performance. `max_block_size` is
256 and `sparsity::min_block_size` is 256 as well, so the clamp returns 256 for
every molecule and the term which follows the size of the molecule never applies.
The driver divides tagrisso and ubiquitin into blocks of the same size. The
measurement says that costs it little, but it is a coincidence of two constants
fitted for different drivers rather than a choice, and refitting the floor for the
overlap moves this driver with it.

### The RI fitting sets of the correlation consistent bases

The two-center Coulomb kernels were extended to angular momentum eight, which is
what the fitting sets of a transition metal reach: cc-pVQZ-rifit reaches k on a
copper complex and cc-pV5Z-rifit reaches l. Measured on the same build and the
same terms as the fitting set table above, with the copper complex added to the
molecules.

Three kinds of case carry no reference and say so in the row rather than being
left out. `ref past i` is where the combination reaches k or l: the reference
dispatches only to angular momentum six and returns an all zero block above it,
silently, so comparing against it would report a false agreement. `dense big` is
where its dense return passes forty thousand basis functions. And copper is absent
from the DZ and 6Z fitting sets, so those combinations do not exist at all.

| molecule | basis | lmax | nao | packed GB | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | cc-pvdz | f | 2534 | 0.02 | 3.92 | 3.53 | 1.1 | 2.03e-13 |
| tagrisso | cc-pvtz | g | 3987 | 0.06 | 5.30 | 6.26 | 0.8 | 2.98e-13 |
| tagrisso | cc-pvqz | h | 6699 | 0.17 | 14.83 | 19.70 | 0.8 | 2.06e-13 |
| tagrisso | cc-pv5z | i | 10144 | 0.38 | 43.81 | 62.21 | 0.7 | 4.83e-13 |
| tagrisso | cc-pv6z | k | 15091 | 0.85 | ref past i | 186.32 | -- | -- |
| tagrisso | aug-cc-pvdz | f | 3423 | 0.04 | 3.63 | 4.48 | 0.8 | 1.34e-12 |
| tagrisso | aug-cc-pvtz | g | 5440 | 0.11 | 8.04 | 10.97 | 0.7 | 8.81e-13 |
| tagrisso | aug-cc-pvqz | h | 8856 | 0.29 | 26.86 | 36.42 | 0.7 | 9.52e-13 |
| tagrisso | aug-cc-pv5z | i | 13145 | 0.64 | 77.58 | 107.76 | 0.7 | 1.22e-12 |
| tagrisso | aug-cc-pv6z | k | 19076 | 1.36 | ref past i | 349.33 | -- | -- |
| c60 | cc-pvdz | f | 3360 | 0.04 | 2.74 | 3.46 | 0.8 | 2.10e-13 |
| c60 | cc-pvtz | g | 4860 | 0.09 | 11.25 | 7.48 | 1.5 | 2.42e-13 |
| c60 | cc-pvqz | h | 7920 | 0.23 | 23.35 | 24.11 | 1.0 | 2.86e-13 |
| c60 | cc-pv5z | i | 11580 | 0.50 | 56.55 | 70.63 | 0.8 | 3.13e-13 |
| c60 | cc-pv6z | k | 16980 | 1.07 | ref past i | 215.34 | -- | -- |
| c60 | aug-cc-pvdz | f | 4320 | 0.07 | 4.60 | 5.55 | 0.8 | 9.38e-13 |
| c60 | aug-cc-pvtz | g | 6360 | 0.15 | 11.81 | 13.07 | 0.9 | 9.95e-13 |
| c60 | aug-cc-pvqz | h | 10080 | 0.38 | 34.70 | 43.79 | 0.8 | 6.25e-13 |
| c60 | aug-cc-pv5z | i | 14520 | 0.79 | 105.05 | 124.75 | 0.8 | 1.19e-12 |
| c60 | aug-cc-pv6z | k | 20820 | 1.61 | ref past i | 368.03 | -- | -- |
| taxol | cc-pvdz | f | 4102 | 0.06 | 3.86 | 5.02 | 0.8 | 2.20e-13 |
| taxol | cc-pvtz | g | 6411 | 0.15 | 10.01 | 11.41 | 0.9 | 2.20e-13 |
| taxol | cc-pvqz | h | 10747 | 0.43 | 34.07 | 39.07 | 0.9 | 2.33e-13 |
| taxol | cc-pv5z | i | 16232 | 0.98 | 102.70 | 120.43 | 0.9 | 6.25e-13 |
| taxol | cc-pv6z | k | 24123 | 2.17 | ref past i | 362.50 | -- | -- |
| taxol | aug-cc-pvdz | f | 5519 | 0.11 | 7.25 | 9.04 | 0.8 | 1.19e-12 |
| taxol | aug-cc-pvtz | g | 8720 | 0.28 | 19.04 | 21.58 | 0.9 | 8.95e-13 |
| taxol | aug-cc-pvqz | h | 14168 | 0.75 | 63.95 | 73.69 | 0.9 | 1.02e-12 |
| taxol | aug-cc-pv5z | i | 20985 | 1.64 | 198.29 | 221.75 | 0.9 | -- |
| taxol | aug-cc-pv6z | k | 30428 | 3.45 | ref past i | 670.20 | -- | -- |
| Cu_PPh3_4_cation | cc-pvtz | i | 8364 | 0.26 | 19.39 | 19.01 | 1.0 | 1.92e-13 |
| Cu_PPh3_4_cation | cc-pvqz | k | 13798 | 0.71 | ref past i | 63.85 | -- | -- |
| Cu_PPh3_4_cation | cc-pv5z | l | 20756 | 1.60 | ref past i | 196.52 | -- | -- |
| Cu_PPh3_4_cation | aug-cc-pvtz | i | 11273 | 0.47 | 37.13 | 34.44 | 1.1 | 8.81e-13 |
| Cu_PPh3_4_cation | aug-cc-pvqz | k | 18098 | 1.22 | ref past i | 118.79 | -- | -- |
| Cu_PPh3_4_cation | aug-cc-pv5z | l | 26721 | 2.66 | ref past i | 360.91 | -- | -- |
| paracetamol_cluster | cc-pvdz | f | 11872 | 0.53 | 30.54 | 33.76 | 0.9 | 1.78e-13 |
| paracetamol_cluster | cc-pvtz | g | 18576 | 1.29 | 84.12 | 89.27 | 0.9 | 2.49e-13 |
| paracetamol_cluster | cc-pvqz | h | 31152 | 3.62 | 626.88 | 354.30 | 1.8 | -- |
| paracetamol_cluster | cc-pv5z | i | 47072 | 8.25 | dense big | 1173.06 | -- | -- |
| paracetamol_cluster | aug-cc-pvdz | f | 15984 | 0.95 | 57.35 | 61.19 | 0.9 | 7.67e-13 |
| paracetamol_cluster | aug-cc-pvtz | g | 25280 | 2.38 | 165.01 | 177.60 | 0.9 | -- |
| paracetamol_cluster | aug-cc-pvqz | h | 41088 | 6.29 | dense big | 692.77 | -- | -- |
| paracetamol_cluster | aug-cc-pv5z | i | 60880 | 13.81 | dense big | 2128.54 | -- | -- |
| crambin | cc-pvdz | f | 22842 | 1.94 | 124.29 | 123.96 | 1.0 | -- |
| crambin | cc-pvtz | g | 36183 | 4.88 | 580.06 | 391.00 | 1.5 | -- |
| crambin | cc-pvqz | h | 60645 | 13.70 | dense big | 1357.36 | -- | -- |
| crambin | aug-cc-pvdz | f | 30909 | 3.56 | 690.78 | 229.29 | 3.0 | -- |
| crambin | aug-cc-pvtz | g | 49398 | 9.09 | dense big | 767.90 | -- | -- |
| ubiquitin | cc-pvdz | f | 42538 | 6.74 | dense big | 506.90 | -- | -- |
| ubiquitin | aug-cc-pvdz | f | 57831 | 12.46 | dense big | 951.00 | -- | -- |

### What these numbers say

The driver is **slower than the reference on most of these**, typically 0.7 to 0.9
times: 26 of the 34 cases which have a reference are below one. That is the same
picture as the def2 fitting sets, where the ratio
hovered at one, and it has the same cause. Nothing screens here, so both sides
compute every atom pair and the comparison is kernel against kernel rather than
work against work.

The wins are concentrated where the matrix is large, crambin in aug-cc-pVDZ at 3.0,
paracetamol in cc-pVQZ at 1.8 and crambin in cc-pVTZ at 1.5. Those are also the
cases where the dense allocation of the reference starts to weigh, so part of that
gain is the packed storage rather than the integrals.

The agreement is 1.8e-13 to 1.3e-12 absolute against integrals of order one
hundred, about 1e-14 relative, wherever both matrices fit at once.

The copper complex is what these kernels were extended for and it runs: 63.9 ms at
13798 functions reaching k, and 196.5 ms at 20756 reaching l. Neither has a
reference to check against, so what stands behind them is the zero separation
limit below.

### What the kernels above angular momentum six are checked against

The reference cannot reach k or l, so the thirty two combinations which involve
them were checked against the value the two-center integral approaches as the two
atoms meet. That limit is the closed formula in the exponents alone which
`one_center_electron_repulsion` carries, and the kernels of the atom pairs do not
use it, so it is an independent reference rather than a restatement.

The largest difference between the diagonal of the (l|l) block and the closed form,
at three interatomic distances, three primitives per function:

| l | closed form | R = 0.01 | R = 0.001 | R = 0.0001 |
| --- | --- | --- | --- | --- |
| s | 24.8977682425 | 2.09e-04 | 2.09e-06 | 2.09e-08 |
| i | 1.9194163469 | 3.16e-04 | 3.16e-06 | 3.16e-08 |
| k | 1.6667778485 | 3.16e-04 | 3.16e-06 | 3.16e-08 |
| l | 1.4733255801 | 3.15e-04 | 3.15e-06 | 3.15e-08 |

Every order converges as the square of the distance with the same coefficient, and
k and l behave exactly as the orders the reference does validate. The rows for s to
i are checked both ways, which is what makes the method trustworthy here rather
than merely self consistent.

## The buffer the kernels write into

The three two-center drivers were reworked so that a kernel no longer allocates
the buffer it computes in, and no longer zeroes all of it. Two changes, in that
order.

The buffer used to be formed and destroyed once per combination of basis
functions. Sampling put the allocator at 56 per cent of the working samples of the
overlap driver and 32 per cent of the two-center Coulomb driver: the shapes repeat
over the blocks, and the allocator of the system serializes the aligned requests
across the threads. It is now one arena per thread, sized from the highest angular
momenta any block carries through the generated `number_of_buffer_rows` tables,
with a borrowed view over it that every combination reshapes to the atom pairs it
reaches. Holding it for the whole loop costs no more memory than a block at a time
did, as each thread held one of them at once in any case.

Sizing the view per combination is not a detail. A view stretched to the pairs of
the whole block leaves every row of a combination which reaches fewer of them a
page away from the next, and that alone made the first version of this change 12
per cent **slower** than what it replaced.

With the allocator gone, `__bzero` became the largest single frame in all three
drivers at close to 14 per cent of the working samples. Only the rows
`contract_primitives` accumulates into with `+=` have to start at zero; every
other row is written outright before it is read. There is no `+=` into a buffer
row anywhere in `simd_t2c_transfer` or `simd_t2c_transform`, and every kernel makes
exactly one contraction call, so the rows to zero are one contiguous span. The
generator now emits its first row and its length, and `prepare_buffer` fills that
span alone -- three rows in a hundred for two-center Coulomb, seven for kinetic
energy, eight and a half for overlap.

Measured over the same twenty four cases in one session, best of five, the two
changes together are 1.21 times: 1.23 for overlap, 1.20 for kinetic energy and
1.21 for two-center Coulomb.

### Where the time goes now

Per cent of the working samples, idle threads excluded, `OMP_NUM_THREADS=14`.

| bucket | overlap, taxol cc-pV6Z | kinetic, crambin cc-pVQZ | Coulomb, taxol cc-pV6Z-RIFIT |
| --- | --- | --- | --- |
| kernels | 17.4 | 60.8 | 77.6 |
| transforms | 56.7 | 10.4 | 5.2 |
| zeroing | 14.1 | 13.8 | 13.9 |
| exp | 3.8 | 7.7 | 0.6 |
| primitives | 2.3 | 3.3 | 0.3 |
| allocator | 2.3 | 1.1 | 0.6 |
| pattern | 1.9 | 1.4 | 0.0 |
| Boys | -- | -- | 1.3 |

That profile is from before the zeroing was narrowed, which is what the row for it
measures. The allocator row is what the arena left behind, down from 56, and from
32 for the Coulomb driver.

Thread scaling on the same three cases is 6.2, 8.4 and 5.7 times on fourteen
threads, which puts the serial part at 5 to 11 per cent by Amdahl. The rest is the
spread between performance and efficiency cores, not a defect of the loop.

### The overlap driver against the reference

`OMP_NUM_THREADS=14`, best of two to five runs. **ref** is `OverlapDriver`, which
computes every atom pair and carries no threshold, so one number serves both
threshold rows. The benchmark also times the driver at an explicitly passed block
size; that column agreed with the default to within a geometric mean of 1.4 per
cent and is left out here.


#### tagrisso

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 683 | 1.60 | 0.82 | 2.0 |
| def2-svp | 1e-12 | 683 | 1.60 | 0.82 | 2.0 |
| def2-svpd | 1e-14 | 1010 | 0.95 | 0.85 | 1.1 |
| def2-svpd | 1e-12 | 1010 | 0.95 | 0.95 | 1.0 |
| def2-tzvp | 1e-14 | 1345 | 1.85 | 1.10 | 1.7 |
| def2-tzvp | 1e-12 | 1345 | 1.85 | 1.06 | 1.7 |
| def2-tzvpp | 1e-14 | 1609 | 1.51 | 0.94 | 1.6 |
| def2-tzvpp | 1e-12 | 1609 | 1.51 | 1.11 | 1.4 |
| def2-tzvpd | 1e-14 | 1672 | 1.26 | 1.23 | 1.0 |
| def2-tzvpd | 1e-12 | 1672 | 1.26 | 1.23 | 1.0 |
| def2-tzvppd | 1e-14 | 1936 | 1.81 | 1.44 | 1.3 |
| def2-tzvppd | 1e-12 | 1936 | 1.81 | 1.41 | 1.3 |
| def2-qzvp | 1e-14 | 3099 | 3.02 | 1.73 | 1.7 |
| def2-qzvp | 1e-12 | 3099 | 3.02 | 2.17 | 1.4 |
| def2-qzvpp | 1e-14 | 3099 | 3.52 | 1.72 | 2.0 |
| def2-qzvpp | 1e-12 | 3099 | 3.52 | 1.83 | 1.9 |
| def2-qzvpd | 1e-14 | 3426 | 3.35 | 2.02 | 1.7 |
| def2-qzvpd | 1e-12 | 3426 | 3.35 | 2.59 | 1.3 |
| def2-qzvppd | 1e-14 | 3426 | 3.23 | 2.04 | 1.6 |
| def2-qzvppd | 1e-12 | 3426 | 3.23 | 2.58 | 1.3 |
| cc-pvdz | 1e-14 | 683 | 1.68 | 0.71 | 2.4 |
| cc-pvdz | 1e-12 | 683 | 1.68 | 0.86 | 2.0 |
| cc-pvtz | 1e-14 | 1572 | 1.89 | 0.97 | 1.9 |
| cc-pvtz | 1e-12 | 1572 | 1.89 | 1.14 | 1.7 |
| cc-pvqz | 1e-14 | 3025 | 3.90 | 1.65 | 2.4 |
| cc-pvqz | 1e-12 | 3025 | 3.90 | 2.01 | 1.9 |
| cc-pv5z | 1e-14 | 5182 | 8.99 | 3.74 | 2.4 |
| cc-pv5z | 1e-12 | 5182 | 8.99 | 3.58 | 2.5 |
| cc-pv6z | 1e-14 | 8183 | 23.66 | 9.12 | 2.6 |
| cc-pv6z | 1e-12 | 8183 | 23.66 | 8.69 | 2.7 |
| aug-cc-pvdz | 1e-14 | 1148 | 1.03 | 0.93 | 1.1 |
| aug-cc-pvdz | 1e-12 | 1148 | 1.03 | 1.11 | 0.9 |
| aug-cc-pvtz | 1e-14 | 2461 | 2.39 | 1.54 | 1.6 |
| aug-cc-pvtz | 1e-12 | 2461 | 2.39 | 2.02 | 1.2 |
| aug-cc-pvqz | 1e-14 | 4478 | 5.79 | 3.48 | 1.7 |
| aug-cc-pvqz | 1e-12 | 4478 | 5.79 | 3.27 | 1.8 |
| aug-cc-pv5z | 1e-14 | 7339 | 16.33 | 8.66 | 1.9 |
| aug-cc-pv5z | 1e-12 | 7339 | 16.33 | 8.21 | 2.0 |
| aug-cc-pv6z | 1e-14 | 11184 | 47.15 | 22.26 | 2.1 |
| aug-cc-pv6z | 1e-12 | 11184 | 47.15 | 21.49 | 2.2 |

#### c60

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 840 | 0.85 | 0.57 | 1.5 |
| def2-svp | 1e-12 | 840 | 0.85 | 0.69 | 1.2 |
| def2-svpd | 1e-14 | 1200 | 1.73 | 0.70 | 2.5 |
| def2-svpd | 1e-12 | 1200 | 1.73 | 0.91 | 1.9 |
| def2-tzvp | 1e-14 | 1860 | 1.80 | 1.39 | 1.3 |
| def2-tzvp | 1e-12 | 1860 | 1.80 | 1.31 | 1.4 |
| def2-tzvpp | 1e-14 | 1860 | 2.44 | 1.07 | 2.3 |
| def2-tzvpp | 1e-12 | 1860 | 2.44 | 1.34 | 1.8 |
| def2-tzvpd | 1e-14 | 2220 | 2.09 | 1.27 | 1.6 |
| def2-tzvpd | 1e-12 | 2220 | 2.09 | 1.32 | 1.6 |
| def2-tzvppd | 1e-14 | 2220 | 1.83 | 1.29 | 1.4 |
| def2-tzvppd | 1e-12 | 2220 | 1.83 | 1.33 | 1.4 |
| def2-qzvp | 1e-14 | 3420 | 3.64 | 2.49 | 1.5 |
| def2-qzvp | 1e-12 | 3420 | 3.64 | 2.76 | 1.3 |
| def2-qzvpp | 1e-14 | 3420 | 3.41 | 2.45 | 1.4 |
| def2-qzvpp | 1e-12 | 3420 | 3.41 | 2.40 | 1.4 |
| def2-qzvpd | 1e-14 | 3780 | 5.05 | 2.91 | 1.7 |
| def2-qzvpd | 1e-12 | 3780 | 5.05 | 2.79 | 1.8 |
| def2-qzvppd | 1e-14 | 3780 | 4.02 | 2.89 | 1.4 |
| def2-qzvppd | 1e-12 | 3780 | 4.02 | 2.80 | 1.4 |
| cc-pvdz | 1e-14 | 840 | 1.18 | 0.67 | 1.8 |
| cc-pvdz | 1e-12 | 840 | 1.18 | 0.86 | 1.4 |
| cc-pvtz | 1e-14 | 1800 | 1.96 | 1.42 | 1.4 |
| cc-pvtz | 1e-12 | 1800 | 1.96 | 1.41 | 1.4 |
| cc-pvqz | 1e-14 | 3300 | 3.67 | 2.34 | 1.6 |
| cc-pvqz | 1e-12 | 3300 | 3.67 | 2.33 | 1.6 |
| cc-pv5z | 1e-14 | 5460 | 10.31 | 6.50 | 1.6 |
| cc-pv5z | 1e-12 | 5460 | 10.31 | 6.20 | 1.7 |
| cc-pv6z | 1e-14 | 8400 | 25.77 | 18.21 | 1.4 |
| cc-pv6z | 1e-12 | 8400 | 25.77 | 17.08 | 1.5 |
| aug-cc-pvdz | 1e-14 | 1380 | 1.23 | 0.92 | 1.3 |
| aug-cc-pvdz | 1e-12 | 1380 | 1.23 | 1.23 | 1.0 |
| aug-cc-pvtz | 1e-14 | 2760 | 2.77 | 1.89 | 1.5 |
| aug-cc-pvtz | 1e-12 | 2760 | 2.77 | 1.94 | 1.4 |
| aug-cc-pvqz | 1e-14 | 4800 | 7.89 | 5.03 | 1.6 |
| aug-cc-pvqz | 1e-12 | 4800 | 7.89 | 4.99 | 1.6 |
| aug-cc-pv5z | 1e-14 | 7620 | 18.57 | 14.39 | 1.3 |
| aug-cc-pv5z | 1e-12 | 7620 | 18.57 | 13.15 | 1.4 |
| aug-cc-pv6z | 1e-14 | 11340 | 51.30 | 43.36 | 1.2 |
| aug-cc-pv6z | 1e-12 | 11340 | 51.30 | 36.97 | 1.4 |

#### taxol

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 1099 | 0.92 | 0.83 | 1.1 |
| def2-svp | 1e-12 | 1099 | 0.92 | 0.85 | 1.1 |
| def2-svpd | 1e-14 | 1657 | 1.73 | 1.28 | 1.4 |
| def2-svpd | 1e-12 | 1657 | 1.73 | 1.26 | 1.4 |
| def2-tzvp | 1e-14 | 2185 | 2.11 | 1.13 | 1.9 |
| def2-tzvp | 1e-12 | 2185 | 2.11 | 1.45 | 1.5 |
| def2-tzvpp | 1e-14 | 2577 | 1.90 | 1.23 | 1.5 |
| def2-tzvpp | 1e-12 | 2577 | 1.90 | 1.58 | 1.2 |
| def2-tzvpd | 1e-14 | 2743 | 2.02 | 1.46 | 1.4 |
| def2-tzvpd | 1e-12 | 2743 | 2.02 | 2.00 | 1.0 |
| def2-tzvppd | 1e-14 | 3135 | 2.57 | 1.50 | 1.7 |
| def2-tzvppd | 1e-12 | 3135 | 2.57 | 2.17 | 1.2 |
| def2-qzvp | 1e-14 | 4947 | 8.15 | 2.36 | 3.5 |
| def2-qzvp | 1e-12 | 4947 | 8.15 | 2.25 | 3.6 |
| def2-qzvpp | 1e-14 | 4947 | 6.23 | 2.32 | 2.7 |
| def2-qzvpp | 1e-12 | 4947 | 6.23 | 2.25 | 2.8 |
| def2-qzvpd | 1e-14 | 5505 | 7.40 | 2.92 | 2.5 |
| def2-qzvpd | 1e-12 | 5505 | 7.40 | 2.79 | 2.7 |
| def2-qzvppd | 1e-14 | 5505 | 7.50 | 2.92 | 2.6 |
| def2-qzvppd | 1e-12 | 5505 | 7.50 | 2.85 | 2.6 |
| cc-pvdz | 1e-14 | 1099 | 1.17 | 0.86 | 1.4 |
| cc-pvdz | 1e-12 | 1099 | 1.17 | 0.92 | 1.3 |
| cc-pvtz | 1e-14 | 2516 | 3.11 | 1.24 | 2.5 |
| cc-pvtz | 1e-12 | 2516 | 3.11 | 1.22 | 2.5 |
| cc-pvqz | 1e-14 | 4825 | 7.04 | 2.25 | 3.1 |
| cc-pvqz | 1e-12 | 4825 | 7.04 | 2.15 | 3.3 |
| cc-pv5z | 1e-14 | 8246 | 19.80 | 5.21 | 3.8 |
| cc-pv5z | 1e-12 | 8246 | 19.80 | 4.88 | 4.1 |
| cc-pv6z | 1e-14 | 12999 | 55.92 | 13.93 | 4.0 |
| cc-pv6z | 1e-12 | 12999 | 55.92 | 12.99 | 4.3 |
| aug-cc-pvdz | 1e-14 | 1844 | 1.49 | 1.13 | 1.3 |
| aug-cc-pvdz | 1e-12 | 1844 | 1.49 | 1.57 | 0.9 |
| aug-cc-pvtz | 1e-14 | 3933 | 4.14 | 2.18 | 1.9 |
| aug-cc-pvtz | 1e-12 | 3933 | 4.14 | 2.09 | 2.0 |
| aug-cc-pvqz | 1e-14 | 7134 | 12.45 | 5.01 | 2.5 |
| aug-cc-pvqz | 1e-12 | 7134 | 12.45 | 4.77 | 2.6 |
| aug-cc-pv5z | 1e-14 | 11667 | 40.61 | 12.82 | 3.2 |
| aug-cc-pv5z | 1e-12 | 11667 | 40.61 | 12.18 | 3.3 |
| aug-cc-pv6z | 1e-14 | 17752 | 109.93 | 34.82 | 3.2 |
| aug-cc-pv6z | 1e-12 | 17752 | 109.93 | 33.46 | 3.3 |

#### paracetamol_cluster

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 3184 | 3.37 | 1.50 | 2.2 |
| def2-svp | 1e-12 | 3184 | 3.37 | 1.49 | 2.3 |
| def2-svpd | 1e-14 | 4768 | 6.20 | 2.28 | 2.7 |
| def2-svpd | 1e-12 | 4768 | 6.20 | 2.29 | 2.7 |
| def2-tzvp | 1e-14 | 6320 | 11.12 | 2.59 | 4.3 |
| def2-tzvp | 1e-12 | 6320 | 11.12 | 2.50 | 4.4 |
| def2-tzvpp | 1e-14 | 7472 | 14.06 | 2.95 | 4.8 |
| def2-tzvpp | 1e-12 | 7472 | 14.06 | 2.70 | 5.2 |
| def2-tzvpd | 1e-14 | 7904 | 15.50 | 4.39 | 3.5 |
| def2-tzvpd | 1e-12 | 7904 | 15.50 | 4.07 | 3.8 |
| def2-tzvppd | 1e-14 | 9056 | 19.97 | 5.05 | 4.0 |
| def2-tzvppd | 1e-12 | 9056 | 19.97 | 4.41 | 4.5 |
| def2-qzvp | 1e-14 | 14352 | 52.84 | 6.93 | 7.6 |
| def2-qzvp | 1e-12 | 14352 | 52.84 | 6.34 | 8.3 |
| def2-qzvpp | 1e-14 | 14352 | 51.66 | 6.85 | 7.5 |
| def2-qzvpp | 1e-12 | 14352 | 51.66 | 6.25 | 8.3 |
| def2-qzvpd | 1e-14 | 15936 | 64.02 | 10.81 | 5.9 |
| def2-qzvpd | 1e-12 | 15936 | 64.02 | 9.88 | 6.5 |
| def2-qzvppd | 1e-14 | 15936 | 63.11 | 10.29 | 6.1 |
| def2-qzvppd | 1e-12 | 15936 | 63.11 | 10.04 | 6.3 |
| cc-pvdz | 1e-14 | 3184 | 5.54 | 1.73 | 3.2 |
| cc-pvdz | 1e-12 | 3184 | 5.54 | 1.69 | 3.3 |
| cc-pvtz | 1e-14 | 7296 | 14.96 | 2.99 | 5.0 |
| cc-pvtz | 1e-12 | 7296 | 14.96 | 2.78 | 5.4 |
| cc-pvqz | 1e-14 | 14000 | 49.97 | 6.35 | 7.9 |
| cc-pvqz | 1e-12 | 14000 | 49.97 | 5.74 | 8.7 |
| cc-pv5z | 1e-14 | 23936 | 169.71 | 16.37 | 10.4 |
| cc-pv5z | 1e-12 | 23936 | 169.71 | 15.37 | 11.0 |
| cc-pv6z | 1e-14 | 37744 | -- | 46.07 | -- |
| cc-pv6z | 1e-12 | 37744 | -- | 40.51 | -- |
| aug-cc-pvdz | 1e-14 | 5344 | 9.57 | 3.39 | 2.8 |
| aug-cc-pvdz | 1e-12 | 5344 | 9.57 | 3.11 | 3.1 |
| aug-cc-pvtz | 1e-14 | 11408 | 32.95 | 7.86 | 4.2 |
| aug-cc-pvtz | 1e-12 | 11408 | 32.95 | 7.23 | 4.6 |
| aug-cc-pvqz | 1e-14 | 20704 | 108.39 | 21.59 | 5.0 |
| aug-cc-pvqz | 1e-12 | 20704 | 108.39 | 19.10 | 5.7 |
| aug-cc-pv5z | 1e-14 | 33872 | -- | 59.81 | -- |
| aug-cc-pv5z | 1e-12 | 33872 | -- | 53.18 | -- |
| aug-cc-pv6z | 1e-14 | 51552 | -- | 163.69 | -- |
| aug-cc-pv6z | 1e-12 | 51552 | -- | 143.55 | -- |

#### crambin

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 6177 | 12.29 | 2.15 | 5.7 |
| def2-svp | 1e-12 | 6177 | 12.29 | 2.03 | 6.1 |
| def2-svpd | 1e-14 | 9294 | 23.56 | 4.22 | 5.6 |
| def2-svpd | 1e-12 | 9294 | 23.56 | 3.80 | 6.2 |
| def2-tzvp | 1e-14 | 12063 | 41.10 | 4.31 | 9.5 |
| def2-tzvp | 1e-12 | 12063 | 41.10 | 3.97 | 10.4 |
| def2-tzvpp | 1e-14 | 14613 | 59.86 | 4.84 | 12.4 |
| def2-tzvpp | 1e-12 | 14613 | 59.86 | 4.45 | 13.5 |
| def2-tzvpd | 1e-14 | 15180 | 60.73 | 8.66 | 7.0 |
| def2-tzvpd | 1e-12 | 15180 | 60.73 | 8.22 | 7.4 |
| def2-tzvppd | 1e-14 | 17730 | 74.90 | 9.41 | 8.0 |
| def2-tzvppd | 1e-12 | 17730 | 74.90 | 8.57 | 8.7 |
| def2-qzvp | 1e-14 | 28167 | 231.79 | 13.14 | 17.6 |
| def2-qzvp | 1e-12 | 28167 | 231.79 | 11.57 | 20.0 |
| def2-qzvpp | 1e-14 | 28167 | 207.46 | 13.15 | 15.8 |
| def2-qzvpp | 1e-12 | 28167 | 207.46 | 11.71 | 17.7 |
| def2-qzvpd | 1e-14 | 31284 | -- | 23.79 | -- |
| def2-qzvpd | 1e-12 | 31284 | -- | 20.82 | -- |
| def2-qzvppd | 1e-14 | 31284 | -- | 23.51 | -- |
| def2-qzvppd | 1e-12 | 31284 | -- | 22.34 | -- |
| cc-pvdz | 1e-14 | 6177 | 19.56 | 2.50 | 7.8 |
| cc-pvdz | 1e-12 | 6177 | 19.56 | 2.30 | 8.5 |
| cc-pvtz | 1e-14 | 14244 | 62.42 | 5.01 | 12.5 |
| cc-pvtz | 1e-12 | 14244 | 62.42 | 4.57 | 13.7 |
| cc-pvqz | 1e-14 | 27459 | 211.35 | 11.96 | 17.7 |
| cc-pvqz | 1e-12 | 27459 | 211.35 | 10.79 | 19.6 |
| cc-pv5z | 1e-14 | 47106 | -- | 34.82 | -- |
| cc-pv5z | 1e-12 | 47106 | -- | 30.15 | -- |
| cc-pv6z | 1e-14 | 74469 | -- | 102.49 | -- |
| cc-pv6z | 1e-12 | 74469 | -- | 87.32 | -- |
| aug-cc-pvdz | 1e-14 | 10380 | 36.96 | 6.71 | 5.5 |
| aug-cc-pvdz | 1e-12 | 10380 | 36.96 | 6.10 | 6.1 |
| aug-cc-pvtz | 1e-14 | 22311 | 135.83 | 18.31 | 7.4 |
| aug-cc-pvtz | 1e-12 | 22311 | 135.83 | 16.17 | 8.4 |
| aug-cc-pvqz | 1e-14 | 40674 | -- | 52.65 | -- |
| aug-cc-pvqz | 1e-12 | 40674 | -- | 45.38 | -- |
| aug-cc-pv5z | 1e-14 | 66753 | -- | 168.01 | -- |
| aug-cc-pv5z | 1e-12 | 66753 | -- | 140.74 | -- |
| aug-cc-pv6z | 1e-14 | 101832 | -- | 485.45 | -- |
| aug-cc-pv6z | 1e-12 | 101832 | -- | 412.13 | -- |

#### ubiquitin

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 11577 | 42.17 | 3.89 | 10.8 |
| def2-svp | 1e-12 | 11577 | 42.17 | 3.53 | 11.9 |
| def2-svpd | 1e-14 | 17433 | 86.04 | 8.72 | 9.9 |
| def2-svpd | 1e-12 | 17433 | 86.04 | 7.72 | 11.1 |
| def2-tzvp | 1e-14 | 22442 | 136.74 | 8.77 | 15.6 |
| def2-tzvp | 1e-12 | 22442 | 136.74 | 7.32 | 18.7 |
| def2-tzvpp | 1e-14 | 27479 | 190.04 | 9.59 | 19.8 |
| def2-tzvpp | 1e-12 | 27479 | 190.04 | 8.41 | 22.6 |
| def2-tzvpd | 1e-14 | 28298 | 197.27 | 19.22 | 10.3 |
| def2-tzvpd | 1e-12 | 28298 | 197.27 | 17.02 | 11.6 |
| def2-tzvppd | 1e-14 | 33335 | -- | 22.75 | -- |
| def2-tzvppd | 1e-12 | 33335 | -- | 19.60 | -- |
| def2-qzvp | 1e-14 | 53197 | -- | 28.42 | -- |
| def2-qzvp | 1e-12 | 53197 | -- | 23.21 | -- |
| def2-qzvpp | 1e-14 | 53197 | -- | 27.99 | -- |
| def2-qzvpp | 1e-12 | 53197 | -- | 24.23 | -- |
| def2-qzvpd | 1e-14 | 59053 | -- | 60.24 | -- |
| def2-qzvpd | 1e-12 | 59053 | -- | 49.85 | -- |
| def2-qzvppd | 1e-14 | 59053 | -- | 58.83 | -- |
| def2-qzvppd | 1e-12 | 59053 | -- | 49.34 | -- |
| cc-pvdz | 1e-14 | 11577 | 64.56 | 4.45 | 14.5 |
| cc-pvdz | 1e-12 | 11577 | 64.56 | 4.27 | 15.1 |
| cc-pvtz | 1e-14 | 26870 | 224.19 | 9.20 | 24.4 |
| cc-pvtz | 1e-12 | 26870 | 224.19 | 8.20 | 27.3 |
| cc-pvqz | 1e-14 | 51984 | -- | 23.94 | -- |
| cc-pvqz | 1e-12 | 51984 | -- | 21.61 | -- |
| cc-pv5z | 1e-14 | 89381 | -- | 76.57 | -- |
| cc-pv5z | 1e-12 | 89381 | -- | 64.86 | -- |
| cc-pv6z | 1e-14 | 141523 | -- | 251.04 | -- |
| cc-pv6z | 1e-12 | 141523 | -- | 201.15 | -- |
| aug-cc-pvdz | 1e-14 | 19511 | 144.51 | 15.29 | 9.5 |
| aug-cc-pvdz | 1e-12 | 19511 | 144.51 | 13.44 | 10.8 |
| aug-cc-pvtz | 1e-14 | 42163 | -- | 46.23 | -- |
| aug-cc-pvtz | 1e-12 | 42163 | -- | 39.13 | -- |
| aug-cc-pvqz | 1e-14 | 77098 | -- | 166.31 | -- |
| aug-cc-pvqz | 1e-12 | 77098 | -- | 132.75 | -- |
| aug-cc-pv5z | 1e-14 | 126778 | -- | 529.75 | -- |
| aug-cc-pv5z | 1e-12 | 126778 | -- | 433.36 | -- |
| aug-cc-pv6z | 1e-14 | 193665 | -- | 1845.16 | -- |
| aug-cc-pv6z | 1e-12 | 193665 | -- | 1281.76 | -- |

The geometric mean over the 196 cases the reference can run is **3.11**, against
2.86 before this change. 44 cases have no reference at all: it returns a dense
matrix, and ubiquitin in aug-cc-pV6Z would be three hundred gigabytes. The
advantage grows with the molecule, from 1.65 on tagrisso to 14.36 on ubiquitin,
because the driver's gain is the atom pairs it does not compute.

### The kinetic energy driver against the reference

The bases are those whose highest angular momentum is g, where the reference
dispatcher stops. `max abs diff` is over the whole matrix, blank where two dense
matrices do not fit at once.


#### tagrisso

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 683 | 1.63 | 0.77 | 2.1 | 8.94e-15 |
| def2-svp | 1e-12 | 2 | 683 | 1.63 | 0.85 | 1.9 | 8.94e-15 |
| def2-svpd | 1e-14 | 2 | 1010 | 1.30 | 0.90 | 1.4 | 8.94e-15 |
| def2-svpd | 1e-12 | 2 | 1010 | 1.30 | 1.05 | 1.2 | 8.94e-15 |
| def2-tzvp | 1e-14 | 3 | 1345 | 1.20 | 1.11 | 1.1 | 9.10e-15 |
| def2-tzvp | 1e-12 | 3 | 1345 | 1.20 | 1.29 | 0.9 | 9.10e-15 |
| def2-tzvpp | 1e-14 | 3 | 1609 | 1.84 | 1.10 | 1.7 | 9.10e-15 |
| def2-tzvpp | 1e-12 | 3 | 1609 | 1.84 | 1.33 | 1.4 | 9.10e-15 |
| def2-tzvpd | 1e-14 | 3 | 1672 | 1.40 | 1.35 | 1.0 | 9.10e-15 |
| def2-tzvpd | 1e-12 | 3 | 1672 | 1.40 | 1.56 | 0.9 | 9.10e-15 |
| def2-tzvppd | 1e-14 | 3 | 1936 | 1.57 | 1.41 | 1.1 | 9.10e-15 |
| def2-tzvppd | 1e-12 | 3 | 1936 | 1.57 | 1.68 | 0.9 | 9.10e-15 |
| def2-qzvp | 1e-14 | 4 | 3099 | 3.90 | 2.46 | 1.6 | 9.16e-15 |
| def2-qzvp | 1e-12 | 4 | 3099 | 3.90 | 2.33 | 1.7 | 9.16e-15 |
| def2-qzvpp | 1e-14 | 4 | 3099 | 3.77 | 2.46 | 1.5 | 9.16e-15 |
| def2-qzvpp | 1e-12 | 4 | 3099 | 3.77 | 2.37 | 1.6 | 9.16e-15 |
| def2-qzvpd | 1e-14 | 4 | 3426 | 4.43 | 2.96 | 1.5 | 9.16e-15 |
| def2-qzvpd | 1e-12 | 4 | 3426 | 4.43 | 2.89 | 1.5 | 9.16e-15 |
| def2-qzvppd | 1e-14 | 4 | 3426 | 5.17 | 3.01 | 1.7 | 9.16e-15 |
| def2-qzvppd | 1e-12 | 4 | 3426 | 5.17 | 2.91 | 1.8 | 9.16e-15 |
| cc-pvdz | 1e-14 | 2 | 683 | 1.02 | 0.73 | 1.4 | 8.23e-15 |
| cc-pvdz | 1e-12 | 2 | 683 | 1.02 | 0.88 | 1.2 | 8.23e-15 |
| cc-pvtz | 1e-14 | 3 | 1572 | 2.30 | 1.13 | 2.0 | 9.15e-15 |
| cc-pvtz | 1e-12 | 3 | 1572 | 2.30 | 1.09 | 2.1 | 9.15e-15 |
| cc-pvqz | 1e-14 | 4 | 3025 | 4.10 | 2.34 | 1.8 | 9.15e-15 |
| cc-pvqz | 1e-12 | 4 | 3025 | 4.10 | 2.20 | 1.9 | 9.15e-15 |
| aug-cc-pvdz | 1e-14 | 2 | 1148 | 1.18 | 1.01 | 1.2 | 8.23e-15 |
| aug-cc-pvdz | 1e-12 | 2 | 1148 | 1.18 | 1.22 | 1.0 | 8.23e-15 |
| aug-cc-pvtz | 1e-14 | 3 | 2461 | 2.76 | 2.03 | 1.4 | 9.15e-15 |
| aug-cc-pvtz | 1e-12 | 3 | 2461 | 2.76 | 2.02 | 1.4 | 9.15e-15 |
| aug-cc-pvqz | 1e-14 | 4 | 4478 | 8.90 | 5.48 | 1.6 | 9.15e-15 |
| aug-cc-pvqz | 1e-12 | 4 | 4478 | 8.90 | 5.25 | 1.7 | 9.15e-15 |

#### c60

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 840 | 1.11 | 0.62 | 1.8 | 1.95e-15 |
| def2-svp | 1e-12 | 2 | 840 | 1.11 | 0.61 | 1.8 | 1.95e-15 |
| def2-svpd | 1e-14 | 2 | 1200 | 0.99 | 0.78 | 1.3 | 1.95e-15 |
| def2-svpd | 1e-12 | 2 | 1200 | 0.99 | 0.78 | 1.3 | 1.95e-15 |
| def2-tzvp | 1e-14 | 3 | 1860 | 2.18 | 1.36 | 1.6 | 8.67e-15 |
| def2-tzvp | 1e-12 | 3 | 1860 | 2.18 | 1.24 | 1.8 | 8.67e-15 |
| def2-tzvpp | 1e-14 | 3 | 1860 | 2.13 | 1.28 | 1.7 | 8.67e-15 |
| def2-tzvpp | 1e-12 | 3 | 1860 | 2.13 | 1.33 | 1.6 | 8.67e-15 |
| def2-tzvpd | 1e-14 | 3 | 2220 | 2.16 | 1.70 | 1.3 | 8.67e-15 |
| def2-tzvpd | 1e-12 | 3 | 2220 | 2.16 | 1.66 | 1.3 | 8.67e-15 |
| def2-tzvppd | 1e-14 | 3 | 2220 | 2.05 | 1.64 | 1.2 | 8.67e-15 |
| def2-tzvppd | 1e-12 | 3 | 2220 | 2.05 | 1.72 | 1.2 | 8.67e-15 |
| def2-qzvp | 1e-14 | 4 | 3420 | 5.61 | 3.96 | 1.4 | 9.17e-15 |
| def2-qzvp | 1e-12 | 4 | 3420 | 5.61 | 3.86 | 1.5 | 9.17e-15 |
| def2-qzvpp | 1e-14 | 4 | 3420 | 6.28 | 3.74 | 1.7 | 9.17e-15 |
| def2-qzvpp | 1e-12 | 4 | 3420 | 6.28 | 3.75 | 1.7 | 9.17e-15 |
| def2-qzvpd | 1e-14 | 4 | 3780 | 5.43 | 4.66 | 1.2 | 9.17e-15 |
| def2-qzvpd | 1e-12 | 4 | 3780 | 5.43 | 4.47 | 1.2 | 9.17e-15 |
| def2-qzvppd | 1e-14 | 4 | 3780 | 6.26 | 4.64 | 1.3 | 9.17e-15 |
| def2-qzvppd | 1e-12 | 4 | 3780 | 6.26 | 4.48 | 1.4 | 9.17e-15 |
| cc-pvdz | 1e-14 | 2 | 840 | 2.21 | 0.79 | 2.8 | 1.78e-15 |
| cc-pvdz | 1e-12 | 2 | 840 | 2.21 | 0.74 | 3.0 | 1.78e-15 |
| cc-pvtz | 1e-14 | 3 | 1800 | 2.23 | 1.44 | 1.5 | 9.16e-15 |
| cc-pvtz | 1e-12 | 3 | 1800 | 2.23 | 1.38 | 1.6 | 9.16e-15 |
| cc-pvqz | 1e-14 | 4 | 3300 | 4.99 | 3.59 | 1.4 | 9.15e-15 |
| cc-pvqz | 1e-12 | 4 | 3300 | 4.99 | 3.66 | 1.4 | 9.15e-15 |
| aug-cc-pvdz | 1e-14 | 2 | 1380 | 1.46 | 1.10 | 1.3 | 1.78e-15 |
| aug-cc-pvdz | 1e-12 | 2 | 1380 | 1.46 | 1.07 | 1.4 | 1.78e-15 |
| aug-cc-pvtz | 1e-14 | 3 | 2760 | 4.32 | 2.72 | 1.6 | 9.16e-15 |
| aug-cc-pvtz | 1e-12 | 3 | 2760 | 4.32 | 2.60 | 1.7 | 9.16e-15 |
| aug-cc-pvqz | 1e-14 | 4 | 4800 | 10.65 | 8.61 | 1.2 | 9.15e-15 |
| aug-cc-pvqz | 1e-12 | 4 | 4800 | 10.65 | 8.38 | 1.3 | 9.15e-15 |

#### taxol

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 1099 | 0.74 | 0.83 | 0.9 | 9.13e-15 |
| def2-svp | 1e-12 | 2 | 1099 | 0.74 | 0.83 | 0.9 | 9.13e-15 |
| def2-svpd | 1e-14 | 2 | 1657 | 1.80 | 1.09 | 1.7 | 9.13e-15 |
| def2-svpd | 1e-12 | 2 | 1657 | 1.80 | 1.06 | 1.7 | 9.13e-15 |
| def2-tzvp | 1e-14 | 3 | 2185 | 2.13 | 1.30 | 1.6 | 9.14e-15 |
| def2-tzvp | 1e-12 | 3 | 2185 | 2.13 | 1.27 | 1.7 | 9.14e-15 |
| def2-tzvpp | 1e-14 | 3 | 2577 | 2.94 | 1.41 | 2.1 | 9.14e-15 |
| def2-tzvpp | 1e-12 | 3 | 2577 | 2.94 | 1.36 | 2.2 | 9.14e-15 |
| def2-tzvpd | 1e-14 | 3 | 2743 | 2.52 | 1.79 | 1.4 | 9.14e-15 |
| def2-tzvpd | 1e-12 | 3 | 2743 | 2.52 | 1.72 | 1.5 | 9.14e-15 |
| def2-tzvppd | 1e-14 | 3 | 3135 | 3.59 | 1.85 | 1.9 | 9.14e-15 |
| def2-tzvppd | 1e-12 | 3 | 3135 | 3.59 | 1.82 | 2.0 | 9.14e-15 |
| def2-qzvp | 1e-14 | 4 | 4947 | 8.62 | 3.39 | 2.5 | 9.12e-15 |
| def2-qzvp | 1e-12 | 4 | 4947 | 8.62 | 3.24 | 2.7 | 9.12e-15 |
| def2-qzvpp | 1e-14 | 4 | 4947 | 8.41 | 3.31 | 2.5 | 9.12e-15 |
| def2-qzvpp | 1e-12 | 4 | 4947 | 8.41 | 3.14 | 2.7 | 9.12e-15 |
| def2-qzvpd | 1e-14 | 4 | 5505 | 9.80 | 4.22 | 2.3 | 9.12e-15 |
| def2-qzvpd | 1e-12 | 4 | 5505 | 9.80 | 4.05 | 2.4 | 9.12e-15 |
| def2-qzvppd | 1e-14 | 4 | 5505 | 10.03 | 4.21 | 2.4 | 9.12e-15 |
| def2-qzvppd | 1e-12 | 4 | 5505 | 10.03 | 4.04 | 2.5 | 9.12e-15 |
| cc-pvdz | 1e-14 | 2 | 1099 | 2.80 | 0.91 | 3.1 | 8.96e-15 |
| cc-pvdz | 1e-12 | 2 | 1099 | 2.80 | 0.91 | 3.1 | 8.96e-15 |
| cc-pvtz | 1e-14 | 3 | 2516 | 2.68 | 1.45 | 1.8 | 9.04e-15 |
| cc-pvtz | 1e-12 | 3 | 2516 | 2.68 | 1.41 | 1.9 | 9.04e-15 |
| cc-pvqz | 1e-14 | 4 | 4825 | 8.62 | 3.24 | 2.7 | 9.17e-15 |
| cc-pvqz | 1e-12 | 4 | 4825 | 8.62 | 3.03 | 2.8 | 9.17e-15 |
| aug-cc-pvdz | 1e-14 | 2 | 1844 | 2.06 | 1.27 | 1.6 | 8.96e-15 |
| aug-cc-pvdz | 1e-12 | 2 | 1844 | 2.06 | 1.25 | 1.6 | 8.96e-15 |
| aug-cc-pvtz | 1e-14 | 3 | 3933 | 4.99 | 2.77 | 1.8 | 9.04e-15 |
| aug-cc-pvtz | 1e-12 | 3 | 3933 | 4.99 | 2.65 | 1.9 | 9.04e-15 |
| aug-cc-pvqz | 1e-14 | 4 | 7134 | 17.54 | 7.94 | 2.2 | 9.17e-15 |
| aug-cc-pvqz | 1e-12 | 4 | 7134 | 17.54 | 7.53 | 2.3 | 9.17e-15 |

#### paracetamol_cluster

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 3184 | 3.95 | 1.66 | 2.4 | 9.13e-15 |
| def2-svp | 1e-12 | 2 | 3184 | 3.95 | 1.60 | 2.5 | 9.13e-15 |
| def2-svpd | 1e-14 | 2 | 4768 | 7.10 | 2.90 | 2.4 | 9.13e-15 |
| def2-svpd | 1e-12 | 2 | 4768 | 7.10 | 2.76 | 2.6 | 9.13e-15 |
| def2-tzvp | 1e-14 | 3 | 6320 | 13.39 | 3.31 | 4.0 | 9.15e-15 |
| def2-tzvp | 1e-12 | 3 | 6320 | 13.39 | 3.07 | 4.4 | 9.15e-15 |
| def2-tzvpp | 1e-14 | 3 | 7472 | 17.26 | 3.76 | 4.6 | 9.15e-15 |
| def2-tzvpp | 1e-12 | 3 | 7472 | 17.26 | 3.52 | 4.9 | 9.15e-15 |
| def2-tzvpd | 1e-14 | 3 | 7904 | 19.24 | 6.12 | 3.1 | 9.15e-15 |
| def2-tzvpd | 1e-12 | 3 | 7904 | 19.24 | 5.68 | 3.4 | 9.15e-15 |
| def2-tzvppd | 1e-14 | 3 | 9056 | 24.44 | 6.89 | 3.5 | -- |
| def2-tzvppd | 1e-12 | 3 | 9056 | 24.44 | 6.35 | 3.8 | -- |
| def2-qzvp | 1e-14 | 4 | 14352 | 69.43 | 10.92 | 6.4 | -- |
| def2-qzvp | 1e-12 | 4 | 14352 | 69.43 | 9.52 | 7.3 | -- |
| def2-qzvpp | 1e-14 | 4 | 14352 | 68.47 | 10.93 | 6.3 | -- |
| def2-qzvpp | 1e-12 | 4 | 14352 | 68.47 | 9.99 | 6.9 | -- |
| def2-qzvpd | 1e-14 | 4 | 15936 | 82.59 | 16.66 | 5.0 | -- |
| def2-qzvpd | 1e-12 | 4 | 15936 | 82.59 | 15.38 | 5.4 | -- |
| def2-qzvppd | 1e-14 | 4 | 15936 | 82.39 | 16.29 | 5.1 | -- |
| def2-qzvppd | 1e-12 | 4 | 15936 | 82.39 | 15.69 | 5.3 | -- |
| cc-pvdz | 1e-14 | 2 | 3184 | 6.75 | 1.92 | 3.5 | 9.11e-15 |
| cc-pvdz | 1e-12 | 2 | 3184 | 6.75 | 1.84 | 3.7 | 9.11e-15 |
| cc-pvtz | 1e-14 | 3 | 7296 | 18.39 | 3.75 | 4.9 | 9.15e-15 |
| cc-pvtz | 1e-12 | 3 | 7296 | 18.39 | 3.50 | 5.3 | 9.15e-15 |
| cc-pvqz | 1e-14 | 4 | 14000 | 67.12 | 9.86 | 6.8 | -- |
| cc-pvqz | 1e-12 | 4 | 14000 | 67.12 | 8.92 | 7.5 | -- |
| aug-cc-pvdz | 1e-14 | 2 | 5344 | 11.74 | 3.94 | 3.0 | 9.11e-15 |
| aug-cc-pvdz | 1e-12 | 2 | 5344 | 11.74 | 3.79 | 3.1 | 9.11e-15 |
| aug-cc-pvtz | 1e-14 | 3 | 11408 | 40.81 | 11.78 | 3.5 | -- |
| aug-cc-pvtz | 1e-12 | 3 | 11408 | 40.81 | 10.36 | 3.9 | -- |
| aug-cc-pvqz | 1e-14 | 4 | 20704 | 147.40 | 36.54 | 4.0 | -- |
| aug-cc-pvqz | 1e-12 | 4 | 20704 | 147.40 | 33.53 | 4.4 | -- |

#### crambin

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 6177 | 16.79 | 2.44 | 6.9 | 9.17e-15 |
| def2-svp | 1e-12 | 2 | 6177 | 16.79 | 2.28 | 7.4 | 9.17e-15 |
| def2-svpd | 1e-14 | 2 | 9294 | 27.13 | 5.39 | 5.0 | -- |
| def2-svpd | 1e-12 | 2 | 9294 | 27.13 | 5.02 | 5.4 | -- |
| def2-tzvp | 1e-14 | 3 | 12063 | 50.85 | 6.18 | 8.2 | -- |
| def2-tzvp | 1e-12 | 3 | 12063 | 50.85 | 5.31 | 9.6 | -- |
| def2-tzvpp | 1e-14 | 3 | 14613 | 69.13 | 6.75 | 10.2 | -- |
| def2-tzvpp | 1e-12 | 3 | 14613 | 69.13 | 6.06 | 11.4 | -- |
| def2-tzvpd | 1e-14 | 3 | 15180 | 75.01 | 13.21 | 5.7 | -- |
| def2-tzvpd | 1e-12 | 3 | 15180 | 75.01 | 12.59 | 6.0 | -- |
| def2-tzvppd | 1e-14 | 3 | 17730 | 92.11 | 14.81 | 6.2 | -- |
| def2-tzvppd | 1e-12 | 3 | 17730 | 92.11 | 13.32 | 6.9 | -- |
| def2-qzvp | 1e-14 | 4 | 28167 | 277.53 | 23.14 | 12.0 | -- |
| def2-qzvp | 1e-12 | 4 | 28167 | 277.53 | 20.44 | 13.6 | -- |
| def2-qzvpp | 1e-14 | 4 | 28167 | 267.94 | 23.96 | 11.2 | -- |
| def2-qzvpp | 1e-12 | 4 | 28167 | 267.94 | 19.81 | 13.5 | -- |
| def2-qzvpd | 1e-14 | 4 | 31284 | -- | 41.59 | -- | -- |
| def2-qzvpd | 1e-12 | 4 | 31284 | -- | 37.76 | -- | -- |
| def2-qzvppd | 1e-14 | 4 | 31284 | -- | 41.92 | -- | -- |
| def2-qzvppd | 1e-12 | 4 | 31284 | -- | 36.95 | -- | -- |
| cc-pvdz | 1e-14 | 2 | 6177 | 22.61 | 3.03 | 7.5 | 1.42e-14 |
| cc-pvdz | 1e-12 | 2 | 6177 | 22.61 | 2.72 | 8.3 | 1.42e-14 |
| cc-pvtz | 1e-14 | 3 | 14244 | 73.95 | 6.81 | 10.9 | -- |
| cc-pvtz | 1e-12 | 3 | 14244 | 73.95 | 6.06 | 12.2 | -- |
| cc-pvqz | 1e-14 | 4 | 27459 | 267.73 | 20.46 | 13.1 | -- |
| cc-pvqz | 1e-12 | 4 | 27459 | 267.73 | 17.77 | 15.1 | -- |
| aug-cc-pvdz | 1e-14 | 2 | 10380 | 44.14 | 8.86 | 5.0 | -- |
| aug-cc-pvdz | 1e-12 | 2 | 10380 | 44.14 | 7.86 | 5.6 | -- |
| aug-cc-pvtz | 1e-14 | 3 | 22311 | 157.70 | 28.54 | 5.5 | -- |
| aug-cc-pvtz | 1e-12 | 3 | 22311 | 157.70 | 24.23 | 6.5 | -- |
| aug-cc-pvqz | 1e-14 | 4 | 40674 | -- | 101.11 | -- | -- |
| aug-cc-pvqz | 1e-12 | 4 | 40674 | -- | 88.17 | -- | -- |

#### ubiquitin

| basis | threshold | lmax | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 2 | 11577 | 48.99 | 4.34 | 11.3 | -- |
| def2-svp | 1e-12 | 2 | 11577 | 48.99 | 4.02 | 12.2 | -- |
| def2-svpd | 1e-14 | 2 | 17433 | 93.69 | 11.92 | 7.9 | -- |
| def2-svpd | 1e-12 | 2 | 17433 | 93.69 | 9.90 | 9.5 | -- |
| def2-tzvp | 1e-14 | 3 | 22442 | 168.41 | 11.77 | 14.3 | -- |
| def2-tzvp | 1e-12 | 3 | 22442 | 168.41 | 10.30 | 16.4 | -- |
| def2-tzvpp | 1e-14 | 3 | 27479 | 232.85 | 13.33 | 17.5 | -- |
| def2-tzvpp | 1e-12 | 3 | 27479 | 232.85 | 12.04 | 19.3 | -- |
| def2-tzvpd | 1e-14 | 3 | 28298 | 237.12 | 30.83 | 7.7 | -- |
| def2-tzvpd | 1e-12 | 3 | 28298 | 237.12 | 27.45 | 8.6 | -- |
| def2-tzvppd | 1e-14 | 3 | 33335 | -- | 35.11 | -- | -- |
| def2-tzvppd | 1e-12 | 3 | 33335 | -- | 31.02 | -- | -- |
| def2-qzvp | 1e-14 | 4 | 53197 | -- | 51.87 | -- | -- |
| def2-qzvp | 1e-12 | 4 | 53197 | -- | 43.43 | -- | -- |
| def2-qzvpp | 1e-14 | 4 | 53197 | -- | 52.79 | -- | -- |
| def2-qzvpp | 1e-12 | 4 | 53197 | -- | 44.14 | -- | -- |
| def2-qzvpd | 1e-14 | 4 | 59053 | -- | 115.53 | -- | -- |
| def2-qzvpd | 1e-12 | 4 | 59053 | -- | 98.85 | -- | -- |
| def2-qzvppd | 1e-14 | 4 | 59053 | -- | 115.37 | -- | -- |
| def2-qzvppd | 1e-12 | 4 | 59053 | -- | 96.26 | -- | -- |
| cc-pvdz | 1e-14 | 2 | 11577 | 76.37 | 5.43 | 14.1 | -- |
| cc-pvdz | 1e-12 | 2 | 11577 | 76.37 | 4.74 | 16.1 | -- |
| cc-pvtz | 1e-14 | 3 | 26870 | 303.21 | 13.36 | 22.7 | -- |
| cc-pvtz | 1e-12 | 3 | 26870 | 303.21 | 11.48 | 26.4 | -- |
| cc-pvqz | 1e-14 | 4 | 51984 | -- | 44.18 | -- | -- |
| cc-pvqz | 1e-12 | 4 | 51984 | -- | 38.17 | -- | -- |
| aug-cc-pvdz | 1e-14 | 2 | 19511 | 147.79 | 21.69 | 6.8 | -- |
| aug-cc-pvdz | 1e-12 | 2 | 19511 | 147.79 | 18.26 | 8.1 | -- |
| aug-cc-pvtz | 1e-14 | 3 | 42163 | -- | 79.26 | -- | -- |
| aug-cc-pvtz | 1e-12 | 3 | 42163 | -- | 68.13 | -- | -- |
| aug-cc-pvqz | 1e-14 | 4 | 77098 | -- | 335.63 | -- | -- |
| aug-cc-pvqz | 1e-12 | 4 | 77098 | -- | 286.71 | -- | -- |

The geometric mean over the 170 comparable cases is **3.03**, against 2.51 before.
The agreement holds at 1e-14 absolute throughout, which is where a driver built on
the same recurrences should sit.

### The two-center Coulomb driver, the def2 fitting sets

| molecule | basis | nao | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- |
| tagrisso | jfit | 2176 | 2.44 | 3.02 | 0.8 | 3.48e-13 |
| tagrisso | jkfit | 3387 | 4.14 | 4.71 | 0.9 | 5.12e-13 |
| c60 | jfit | 2940 | 4.02 | 4.32 | 0.9 | 1.28e-13 |
| c60 | jkfit | 4500 | 6.15 | 7.07 | 0.9 | 5.51e-13 |
| taxol | jfit | 3528 | 6.82 | 5.19 | 1.3 | 4.12e-13 |
| taxol | jkfit | 5489 | 8.72 | 8.96 | 1.0 | 5.97e-13 |
| paracetamol_cluster | jfit | 10208 | 35.37 | 39.12 | 0.9 | 4.05e-13 |
| paracetamol_cluster | jkfit | 15888 | 76.42 | 72.43 | 1.1 | 4.97e-13 |
| crambin | jfit | 19500 | 135.20 | 141.94 | 1.0 | 3.55e-13 |
| crambin | jkfit | 30751 | 613.97 | 290.86 | 2.1 | -- |
| ubiquitin | jfit | 36419 | 882.18 | 569.70 | 1.5 | -- |
| ubiquitin | jkfit | 56971 | -- | 1112.72 | -- | -- |

The geometric mean is **1.07**, against 0.84 before. It crosses over with size:
0.84 on tagrisso, 1.42 on crambin, 1.55 on ubiquitin. Nothing screens here, so both
sides compute every atom pair and the comparison is arithmetic against arithmetic.
What moves the large cases is the packed storage as much as the kernels.

### The two-center Coulomb driver, the RI fitting sets

| molecule | basis | lmax | nao | GB | ref ms | simd ms | x ref | max abs diff |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tagrisso | cc-pvdz | 3 | 2534 | 0.02 | 2.16 | 2.28 | 0.9 | 2.03e-13 |
| tagrisso | cc-pvtz | 4 | 3987 | 0.06 | 7.03 | 4.76 | 1.5 | 2.98e-13 |
| tagrisso | cc-pvqz | 5 | 6699 | 0.17 | 15.66 | 16.85 | 0.9 | 2.06e-13 |
| tagrisso | cc-pv5z | 6 | 10144 | 0.38 | 42.39 | 60.38 | 0.7 | 4.83e-13 |
| tagrisso | cc-pv6z | 7 | 15091 | 0.85 | -- | 260.07 | -- | -- |
| tagrisso | aug-cc-pvdz | 3 | 3423 | 0.04 | 10.73 | 3.45 | 3.1 | 1.34e-12 |
| tagrisso | aug-cc-pvtz | 4 | 5440 | 0.11 | 11.50 | 8.70 | 1.3 | 8.81e-13 |
| tagrisso | aug-cc-pvqz | 5 | 8856 | 0.29 | 29.35 | 36.85 | 0.8 | 9.52e-13 |
| tagrisso | aug-cc-pv5z | 6 | 13145 | 0.64 | 81.55 | 127.40 | 0.6 | 1.22e-12 |
| tagrisso | aug-cc-pv6z | 7 | 19076 | 1.36 | -- | 500.53 | -- | -- |
| c60 | cc-pvdz | 3 | 3360 | 0.04 | 2.71 | 3.60 | 0.8 | 2.10e-13 |
| c60 | cc-pvtz | 4 | 4860 | 0.09 | 5.95 | 8.11 | 0.7 | 2.42e-13 |
| c60 | cc-pvqz | 5 | 7920 | 0.23 | 21.05 | 24.50 | 0.9 | 2.86e-13 |
| c60 | cc-pv5z | 6 | 11580 | 0.50 | 54.95 | 80.65 | 0.7 | 3.13e-13 |
| c60 | cc-pv6z | 7 | 16980 | 1.07 | -- | 356.43 | -- | -- |
| c60 | aug-cc-pvdz | 3 | 4320 | 0.07 | 4.49 | 5.72 | 0.8 | 9.38e-13 |
| c60 | aug-cc-pvtz | 4 | 6360 | 0.15 | 10.61 | 14.38 | 0.7 | 9.95e-13 |
| c60 | aug-cc-pvqz | 5 | 10080 | 0.38 | 35.25 | 42.43 | 0.8 | 6.25e-13 |
| c60 | aug-cc-pv5z | 6 | 14520 | 0.79 | 96.30 | 158.16 | 0.6 | 1.19e-12 |
| c60 | aug-cc-pv6z | 7 | 20820 | 1.61 | -- | 629.01 | -- | -- |
| taxol | cc-pvdz | 3 | 4102 | 0.06 | 4.25 | 4.20 | 1.0 | 2.20e-13 |
| taxol | cc-pvtz | 4 | 6411 | 0.15 | 10.17 | 10.34 | 1.0 | 2.20e-13 |
| taxol | cc-pvqz | 5 | 10747 | 0.43 | 35.46 | 40.39 | 0.9 | 2.33e-13 |
| taxol | cc-pv5z | 6 | 16232 | 0.98 | 103.13 | 142.87 | 0.7 | 6.25e-13 |
| taxol | cc-pv6z | 7 | 24123 | 2.17 | -- | 595.88 | -- | -- |
| taxol | aug-cc-pvdz | 3 | 5519 | 0.11 | 7.37 | 6.83 | 1.1 | 1.19e-12 |
| taxol | aug-cc-pvtz | 4 | 8720 | 0.28 | 19.34 | 18.71 | 1.0 | 8.95e-13 |
| taxol | aug-cc-pvqz | 5 | 14168 | 0.75 | 64.11 | 79.42 | 0.8 | 1.02e-12 |
| taxol | aug-cc-pv5z | 6 | 20985 | 1.64 | 186.79 | 300.83 | 0.6 | -- |
| taxol | aug-cc-pv6z | 7 | 30428 | 3.45 | -- | 1221.81 | -- | -- |
| Cu_PPh3_4_cation | cc-pvtz | 6 | 8364 | 0.26 | 20.24 | 17.43 | 1.2 | 1.92e-13 |
| Cu_PPh3_4_cation | cc-pvqz | 7 | 13798 | 0.71 | -- | 66.40 | -- | -- |
| Cu_PPh3_4_cation | cc-pv5z | 8 | 20756 | 1.60 | -- | 253.60 | -- | -- |
| Cu_PPh3_4_cation | aug-cc-pvtz | 6 | 11273 | 0.47 | 38.91 | 35.10 | 1.1 | 8.81e-13 |
| Cu_PPh3_4_cation | aug-cc-pvqz | 7 | 18098 | 1.22 | -- | 132.02 | -- | -- |
| Cu_PPh3_4_cation | aug-cc-pv5z | 8 | 26721 | 2.66 | -- | 510.58 | -- | -- |
| paracetamol_cluster | cc-pvdz | 3 | 11872 | 0.53 | 31.88 | 32.77 | 1.0 | 1.78e-13 |
| paracetamol_cluster | cc-pvtz | 4 | 18576 | 1.29 | 88.23 | 86.73 | 1.0 | 2.49e-13 |
| paracetamol_cluster | cc-pvqz | 5 | 31152 | 3.62 | 621.83 | 348.95 | 1.8 | -- |
| paracetamol_cluster | cc-pv5z | 6 | 47072 | 8.25 | -- | 1377.96 | -- | -- |
| paracetamol_cluster | aug-cc-pvdz | 3 | 15984 | 0.95 | 58.59 | 57.07 | 1.0 | 7.67e-13 |
| paracetamol_cluster | aug-cc-pvtz | 4 | 25280 | 2.38 | 166.40 | 177.81 | 0.9 | -- |
| paracetamol_cluster | aug-cc-pvqz | 5 | 41088 | 6.29 | -- | 778.62 | -- | -- |
| paracetamol_cluster | aug-cc-pv5z | 6 | 60880 | 13.81 | -- | 2874.28 | -- | -- |
| crambin | cc-pvdz | 3 | 22842 | 1.94 | 119.15 | 118.70 | 1.0 | -- |
| crambin | cc-pvtz | 4 | 36183 | 4.88 | 693.31 | 387.52 | 1.8 | -- |
| crambin | cc-pvqz | 5 | 60645 | 13.70 | -- | 1507.05 | -- | -- |
| crambin | aug-cc-pvdz | 3 | 30909 | 3.56 | 215.39 | 217.00 | 1.0 | -- |
| crambin | aug-cc-pvtz | 4 | 49398 | 9.09 | -- | 760.04 | -- | -- |
| ubiquitin | cc-pvdz | 3 | 42538 | 6.74 | -- | 511.81 | -- | -- |
| ubiquitin | aug-cc-pvdz | 3 | 57831 | 12.46 | -- | 922.74 | -- | -- |

The geometric mean over the 34 cases with a reference is **0.96**, against 0.72
before. This is still the one place the reference is ahead, and it is ahead in the
middle of the range: c60 in aug-cc-pV5Z-RIFIT at 0.6, taxol in cc-pV5Z-RIFIT at
0.7. It falls behind at both ends -- tagrisso in aug-cc-pVDZ-RIFIT at 3.1,
paracetamol in cc-pVQZ-RIFIT at 1.8, crambin in cc-pVTZ-RIFIT at 1.8.

17 of the 51 cases have no reference: everything reaching k or l, which the
reference dispatcher does not implement. That includes all six cases of the copper
complex these kernels were extended for, which now run at 66.4 ms for 13798
functions reaching k and 253.6 ms for 20756 reaching l.

### What these numbers say

The two changes moved every driver by about the same fifth, which is what a change
to shared plumbing should do. Where they leave the three is not the same place.

The overlap and the kinetic energy drivers are far ahead of the reference and pull
further ahead with the molecule, because they skip atom pairs the reference
computes. Their ratios rose from 2.86 and 2.51 to 3.11 and 3.03 without any change
to the integrals themselves.

The two-center Coulomb driver has no screening to win with, and its ratio moved
from 0.84 to 1.07 on the def2 fitting sets and from 0.72 to 0.96 on the RI fitting
sets. It has drawn level, and it reaches angular momenta the reference cannot.

The overlap driver is now transform bound: 57 per cent of its working samples are
in `simdtrf::`, spread over some twenty horizontal transfer and spherical transform
routines with no single one above 1.9 per cent. There is no hot spot left there,
only volume.

## The unit of work the threads draw on

Three changes, in the order they were made and measured.

The **unit of work is now one combination of basis functions of one block** and no
longer one block. A block is what the sparsity and the storage are described in and it
carries a fixed cost, so it cannot be made small enough to feed a large machine: an
ordinary molecule holds tens of blocks whatever the number of the threads, as the
target size of a block is bounded from below. c60 forms seven blocks at any core
count, taxol twenty eight, the copper complex forty one. The combinations of a block
are tens to hundreds, cost nothing to enumerate, and are independent of one another,
so a flat loop over them feeds a machine the loop over blocks cannot. Their ceiling,
the total work divided by the cost of the fattest task, is 67 to 142 for a small
molecule in an ordinary basis and two to fourteen thousand for the fitting sets.

A **loop of too few iterations is left to the encountering thread**. Opening a region
costs a fork and a join whose price grows with the number of the threads while the work
of a loop with a handful of iterations does not. The chunk loop of make_pair_groups
runs a single chunk for c60 and cost 0.008 ms on one thread against 0.060 ms on
sixteen, for the same work. The bound applies only to a loop whose iteration carries
little work: gating the ordering and the screening of the blocks the same way was
measured and made c60 half again slower, as seven blocks over seven threads beats seven
blocks over one by far more than a fork costs.

The **keys of a small group are ordered by comparison and not by radix**. This was the
largest of the three by a wide margin.

### Eight microseconds a block, and where they were

Sweeping the target size of a block at a fixed number of threads shows the time of a
driver to be linear in the number of the blocks, at about eight microseconds each:
crambin gave 8.4 and taxol 8.2, measured independently. That constant caps the cores a
driver can use at roughly its single thread time divided by forty times it -- ten cores
for taxol, eighty for crambin -- and it is why the smallest useful block is bounded
from below, which is in turn why an ordinary molecule cannot fill a large machine.

Timing the three stages inside sort_by_distance found all of it in one place.

| stage | taxol, 214 pairs a block | crambin, 3170 pairs a block |
| --- | --- | --- |
| _make_keys | 0.00034 ms | 0.00428 ms |
| _order_keys, the radix | 0.0988 ms | 0.1252 ms |
| _move_pairs | 0.00024 ms | 0.00322 ms |

The radix cost two hundred and ninety times what forming the keys cost, and it barely
moved with the data: fifteen times the atom pairs cost 1.27 times the time. It ran four
passes over sixty five thousand counters whatever the number of the keys, so a block of
two hundred atom pairs paid half a million counter operations for two hundred elements
of data. Divided by the threads, that is the eight microseconds the sweep measured.

The keys are now ordered by comparison below four thousand and ninety six atom pairs
and by radix above. The two give the same order: the second half of a key is the index
of the atom pair and increases with it, so ordering the keys as pairs is the stable
radix on their first half. The bound earns its place at both ends. taxol runs 1.084 ms
with the radix always and 0.819 ms with the bound anywhere between 1024 and 8192, while
ubiquitin, whose blocks hold thirteen thousand atom pairs, runs 7.47 ms with the bound
and 8.72 ms with the comparison always.

The fixed cost of a block fell from 8.4 to 2.6 microseconds on crambin and from 8.2 to
4.2 on taxol. crambin at a target size of eight atom pairs went from 214.7 ms to 71.8.

### What was tried first and did not work

Four rounds of work before this one measured as noise, and the reason is the same each
time.

Making the hot assertions lazy removed two heap allocations from every call of
`_cell_index`, which a profile had put at 39.7 per cent of the overlap driver. It moved
the single thread time from 3.32 ms to 3.23. **macOS `sample` over-reports allocator
frames by more than an order of magnitude**, as it catches threads holding the malloc
lock out of proportion, and the ranking it gave was wrong.

Collapsing the parallel regions of make_pattern from three into one bought nothing. The
cost is per thread wake-up, not per region, and a barrier inside one region wakes
sleeping threads exactly as a fork does. Holding the threads awake proves it: at eight
threads `KMP_BLOCKTIME=1` is worth 1.17 times, and past twelve, which is this machine's
performance core count, it inverts violently.

Reusing the sparsity pattern between calls was abandoned before it was written. Most
two-center quantities are computed once per calculation, so there is nothing to reuse
across, and removing work does not make the remaining work scale.

Folding the cost of a block into the description of the block, which is inside a
parallel region, cut the serial pass that recomputed it from 0.140 ms to 0.0047 and
stopped it growing with the thread count. Correct, and worth about one per cent.

What found the radix was none of these. It was sweeping the size of a block and
watching the time move.


### The overlap driver against the reference

`OMP_NUM_THREADS=14`, best of two to five runs. **ref** is `OverlapDriver`, which
computes every atom pair and carries no threshold, so one number serves both threshold
rows. A dash in the ref column marks a case it cannot run: it returns a dense matrix,
and ubiquitin in aug-cc-pV6Z would be three hundred gigabytes.

The geometric mean over the 196 cases which have a reference is **3.89**, against 3.11
before these changes, on the same grid.


#### tagrisso

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 683 | 1.56 | 0.61 | 2.6 |
| def2-svp | 1e-12 | 683 | 1.56 | 0.62 | 2.5 |
| def2-svpd | 1e-14 | 1010 | 1.59 | 0.73 | 2.2 |
| def2-svpd | 1e-12 | 1010 | 1.59 | 0.71 | 2.2 |
| def2-tzvp | 1e-14 | 1345 | 1.09 | 0.72 | 1.5 |
| def2-tzvp | 1e-12 | 1345 | 1.09 | 0.66 | 1.7 |
| def2-tzvpp | 1e-14 | 1609 | 1.00 | 0.74 | 1.4 |
| def2-tzvpp | 1e-12 | 1609 | 1.00 | 0.73 | 1.4 |
| def2-tzvpd | 1e-14 | 1672 | 1.79 | 0.81 | 2.2 |
| def2-tzvpd | 1e-12 | 1672 | 1.79 | 1.04 | 1.7 |
| def2-tzvppd | 1e-14 | 1936 | 1.92 | 1.15 | 1.7 |
| def2-tzvppd | 1e-12 | 1936 | 1.92 | 1.11 | 1.7 |
| def2-qzvp | 1e-14 | 3099 | 3.18 | 1.24 | 2.6 |
| def2-qzvp | 1e-12 | 3099 | 3.18 | 1.17 | 2.7 |
| def2-qzvpp | 1e-14 | 3099 | 3.05 | 1.25 | 2.4 |
| def2-qzvpp | 1e-12 | 3099 | 3.05 | 1.23 | 2.5 |
| def2-qzvpd | 1e-14 | 3426 | 3.84 | 1.49 | 2.6 |
| def2-qzvpd | 1e-12 | 3426 | 3.84 | 1.46 | 2.6 |
| def2-qzvppd | 1e-14 | 3426 | 3.21 | 1.49 | 2.2 |
| def2-qzvppd | 1e-12 | 3426 | 3.21 | 1.45 | 2.2 |
| cc-pvdz | 1e-14 | 683 | 1.59 | 0.62 | 2.6 |
| cc-pvdz | 1e-12 | 683 | 1.59 | 0.63 | 2.5 |
| cc-pvtz | 1e-14 | 1572 | 1.71 | 0.92 | 1.9 |
| cc-pvtz | 1e-12 | 1572 | 1.71 | 0.90 | 1.9 |
| cc-pvqz | 1e-14 | 3025 | 3.05 | 1.17 | 2.6 |
| cc-pvqz | 1e-12 | 3025 | 3.05 | 1.15 | 2.7 |
| cc-pv5z | 1e-14 | 5182 | 8.76 | 2.55 | 3.4 |
| cc-pv5z | 1e-12 | 5182 | 8.76 | 2.41 | 3.6 |
| cc-pv6z | 1e-14 | 8183 | 23.78 | 6.26 | 3.8 |
| cc-pv6z | 1e-12 | 8183 | 23.78 | 5.94 | 4.0 |
| aug-cc-pvdz | 1e-14 | 1148 | 3.33 | 0.65 | 5.1 |
| aug-cc-pvdz | 1e-12 | 1148 | 3.33 | 0.65 | 5.1 |
| aug-cc-pvtz | 1e-14 | 2461 | 2.06 | 1.17 | 1.8 |
| aug-cc-pvtz | 1e-12 | 2461 | 2.06 | 1.50 | 1.4 |
| aug-cc-pvqz | 1e-14 | 4478 | 5.71 | 2.39 | 2.4 |
| aug-cc-pvqz | 1e-12 | 4478 | 5.71 | 2.30 | 2.5 |
| aug-cc-pv5z | 1e-14 | 7339 | 16.84 | 5.93 | 2.8 |
| aug-cc-pv5z | 1e-12 | 7339 | 16.84 | 5.64 | 3.0 |
| aug-cc-pv6z | 1e-14 | 11184 | 45.59 | 15.45 | 3.0 |
| aug-cc-pv6z | 1e-12 | 11184 | 45.59 | 14.69 | 3.1 |

#### c60

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 840 | 0.45 | 0.45 | 1.0 |
| def2-svp | 1e-12 | 840 | 0.45 | 0.46 | 1.0 |
| def2-svpd | 1e-14 | 1200 | 1.75 | 0.69 | 2.5 |
| def2-svpd | 1e-12 | 1200 | 1.75 | 0.69 | 2.5 |
| def2-tzvp | 1e-14 | 1860 | 1.18 | 0.80 | 1.5 |
| def2-tzvp | 1e-12 | 1860 | 1.18 | 1.01 | 1.2 |
| def2-tzvpp | 1e-14 | 1860 | 1.63 | 0.78 | 2.1 |
| def2-tzvpp | 1e-12 | 1860 | 1.63 | 1.01 | 1.6 |
| def2-tzvpd | 1e-14 | 2220 | 1.60 | 0.92 | 1.7 |
| def2-tzvpd | 1e-12 | 2220 | 1.60 | 0.94 | 1.7 |
| def2-tzvppd | 1e-14 | 2220 | 1.85 | 0.92 | 2.0 |
| def2-tzvppd | 1e-12 | 2220 | 1.85 | 0.91 | 2.0 |
| def2-qzvp | 1e-14 | 3420 | 3.43 | 1.64 | 2.1 |
| def2-qzvp | 1e-12 | 3420 | 3.43 | 1.53 | 2.2 |
| def2-qzvpp | 1e-14 | 3420 | 3.80 | 1.61 | 2.4 |
| def2-qzvpp | 1e-12 | 3420 | 3.80 | 1.52 | 2.5 |
| def2-qzvpd | 1e-14 | 3780 | 3.87 | 1.88 | 2.1 |
| def2-qzvpd | 1e-12 | 3780 | 3.87 | 1.81 | 2.1 |
| def2-qzvppd | 1e-14 | 3780 | 4.84 | 1.91 | 2.5 |
| def2-qzvppd | 1e-12 | 3780 | 4.84 | 1.83 | 2.6 |
| cc-pvdz | 1e-14 | 840 | 1.62 | 0.67 | 2.4 |
| cc-pvdz | 1e-12 | 840 | 1.62 | 0.65 | 2.5 |
| cc-pvtz | 1e-14 | 1800 | 1.52 | 0.81 | 1.9 |
| cc-pvtz | 1e-12 | 1800 | 1.52 | 0.75 | 2.0 |
| cc-pvqz | 1e-14 | 3300 | 3.60 | 1.57 | 2.3 |
| cc-pvqz | 1e-12 | 3300 | 3.60 | 1.53 | 2.4 |
| cc-pv5z | 1e-14 | 5460 | 9.41 | 3.81 | 2.5 |
| cc-pv5z | 1e-12 | 5460 | 9.41 | 3.57 | 2.6 |
| cc-pv6z | 1e-14 | 8400 | 26.02 | 10.04 | 2.6 |
| cc-pv6z | 1e-12 | 8400 | 26.02 | 9.41 | 2.8 |
| aug-cc-pvdz | 1e-14 | 1380 | 1.24 | 0.69 | 1.8 |
| aug-cc-pvdz | 1e-12 | 1380 | 1.24 | 0.68 | 1.8 |
| aug-cc-pvtz | 1e-14 | 2760 | 3.37 | 1.34 | 2.5 |
| aug-cc-pvtz | 1e-12 | 2760 | 3.37 | 1.91 | 1.8 |
| aug-cc-pvqz | 1e-14 | 4800 | 6.73 | 3.27 | 2.1 |
| aug-cc-pvqz | 1e-12 | 4800 | 6.73 | 3.18 | 2.1 |
| aug-cc-pv5z | 1e-14 | 7620 | 19.01 | 9.04 | 2.1 |
| aug-cc-pv5z | 1e-12 | 7620 | 19.01 | 8.75 | 2.2 |
| aug-cc-pv6z | 1e-14 | 11340 | 51.09 | 25.65 | 2.0 |
| aug-cc-pv6z | 1e-12 | 11340 | 51.09 | 24.37 | 2.1 |

#### taxol

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 1099 | 0.95 | 0.61 | 1.6 |
| def2-svp | 1e-12 | 1099 | 0.95 | 0.59 | 1.6 |
| def2-svpd | 1e-14 | 1657 | 1.86 | 0.99 | 1.9 |
| def2-svpd | 1e-12 | 1657 | 1.86 | 0.96 | 1.9 |
| def2-tzvp | 1e-14 | 2185 | 2.12 | 0.89 | 2.4 |
| def2-tzvp | 1e-12 | 2185 | 2.12 | 1.04 | 2.0 |
| def2-tzvpp | 1e-14 | 2577 | 2.53 | 0.97 | 2.6 |
| def2-tzvpp | 1e-12 | 2577 | 2.53 | 1.27 | 2.0 |
| def2-tzvpd | 1e-14 | 2743 | 2.67 | 1.66 | 1.6 |
| def2-tzvpd | 1e-12 | 2743 | 2.67 | 1.59 | 1.7 |
| def2-tzvppd | 1e-14 | 3135 | 2.50 | 1.30 | 1.9 |
| def2-tzvppd | 1e-12 | 3135 | 2.50 | 1.78 | 1.4 |
| def2-qzvp | 1e-14 | 4947 | 6.18 | 2.00 | 3.1 |
| def2-qzvp | 1e-12 | 4947 | 6.18 | 1.93 | 3.2 |
| def2-qzvpp | 1e-14 | 4947 | 6.36 | 1.97 | 3.2 |
| def2-qzvpp | 1e-12 | 4947 | 6.36 | 1.92 | 3.3 |
| def2-qzvpd | 1e-14 | 5505 | 8.34 | 2.60 | 3.2 |
| def2-qzvpd | 1e-12 | 5505 | 8.34 | 2.49 | 3.3 |
| def2-qzvppd | 1e-14 | 5505 | 12.27 | 2.78 | 4.4 |
| def2-qzvppd | 1e-12 | 5505 | 12.27 | 2.58 | 4.8 |
| cc-pvdz | 1e-14 | 1099 | 1.29 | 0.66 | 2.0 |
| cc-pvdz | 1e-12 | 1099 | 1.29 | 0.65 | 2.0 |
| cc-pvtz | 1e-14 | 2516 | 2.48 | 0.99 | 2.5 |
| cc-pvtz | 1e-12 | 2516 | 2.48 | 1.28 | 1.9 |
| cc-pvqz | 1e-14 | 4825 | 6.98 | 1.90 | 3.7 |
| cc-pvqz | 1e-12 | 4825 | 6.98 | 1.82 | 3.8 |
| cc-pv5z | 1e-14 | 8246 | 21.64 | 4.59 | 4.7 |
| cc-pv5z | 1e-12 | 8246 | 21.64 | 4.27 | 5.1 |
| cc-pv6z | 1e-14 | 12999 | 57.31 | 11.99 | 4.8 |
| cc-pv6z | 1e-12 | 12999 | 57.31 | 11.31 | 5.1 |
| aug-cc-pvdz | 1e-14 | 1844 | 1.47 | 0.93 | 1.6 |
| aug-cc-pvdz | 1e-12 | 1844 | 1.47 | 0.90 | 1.6 |
| aug-cc-pvtz | 1e-14 | 3933 | 3.98 | 1.84 | 2.2 |
| aug-cc-pvtz | 1e-12 | 3933 | 3.98 | 1.78 | 2.2 |
| aug-cc-pvqz | 1e-14 | 7134 | 12.72 | 4.45 | 2.9 |
| aug-cc-pvqz | 1e-12 | 7134 | 12.72 | 4.18 | 3.0 |
| aug-cc-pv5z | 1e-14 | 11667 | 41.72 | 11.84 | 3.5 |
| aug-cc-pv5z | 1e-12 | 11667 | 41.72 | 11.11 | 3.8 |
| aug-cc-pv6z | 1e-14 | 17752 | 109.99 | 31.93 | 3.4 |
| aug-cc-pv6z | 1e-12 | 17752 | 109.99 | 30.16 | 3.6 |

#### paracetamol_cluster

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 3184 | 3.30 | 1.13 | 2.9 |
| def2-svp | 1e-12 | 3184 | 3.30 | 1.04 | 3.2 |
| def2-svpd | 1e-14 | 4768 | 6.21 | 1.94 | 3.2 |
| def2-svpd | 1e-12 | 4768 | 6.21 | 1.79 | 3.5 |
| def2-tzvp | 1e-14 | 6320 | 10.77 | 2.09 | 5.2 |
| def2-tzvp | 1e-12 | 6320 | 10.77 | 1.94 | 5.6 |
| def2-tzvpp | 1e-14 | 7472 | 13.96 | 2.44 | 5.7 |
| def2-tzvpp | 1e-12 | 7472 | 13.96 | 2.23 | 6.3 |
| def2-tzvpd | 1e-14 | 7904 | 15.66 | 3.66 | 4.3 |
| def2-tzvpd | 1e-12 | 7904 | 15.66 | 3.40 | 4.6 |
| def2-tzvppd | 1e-14 | 9056 | 19.70 | 4.21 | 4.7 |
| def2-tzvppd | 1e-12 | 9056 | 19.70 | 3.91 | 5.0 |
| def2-qzvp | 1e-14 | 14352 | 51.03 | 6.19 | 8.2 |
| def2-qzvp | 1e-12 | 14352 | 51.03 | 5.70 | 9.0 |
| def2-qzvpp | 1e-14 | 14352 | 50.86 | 6.16 | 8.3 |
| def2-qzvpp | 1e-12 | 14352 | 50.86 | 5.71 | 8.9 |
| def2-qzvpd | 1e-14 | 15936 | 64.12 | 9.74 | 6.6 |
| def2-qzvpd | 1e-12 | 15936 | 64.12 | 8.91 | 7.2 |
| def2-qzvppd | 1e-14 | 15936 | 64.27 | 9.91 | 6.5 |
| def2-qzvppd | 1e-12 | 15936 | 64.27 | 9.10 | 7.1 |
| cc-pvdz | 1e-14 | 3184 | 5.41 | 1.26 | 4.3 |
| cc-pvdz | 1e-12 | 3184 | 5.41 | 1.20 | 4.5 |
| cc-pvtz | 1e-14 | 7296 | 14.87 | 2.42 | 6.1 |
| cc-pvtz | 1e-12 | 7296 | 14.87 | 2.28 | 6.5 |
| cc-pvqz | 1e-14 | 14000 | 51.09 | 5.64 | 9.1 |
| cc-pvqz | 1e-12 | 14000 | 51.09 | 5.17 | 9.9 |
| cc-pv5z | 1e-14 | 23936 | 164.73 | 15.58 | 10.6 |
| cc-pv5z | 1e-12 | 23936 | 164.73 | 14.08 | 11.7 |
| cc-pv6z | 1e-14 | 37744 | -- | 42.30 | -- |
| cc-pv6z | 1e-12 | 37744 | -- | 38.70 | -- |
| aug-cc-pvdz | 1e-14 | 5344 | 9.84 | 2.71 | 3.6 |
| aug-cc-pvdz | 1e-12 | 5344 | 9.84 | 2.57 | 3.8 |
| aug-cc-pvtz | 1e-14 | 11408 | 35.36 | 7.11 | 5.0 |
| aug-cc-pvtz | 1e-12 | 11408 | 35.36 | 6.51 | 5.4 |
| aug-cc-pvqz | 1e-14 | 20704 | 116.58 | 19.88 | 5.9 |
| aug-cc-pvqz | 1e-12 | 20704 | 116.58 | 18.16 | 6.4 |
| aug-cc-pv5z | 1e-14 | 33872 | -- | 54.78 | -- |
| aug-cc-pv5z | 1e-12 | 33872 | -- | 50.41 | -- |
| aug-cc-pv6z | 1e-14 | 51552 | -- | 159.66 | -- |
| aug-cc-pv6z | 1e-12 | 51552 | -- | 136.13 | -- |

#### crambin

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 6177 | 12.19 | 2.03 | 6.0 |
| def2-svp | 1e-12 | 6177 | 12.19 | 1.92 | 6.3 |
| def2-svpd | 1e-14 | 9294 | 23.20 | 3.99 | 5.8 |
| def2-svpd | 1e-12 | 9294 | 23.20 | 3.63 | 6.4 |
| def2-tzvp | 1e-14 | 12063 | 45.55 | 4.01 | 11.4 |
| def2-tzvp | 1e-12 | 12063 | 45.55 | 3.65 | 12.5 |
| def2-tzvpp | 1e-14 | 14613 | 56.24 | 4.69 | 12.0 |
| def2-tzvpp | 1e-12 | 14613 | 56.24 | 4.25 | 13.2 |
| def2-tzvpd | 1e-14 | 15180 | 61.28 | 7.98 | 7.7 |
| def2-tzvpd | 1e-12 | 15180 | 61.28 | 7.30 | 8.4 |
| def2-tzvppd | 1e-14 | 17730 | 79.75 | 9.40 | 8.5 |
| def2-tzvppd | 1e-12 | 17730 | 79.75 | 8.47 | 9.4 |
| def2-qzvp | 1e-14 | 28167 | 210.30 | 12.91 | 16.3 |
| def2-qzvp | 1e-12 | 28167 | 210.30 | 11.39 | 18.5 |
| def2-qzvpp | 1e-14 | 28167 | 204.83 | 12.86 | 15.9 |
| def2-qzvpp | 1e-12 | 28167 | 204.83 | 11.35 | 18.0 |
| def2-qzvpd | 1e-14 | 31284 | -- | 22.73 | -- |
| def2-qzvpd | 1e-12 | 31284 | -- | 20.27 | -- |
| def2-qzvppd | 1e-14 | 31284 | -- | 22.64 | -- |
| def2-qzvppd | 1e-12 | 31284 | -- | 20.28 | -- |
| cc-pvdz | 1e-14 | 6177 | 18.82 | 2.36 | 8.0 |
| cc-pvdz | 1e-12 | 6177 | 18.82 | 2.21 | 8.5 |
| cc-pvtz | 1e-14 | 14244 | 59.96 | 4.63 | 13.0 |
| cc-pvtz | 1e-12 | 14244 | 59.96 | 4.22 | 14.2 |
| cc-pvqz | 1e-14 | 27459 | 203.96 | 11.42 | 17.9 |
| cc-pvqz | 1e-12 | 27459 | 203.96 | 10.23 | 19.9 |
| cc-pv5z | 1e-14 | 47106 | -- | 33.35 | -- |
| cc-pv5z | 1e-12 | 47106 | -- | 29.51 | -- |
| cc-pv6z | 1e-14 | 74469 | -- | 95.06 | -- |
| cc-pv6z | 1e-12 | 74469 | -- | 83.53 | -- |
| aug-cc-pvdz | 1e-14 | 10380 | 37.52 | 6.11 | 6.1 |
| aug-cc-pvdz | 1e-12 | 10380 | 37.52 | 5.63 | 6.7 |
| aug-cc-pvtz | 1e-14 | 22311 | 134.73 | 17.08 | 7.9 |
| aug-cc-pvtz | 1e-12 | 22311 | 134.73 | 15.74 | 8.6 |
| aug-cc-pvqz | 1e-14 | 40674 | -- | 49.63 | -- |
| aug-cc-pvqz | 1e-12 | 40674 | -- | 44.19 | -- |
| aug-cc-pv5z | 1e-14 | 66753 | -- | 160.71 | -- |
| aug-cc-pv5z | 1e-12 | 66753 | -- | 131.36 | -- |
| aug-cc-pv6z | 1e-14 | 101832 | -- | 486.54 | -- |
| aug-cc-pv6z | 1e-12 | 101832 | -- | 393.91 | -- |

#### ubiquitin

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 11577 | 41.95 | 3.65 | 11.5 |
| def2-svp | 1e-12 | 11577 | 41.95 | 3.45 | 12.2 |
| def2-svpd | 1e-14 | 17433 | 78.24 | 8.28 | 9.4 |
| def2-svpd | 1e-12 | 17433 | 78.24 | 7.39 | 10.6 |
| def2-tzvp | 1e-14 | 22442 | 137.59 | 7.65 | 18.0 |
| def2-tzvp | 1e-12 | 22442 | 137.59 | 6.72 | 20.5 |
| def2-tzvpp | 1e-14 | 27479 | 183.59 | 8.94 | 20.5 |
| def2-tzvpp | 1e-12 | 27479 | 183.59 | 8.06 | 22.8 |
| def2-tzvpd | 1e-14 | 28298 | 201.79 | 17.68 | 11.4 |
| def2-tzvpd | 1e-12 | 28298 | 201.79 | 15.42 | 13.1 |
| def2-tzvppd | 1e-14 | 33335 | -- | 21.04 | -- |
| def2-tzvppd | 1e-12 | 33335 | -- | 18.31 | -- |
| def2-qzvp | 1e-14 | 53197 | -- | 26.03 | -- |
| def2-qzvp | 1e-12 | 53197 | -- | 22.70 | -- |
| def2-qzvpp | 1e-14 | 53197 | -- | 26.05 | -- |
| def2-qzvpp | 1e-12 | 53197 | -- | 22.46 | -- |
| def2-qzvpd | 1e-14 | 59053 | -- | 55.36 | -- |
| def2-qzvpd | 1e-12 | 59053 | -- | 47.36 | -- |
| def2-qzvppd | 1e-14 | 59053 | -- | 54.81 | -- |
| def2-qzvppd | 1e-12 | 59053 | -- | 47.11 | -- |
| cc-pvdz | 1e-14 | 11577 | 65.24 | 4.16 | 15.7 |
| cc-pvdz | 1e-12 | 11577 | 65.24 | 3.88 | 16.8 |
| cc-pvtz | 1e-14 | 26870 | 218.55 | 8.71 | 25.1 |
| cc-pvtz | 1e-12 | 26870 | 218.55 | 7.60 | 28.8 |
| cc-pvqz | 1e-14 | 51984 | -- | 22.62 | -- |
| cc-pvqz | 1e-12 | 51984 | -- | 19.81 | -- |
| cc-pv5z | 1e-14 | 89381 | -- | 73.71 | -- |
| cc-pv5z | 1e-12 | 89381 | -- | 62.30 | -- |
| cc-pv6z | 1e-14 | 141523 | -- | 245.58 | -- |
| cc-pv6z | 1e-12 | 141523 | -- | 192.29 | -- |
| aug-cc-pvdz | 1e-14 | 19511 | 143.70 | 14.10 | 10.2 |
| aug-cc-pvdz | 1e-12 | 19511 | 143.70 | 12.60 | 11.4 |
| aug-cc-pvtz | 1e-14 | 42163 | -- | 43.20 | -- |
| aug-cc-pvtz | 1e-12 | 42163 | -- | 37.46 | -- |
| aug-cc-pvqz | 1e-14 | 77098 | -- | 163.80 | -- |
| aug-cc-pvqz | 1e-12 | 77098 | -- | 127.99 | -- |
| aug-cc-pv5z | 1e-14 | 126778 | -- | 555.45 | -- |
| aug-cc-pv5z | 1e-12 | 126778 | -- | 432.38 | -- |
| aug-cc-pv6z | 1e-14 | 193665 | -- | 1784.11 | -- |
| aug-cc-pv6z | 1e-12 | 193665 | -- | 1267.84 | -- |

### The kinetic energy driver against the reference

The bases are those whose highest angular momentum is g, where the reference dispatcher
stops. The geometric mean over the 170 comparable cases is **3.94**, against 2.51
before these changes, on the same grid.


#### tagrisso

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 683 | 1.20 | 0.55 | 2.2 |
| def2-svp | 1e-12 | 683 | 1.20 | 0.58 | 2.1 |
| def2-svpd | 1e-14 | 1010 | 0.71 | 0.68 | 1.0 |
| def2-svpd | 1e-12 | 1010 | 0.71 | 0.66 | 1.1 |
| def2-tzvp | 1e-14 | 1345 | 2.66 | 0.77 | 3.5 |
| def2-tzvp | 1e-12 | 1345 | 2.66 | 0.74 | 3.6 |
| def2-tzvpp | 1e-14 | 1609 | 1.75 | 0.82 | 2.1 |
| def2-tzvpp | 1e-12 | 1609 | 1.75 | 0.78 | 2.2 |
| def2-tzvpd | 1e-14 | 1672 | 1.31 | 0.96 | 1.4 |
| def2-tzvpd | 1e-12 | 1672 | 1.31 | 0.94 | 1.4 |
| def2-tzvppd | 1e-14 | 1936 | 1.72 | 1.02 | 1.7 |
| def2-tzvppd | 1e-12 | 1936 | 1.72 | 0.99 | 1.7 |
| def2-qzvp | 1e-14 | 3099 | 5.02 | 1.64 | 3.1 |
| def2-qzvp | 1e-12 | 3099 | 5.02 | 1.61 | 3.1 |
| def2-qzvpp | 1e-14 | 3099 | 4.07 | 1.67 | 2.4 |
| def2-qzvpp | 1e-12 | 3099 | 4.07 | 1.57 | 2.6 |
| def2-qzvpd | 1e-14 | 3426 | 8.98 | 1.98 | 4.5 |
| def2-qzvpd | 1e-12 | 3426 | 8.98 | 1.95 | 4.6 |
| def2-qzvppd | 1e-14 | 3426 | 8.18 | 2.01 | 4.1 |
| def2-qzvppd | 1e-12 | 3426 | 8.18 | 1.94 | 4.2 |
| cc-pvdz | 1e-14 | 683 | 1.02 | 0.54 | 1.9 |
| cc-pvdz | 1e-12 | 683 | 1.02 | 0.57 | 1.8 |
| cc-pvtz | 1e-14 | 1572 | 1.77 | 0.80 | 2.2 |
| cc-pvtz | 1e-12 | 1572 | 1.77 | 0.81 | 2.2 |
| cc-pvqz | 1e-14 | 3025 | 6.81 | 1.56 | 4.4 |
| cc-pvqz | 1e-12 | 3025 | 6.81 | 1.51 | 4.5 |
| aug-cc-pvdz | 1e-14 | 1148 | 1.18 | 0.73 | 1.6 |
| aug-cc-pvdz | 1e-12 | 1148 | 1.18 | 0.79 | 1.5 |
| aug-cc-pvtz | 1e-14 | 2461 | 2.84 | 1.38 | 2.1 |
| aug-cc-pvtz | 1e-12 | 2461 | 2.84 | 1.34 | 2.1 |
| aug-cc-pvqz | 1e-14 | 4478 | 8.98 | 3.56 | 2.5 |
| aug-cc-pvqz | 1e-12 | 4478 | 8.98 | 3.44 | 2.6 |

#### c60

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 840 | 0.69 | 0.47 | 1.5 |
| def2-svp | 1e-12 | 840 | 0.69 | 0.45 | 1.5 |
| def2-svpd | 1e-14 | 1200 | 1.04 | 0.61 | 1.7 |
| def2-svpd | 1e-12 | 1200 | 1.04 | 0.60 | 1.7 |
| def2-tzvp | 1e-14 | 1860 | 2.06 | 0.95 | 2.2 |
| def2-tzvp | 1e-12 | 1860 | 2.06 | 0.92 | 2.2 |
| def2-tzvpp | 1e-14 | 1860 | 1.94 | 0.90 | 2.2 |
| def2-tzvpp | 1e-12 | 1860 | 1.94 | 0.94 | 2.1 |
| def2-tzvpd | 1e-14 | 2220 | 2.51 | 1.21 | 2.1 |
| def2-tzvpd | 1e-12 | 2220 | 2.51 | 1.15 | 2.2 |
| def2-tzvppd | 1e-14 | 2220 | 1.82 | 1.22 | 1.5 |
| def2-tzvppd | 1e-12 | 2220 | 1.82 | 1.14 | 1.6 |
| def2-qzvp | 1e-14 | 3420 | 6.77 | 2.48 | 2.7 |
| def2-qzvp | 1e-12 | 3420 | 6.77 | 2.31 | 2.9 |
| def2-qzvpp | 1e-14 | 3420 | 6.10 | 2.45 | 2.5 |
| def2-qzvpp | 1e-12 | 3420 | 6.10 | 2.34 | 2.6 |
| def2-qzvpd | 1e-14 | 3780 | 6.43 | 2.89 | 2.2 |
| def2-qzvpd | 1e-12 | 3780 | 6.43 | 2.79 | 2.3 |
| def2-qzvppd | 1e-14 | 3780 | 7.62 | 2.94 | 2.6 |
| def2-qzvppd | 1e-12 | 3780 | 7.62 | 2.82 | 2.7 |
| cc-pvdz | 1e-14 | 840 | 2.09 | 0.59 | 3.5 |
| cc-pvdz | 1e-12 | 840 | 2.09 | 0.72 | 2.9 |
| cc-pvtz | 1e-14 | 1800 | 2.57 | 1.01 | 2.5 |
| cc-pvtz | 1e-12 | 1800 | 2.57 | 0.92 | 2.8 |
| cc-pvqz | 1e-14 | 3300 | 7.44 | 2.36 | 3.2 |
| cc-pvqz | 1e-12 | 3300 | 7.44 | 2.23 | 3.3 |
| aug-cc-pvdz | 1e-14 | 1380 | 1.44 | 0.79 | 1.8 |
| aug-cc-pvdz | 1e-12 | 1380 | 1.44 | 0.77 | 1.9 |
| aug-cc-pvtz | 1e-14 | 2760 | 4.00 | 1.82 | 2.2 |
| aug-cc-pvtz | 1e-12 | 2760 | 4.00 | 1.77 | 2.3 |
| aug-cc-pvqz | 1e-14 | 4800 | 10.29 | 5.41 | 1.9 |
| aug-cc-pvqz | 1e-12 | 4800 | 10.29 | 5.31 | 1.9 |

#### taxol

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 1099 | 1.09 | 0.61 | 1.8 |
| def2-svp | 1e-12 | 1099 | 1.09 | 0.62 | 1.8 |
| def2-svpd | 1e-14 | 1657 | 1.12 | 1.11 | 1.0 |
| def2-svpd | 1e-12 | 1657 | 1.12 | 1.07 | 1.0 |
| def2-tzvp | 1e-14 | 2185 | 2.45 | 1.04 | 2.4 |
| def2-tzvp | 1e-12 | 2185 | 2.45 | 1.00 | 2.5 |
| def2-tzvpp | 1e-14 | 2577 | 3.04 | 1.13 | 2.7 |
| def2-tzvpp | 1e-12 | 2577 | 3.04 | 1.11 | 2.7 |
| def2-tzvpd | 1e-14 | 2743 | 2.94 | 1.40 | 2.1 |
| def2-tzvpd | 1e-12 | 2743 | 2.94 | 1.42 | 2.1 |
| def2-tzvppd | 1e-14 | 3135 | 3.31 | 1.59 | 2.1 |
| def2-tzvppd | 1e-12 | 3135 | 3.31 | 1.57 | 2.1 |
| def2-qzvp | 1e-14 | 4947 | 9.10 | 2.85 | 3.2 |
| def2-qzvp | 1e-12 | 4947 | 9.10 | 2.66 | 3.4 |
| def2-qzvpp | 1e-14 | 4947 | 8.77 | 2.89 | 3.0 |
| def2-qzvpp | 1e-12 | 4947 | 8.77 | 2.67 | 3.3 |
| def2-qzvpd | 1e-14 | 5505 | 10.00 | 3.72 | 2.7 |
| def2-qzvpd | 1e-12 | 5505 | 10.00 | 3.57 | 2.8 |
| def2-qzvppd | 1e-14 | 5505 | 11.66 | 3.76 | 3.1 |
| def2-qzvppd | 1e-12 | 5505 | 11.66 | 3.53 | 3.3 |
| cc-pvdz | 1e-14 | 1099 | 2.52 | 0.69 | 3.7 |
| cc-pvdz | 1e-12 | 1099 | 2.52 | 0.69 | 3.7 |
| cc-pvtz | 1e-14 | 2516 | 3.09 | 1.16 | 2.7 |
| cc-pvtz | 1e-12 | 2516 | 3.09 | 1.13 | 2.7 |
| cc-pvqz | 1e-14 | 4825 | 8.97 | 2.71 | 3.3 |
| cc-pvqz | 1e-12 | 4825 | 8.97 | 2.50 | 3.6 |
| aug-cc-pvdz | 1e-14 | 1844 | 1.83 | 1.03 | 1.8 |
| aug-cc-pvdz | 1e-12 | 1844 | 1.83 | 1.01 | 1.8 |
| aug-cc-pvtz | 1e-14 | 3933 | 6.03 | 2.40 | 2.5 |
| aug-cc-pvtz | 1e-12 | 3933 | 6.03 | 2.32 | 2.6 |
| aug-cc-pvqz | 1e-14 | 7134 | 18.22 | 6.93 | 2.6 |
| aug-cc-pvqz | 1e-12 | 7134 | 18.22 | 6.67 | 2.7 |

#### paracetamol_cluster

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 3184 | 4.07 | 1.25 | 3.3 |
| def2-svp | 1e-12 | 3184 | 4.07 | 1.15 | 3.5 |
| def2-svpd | 1e-14 | 4768 | 7.39 | 2.44 | 3.0 |
| def2-svpd | 1e-12 | 4768 | 7.39 | 2.22 | 3.3 |
| def2-tzvp | 1e-14 | 6320 | 13.47 | 2.77 | 4.9 |
| def2-tzvp | 1e-12 | 6320 | 13.47 | 2.54 | 5.3 |
| def2-tzvpp | 1e-14 | 7472 | 18.48 | 3.22 | 5.7 |
| def2-tzvpp | 1e-12 | 7472 | 18.48 | 3.00 | 6.2 |
| def2-tzvpd | 1e-14 | 7904 | 20.76 | 5.23 | 4.0 |
| def2-tzvpd | 1e-12 | 7904 | 20.76 | 4.86 | 4.3 |
| def2-tzvppd | 1e-14 | 9056 | 27.82 | 6.01 | 4.6 |
| def2-tzvppd | 1e-12 | 9056 | 27.82 | 5.55 | 5.0 |
| def2-qzvp | 1e-14 | 14352 | 69.87 | 9.86 | 7.1 |
| def2-qzvp | 1e-12 | 14352 | 69.87 | 9.10 | 7.7 |
| def2-qzvpp | 1e-14 | 14352 | 72.36 | 9.92 | 7.3 |
| def2-qzvpp | 1e-12 | 14352 | 72.36 | 8.90 | 8.1 |
| def2-qzvpd | 1e-14 | 15936 | 84.76 | 16.00 | 5.3 |
| def2-qzvpd | 1e-12 | 15936 | 84.76 | 14.78 | 5.7 |
| def2-qzvppd | 1e-14 | 15936 | 85.76 | 15.98 | 5.4 |
| def2-qzvppd | 1e-12 | 15936 | 85.76 | 14.53 | 5.9 |
| cc-pvdz | 1e-14 | 3184 | 7.52 | 1.47 | 5.1 |
| cc-pvdz | 1e-12 | 3184 | 7.52 | 1.41 | 5.3 |
| cc-pvtz | 1e-14 | 7296 | 19.20 | 3.22 | 6.0 |
| cc-pvtz | 1e-12 | 7296 | 19.20 | 2.96 | 6.5 |
| cc-pvqz | 1e-14 | 14000 | 71.48 | 9.01 | 7.9 |
| cc-pvqz | 1e-12 | 14000 | 71.48 | 8.14 | 8.8 |
| aug-cc-pvdz | 1e-14 | 5344 | 11.89 | 3.49 | 3.4 |
| aug-cc-pvdz | 1e-12 | 5344 | 11.89 | 3.28 | 3.6 |
| aug-cc-pvtz | 1e-14 | 11408 | 41.62 | 10.68 | 3.9 |
| aug-cc-pvtz | 1e-12 | 11408 | 41.62 | 9.79 | 4.3 |
| aug-cc-pvqz | 1e-14 | 20704 | 150.52 | 34.92 | 4.3 |
| aug-cc-pvqz | 1e-12 | 20704 | 150.52 | 31.52 | 4.8 |

#### crambin

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 6177 | 17.18 | 2.34 | 7.3 |
| def2-svp | 1e-12 | 6177 | 17.18 | 2.14 | 8.0 |
| def2-svpd | 1e-14 | 9294 | 28.03 | 5.32 | 5.3 |
| def2-svpd | 1e-12 | 9294 | 28.03 | 4.83 | 5.8 |
| def2-tzvp | 1e-14 | 12063 | 53.59 | 5.43 | 9.9 |
| def2-tzvp | 1e-12 | 12063 | 53.59 | 4.97 | 10.8 |
| def2-tzvpp | 1e-14 | 14613 | 71.48 | 6.54 | 10.9 |
| def2-tzvpp | 1e-12 | 14613 | 71.48 | 5.90 | 12.1 |
| def2-tzvpd | 1e-14 | 15180 | 75.51 | 12.18 | 6.2 |
| def2-tzvpd | 1e-12 | 15180 | 75.51 | 10.74 | 7.0 |
| def2-tzvppd | 1e-14 | 17730 | 93.71 | 14.36 | 6.5 |
| def2-tzvppd | 1e-12 | 17730 | 93.71 | 12.59 | 7.4 |
| def2-qzvp | 1e-14 | 28167 | 279.01 | 22.00 | 12.7 |
| def2-qzvp | 1e-12 | 28167 | 279.01 | 19.44 | 14.4 |
| def2-qzvpp | 1e-14 | 28167 | 274.63 | 21.85 | 12.6 |
| def2-qzvpp | 1e-12 | 28167 | 274.63 | 19.14 | 14.3 |
| def2-qzvpd | 1e-14 | 31284 | -- | 40.69 | -- |
| def2-qzvpd | 1e-12 | 31284 | -- | 35.72 | -- |
| def2-qzvppd | 1e-14 | 31284 | -- | 39.69 | -- |
| def2-qzvppd | 1e-12 | 31284 | -- | 35.42 | -- |
| cc-pvdz | 1e-14 | 6177 | 24.32 | 2.81 | 8.7 |
| cc-pvdz | 1e-12 | 6177 | 24.32 | 2.58 | 9.4 |
| cc-pvtz | 1e-14 | 14244 | 77.38 | 6.38 | 12.1 |
| cc-pvtz | 1e-12 | 14244 | 77.38 | 5.78 | 13.4 |
| cc-pvqz | 1e-14 | 27459 | 276.88 | 19.75 | 14.0 |
| cc-pvqz | 1e-12 | 27459 | 276.88 | 17.30 | 16.0 |
| aug-cc-pvdz | 1e-14 | 10380 | 46.86 | 8.18 | 5.7 |
| aug-cc-pvdz | 1e-12 | 10380 | 46.86 | 7.60 | 6.2 |
| aug-cc-pvtz | 1e-14 | 22311 | 157.03 | 27.15 | 5.8 |
| aug-cc-pvtz | 1e-12 | 22311 | 157.03 | 24.43 | 6.4 |
| aug-cc-pvqz | 1e-14 | 40674 | -- | 102.23 | -- |
| aug-cc-pvqz | 1e-12 | 40674 | -- | 88.48 | -- |

#### ubiquitin

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 11577 | 50.14 | 4.23 | 11.9 |
| def2-svp | 1e-12 | 11577 | 50.14 | 3.87 | 13.0 |
| def2-svpd | 1e-14 | 17433 | 109.38 | 11.52 | 9.5 |
| def2-svpd | 1e-12 | 17433 | 109.38 | 9.99 | 10.9 |
| def2-tzvp | 1e-14 | 22442 | 171.13 | 10.82 | 15.8 |
| def2-tzvp | 1e-12 | 22442 | 171.13 | 9.23 | 18.5 |
| def2-tzvpp | 1e-14 | 27479 | 230.49 | 13.22 | 17.4 |
| def2-tzvpp | 1e-12 | 27479 | 230.49 | 11.31 | 20.4 |
| def2-tzvpd | 1e-14 | 28298 | 255.80 | 29.35 | 8.7 |
| def2-tzvpd | 1e-12 | 28298 | 255.80 | 25.19 | 10.2 |
| def2-tzvppd | 1e-14 | 33335 | -- | 34.99 | -- |
| def2-tzvppd | 1e-12 | 33335 | -- | 30.13 | -- |
| def2-qzvp | 1e-14 | 53197 | -- | 50.19 | -- |
| def2-qzvp | 1e-12 | 53197 | -- | 42.43 | -- |
| def2-qzvpp | 1e-14 | 53197 | -- | 50.24 | -- |
| def2-qzvpp | 1e-12 | 53197 | -- | 42.38 | -- |
| def2-qzvpd | 1e-14 | 59053 | -- | 112.16 | -- |
| def2-qzvpd | 1e-12 | 59053 | -- | 93.82 | -- |
| def2-qzvppd | 1e-14 | 59053 | -- | 112.31 | -- |
| def2-qzvppd | 1e-12 | 59053 | -- | 94.45 | -- |
| cc-pvdz | 1e-14 | 11577 | 77.93 | 5.13 | 15.2 |
| cc-pvdz | 1e-12 | 11577 | 77.93 | 4.59 | 17.0 |
| cc-pvtz | 1e-14 | 26870 | 253.76 | 12.65 | 20.1 |
| cc-pvtz | 1e-12 | 26870 | 253.76 | 10.85 | 23.4 |
| cc-pvqz | 1e-14 | 51984 | -- | 43.52 | -- |
| cc-pvqz | 1e-12 | 51984 | -- | 36.36 | -- |
| aug-cc-pvdz | 1e-14 | 19511 | 160.83 | 19.37 | 8.3 |
| aug-cc-pvdz | 1e-12 | 19511 | 160.83 | 17.14 | 9.4 |
| aug-cc-pvtz | 1e-14 | 42163 | -- | 77.89 | -- |
| aug-cc-pvtz | 1e-12 | 42163 | -- | 66.24 | -- |
| aug-cc-pvqz | 1e-14 | 77098 | -- | 348.78 | -- |
| aug-cc-pvqz | 1e-12 | 77098 | -- | 294.58 | -- |

### The nuclear attraction driver against the reference

`OMP_NUM_THREADS=14`, best of five runs inside a two second budget, so the large cases
run once. Measured on 2026-09-15, after the driver was finished to angular momentum
six; the overlap and kinetic numbers above are an A/B across a set of changes and this
one is not, so read it as where the driver stands rather than as a movement.

**ref** is `CNuclearPotentialDriver`, which computes every atom pair and carries no
threshold, so one number serves both threshold rows. A dash marks a case it cannot run:
it returns a dense matrix, and the cut is at 30000 functions, six gigabytes. The bases
are those whose highest angular momentum is g, where the reference stops -- above that
it returns zeros, which is not a reference.

The geometric mean over the 170 comparable cases is **2.20**, and it is not one number:
it rises with the molecule, 1.45 on tagrisso to 4.93 on ubiquitin.

**This integral is the overlap times the number of nuclei.** Every kept pair is
evaluated against every charge, so tagrisso at def2-qzvp is 85 ms here against 1.24 ms
for the overlap, a factor of 68 on 70 atoms, and ubiquitin reaches 182 seconds. The
grid is otherwise the one the overlap and kinetic sections use.

The tagrisso, c60, taxol and paracetamol cluster tables were measured again after the
block size floor below; crambin and ubiquitin were not, because the floor cannot reach
them -- they choose 3685 and 13541 atom pairs a block, far above either value -- and a
spot check agrees, crambin at def2-tzvp giving 3459 ms against 3429 before.


#### tagrisso, 70 atoms

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 683 | 8.69 | 10.54 | 0.8 |
| def2-svp | 1e-12 | 683 | 8.69 | 10.09 | 0.9 |
| def2-svpd | 1e-14 | 1010 | 17.85 | 18.84 | 0.9 |
| def2-svpd | 1e-12 | 1010 | 17.85 | 18.35 | 1.0 |
| def2-tzvp | 1e-14 | 1345 | 46.00 | 29.21 | 1.6 |
| def2-tzvp | 1e-12 | 1345 | 46.00 | 28.43 | 1.6 |
| def2-tzvpp | 1e-14 | 1609 | 46.05 | 33.17 | 1.4 |
| def2-tzvpp | 1e-12 | 1609 | 46.05 | 31.64 | 1.5 |
| def2-tzvpd | 1e-14 | 1672 | 43.31 | 43.66 | 1.0 |
| def2-tzvpd | 1e-12 | 1672 | 43.31 | 42.38 | 1.0 |
| def2-tzvppd | 1e-14 | 1936 | 43.72 | 48.66 | 0.9 |
| def2-tzvppd | 1e-12 | 1936 | 43.72 | 46.64 | 0.9 |
| def2-qzvp | 1e-14 | 3099 | 186.46 | 87.89 | 2.1 |
| def2-qzvp | 1e-12 | 3099 | 186.46 | 84.68 | 2.2 |
| def2-qzvpp | 1e-14 | 3099 | 183.76 | 88.23 | 2.1 |
| def2-qzvpp | 1e-12 | 3099 | 183.76 | 84.44 | 2.2 |
| def2-qzvpd | 1e-14 | 3426 | 184.57 | 113.32 | 1.6 |
| def2-qzvpd | 1e-12 | 3426 | 184.57 | 109.68 | 1.7 |
| def2-qzvppd | 1e-14 | 3426 | 187.54 | 114.07 | 1.6 |
| def2-qzvppd | 1e-12 | 3426 | 187.54 | 110.23 | 1.7 |
| cc-pvdz | 1e-14 | 683 | 34.74 | 21.96 | 1.6 |
| cc-pvdz | 1e-12 | 683 | 34.74 | 21.44 | 1.6 |
| cc-pvtz | 1e-14 | 1572 | 51.72 | 41.48 | 1.2 |
| cc-pvtz | 1e-12 | 1572 | 51.72 | 40.19 | 1.3 |
| cc-pvqz | 1e-14 | 3025 | 191.37 | 91.43 | 2.1 |
| cc-pvqz | 1e-12 | 3025 | 191.37 | 87.45 | 2.2 |
| aug-cc-pvdz | 1e-14 | 1148 | 45.05 | 39.41 | 1.1 |
| aug-cc-pvdz | 1e-12 | 1148 | 45.05 | 38.31 | 1.2 |
| aug-cc-pvtz | 1e-14 | 2461 | 128.13 | 85.86 | 1.5 |
| aug-cc-pvtz | 1e-12 | 2461 | 128.13 | 83.46 | 1.5 |
| aug-cc-pvqz | 1e-14 | 4478 | 509.61 | 200.65 | 2.5 |
| aug-cc-pvqz | 1e-12 | 4478 | 509.61 | 193.55 | 2.6 |

#### c60, 60 atoms

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 840 | 16.33 | 15.03 | 1.1 |
| def2-svp | 1e-12 | 840 | 16.33 | 14.40 | 1.1 |
| def2-svpd | 1e-14 | 1200 | 22.92 | 25.37 | 0.9 |
| def2-svpd | 1e-12 | 1200 | 22.92 | 24.84 | 0.9 |
| def2-tzvp | 1e-14 | 1860 | 96.92 | 50.07 | 1.9 |
| def2-tzvp | 1e-12 | 1860 | 96.92 | 47.53 | 2.0 |
| def2-tzvpp | 1e-14 | 1860 | 94.27 | 49.44 | 1.9 |
| def2-tzvpp | 1e-12 | 1860 | 94.27 | 47.65 | 2.0 |
| def2-tzvpd | 1e-14 | 2220 | 96.39 | 67.91 | 1.4 |
| def2-tzvpd | 1e-12 | 2220 | 96.39 | 66.12 | 1.5 |
| def2-tzvppd | 1e-14 | 2220 | 91.19 | 68.32 | 1.3 |
| def2-tzvppd | 1e-12 | 2220 | 91.19 | 65.58 | 1.4 |
| def2-qzvp | 1e-14 | 3420 | 389.30 | 130.37 | 3.0 |
| def2-qzvp | 1e-12 | 3420 | 389.30 | 123.95 | 3.1 |
| def2-qzvpp | 1e-14 | 3420 | 391.52 | 129.68 | 3.0 |
| def2-qzvpp | 1e-12 | 3420 | 391.52 | 123.71 | 3.2 |
| def2-qzvpd | 1e-14 | 3780 | 381.38 | 159.09 | 2.4 |
| def2-qzvpd | 1e-12 | 3780 | 381.38 | 152.41 | 2.5 |
| def2-qzvppd | 1e-14 | 3780 | 368.34 | 158.53 | 2.3 |
| def2-qzvppd | 1e-12 | 3780 | 368.34 | 153.42 | 2.4 |
| cc-pvdz | 1e-14 | 840 | 35.61 | 33.12 | 1.1 |
| cc-pvdz | 1e-12 | 840 | 35.61 | 31.98 | 1.1 |
| cc-pvtz | 1e-14 | 1800 | 91.20 | 63.43 | 1.4 |
| cc-pvtz | 1e-12 | 1800 | 91.20 | 60.93 | 1.5 |
| cc-pvqz | 1e-14 | 3300 | 376.02 | 134.58 | 2.8 |
| cc-pvqz | 1e-12 | 3300 | 376.02 | 129.60 | 2.9 |
| aug-cc-pvdz | 1e-14 | 1380 | 48.98 | 56.46 | 0.9 |
| aug-cc-pvdz | 1e-12 | 1380 | 48.98 | 55.13 | 0.9 |
| aug-cc-pvtz | 1e-14 | 2760 | 120.06 | 121.65 | 1.0 |
| aug-cc-pvtz | 1e-12 | 2760 | 120.06 | 117.95 | 1.0 |
| aug-cc-pvqz | 1e-14 | 4800 | 453.16 | 276.83 | 1.6 |
| aug-cc-pvqz | 1e-12 | 4800 | 453.16 | 266.08 | 1.7 |

#### taxol, 110 atoms

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 1099 | 36.31 | 29.33 | 1.2 |
| def2-svp | 1e-12 | 1099 | 36.31 | 27.63 | 1.3 |
| def2-svpd | 1e-14 | 1657 | 46.65 | 58.24 | 0.8 |
| def2-svpd | 1e-12 | 1657 | 46.65 | 56.33 | 0.8 |
| def2-tzvp | 1e-14 | 2185 | 186.86 | 84.13 | 2.2 |
| def2-tzvp | 1e-12 | 2185 | 186.86 | 80.14 | 2.3 |
| def2-tzvpp | 1e-14 | 2577 | 172.59 | 94.41 | 1.8 |
| def2-tzvpp | 1e-12 | 2577 | 172.59 | 90.37 | 1.9 |
| def2-tzvpd | 1e-14 | 2743 | 182.34 | 134.05 | 1.4 |
| def2-tzvpd | 1e-12 | 2743 | 182.34 | 127.99 | 1.4 |
| def2-tzvppd | 1e-14 | 3135 | 186.80 | 148.95 | 1.3 |
| def2-tzvppd | 1e-12 | 3135 | 186.80 | 142.40 | 1.3 |
| def2-qzvp | 1e-14 | 4947 | 766.25 | 255.82 | 3.0 |
| def2-qzvp | 1e-12 | 4947 | 766.25 | 243.59 | 3.1 |
| def2-qzvpp | 1e-14 | 4947 | 776.72 | 256.66 | 3.0 |
| def2-qzvpp | 1e-12 | 4947 | 776.72 | 242.07 | 3.2 |
| def2-qzvpd | 1e-14 | 5505 | 768.16 | 349.44 | 2.2 |
| def2-qzvpd | 1e-12 | 5505 | 768.16 | 334.23 | 2.3 |
| def2-qzvppd | 1e-14 | 5505 | 760.93 | 348.85 | 2.2 |
| def2-qzvppd | 1e-12 | 5505 | 760.93 | 335.48 | 2.3 |
| cc-pvdz | 1e-14 | 1099 | 69.31 | 58.92 | 1.2 |
| cc-pvdz | 1e-12 | 1099 | 69.31 | 55.73 | 1.2 |
| cc-pvtz | 1e-14 | 2516 | 185.38 | 116.98 | 1.6 |
| cc-pvtz | 1e-12 | 2516 | 185.38 | 110.07 | 1.7 |
| cc-pvqz | 1e-14 | 4825 | 779.12 | 260.54 | 3.0 |
| cc-pvqz | 1e-12 | 4825 | 779.12 | 247.23 | 3.2 |
| aug-cc-pvdz | 1e-14 | 1844 | 98.85 | 116.55 | 0.8 |
| aug-cc-pvdz | 1e-12 | 1844 | 98.85 | 112.56 | 0.9 |
| aug-cc-pvtz | 1e-14 | 3933 | 309.31 | 261.33 | 1.2 |
| aug-cc-pvtz | 1e-12 | 3933 | 309.31 | 250.47 | 1.2 |
| aug-cc-pvqz | 1e-14 | 7134 | 1379.48 | 617.56 | 2.2 |
| aug-cc-pvqz | 1e-12 | 7134 | 1379.48 | 591.25 | 2.3 |

#### paracetamol_cluster, 320 atoms

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 3184 | 512.15 | 316.39 | 1.6 |
| def2-svp | 1e-12 | 3184 | 512.15 | 287.20 | 1.8 |
| def2-svpd | 1e-14 | 4768 | 923.26 | 801.14 | 1.2 |
| def2-svpd | 1e-12 | 4768 | 923.26 | 728.53 | 1.3 |
| def2-tzvp | 1e-14 | 6320 | 2149.43 | 942.27 | 2.3 |
| def2-tzvp | 1e-12 | 6320 | 2149.43 | 857.68 | 2.5 |
| def2-tzvpp | 1e-14 | 7472 | 2695.89 | 1081.83 | 2.5 |
| def2-tzvpp | 1e-12 | 7472 | 2695.89 | 983.99 | 2.7 |
| def2-tzvpd | 1e-14 | 7904 | 2941.29 | 1834.47 | 1.6 |
| def2-tzvpd | 1e-12 | 7904 | 2941.29 | 1684.57 | 1.7 |
| def2-tzvppd | 1e-14 | 9056 | 3586.70 | 2056.64 | 1.7 |
| def2-tzvppd | 1e-12 | 9056 | 3586.70 | 1880.71 | 1.9 |
| def2-qzvp | 1e-14 | 14352 | 13601.73 | 2930.58 | 4.6 |
| def2-qzvp | 1e-12 | 14352 | 13601.73 | 2686.37 | 5.1 |
| def2-qzvpp | 1e-14 | 14352 | 13714.12 | 2939.45 | 4.7 |
| def2-qzvpp | 1e-12 | 14352 | 13714.12 | 2679.86 | 5.1 |
| def2-qzvpd | 1e-14 | 15936 | 15498.40 | 4696.47 | 3.3 |
| def2-qzvpd | 1e-12 | 15936 | 15498.40 | 4383.05 | 3.5 |
| def2-qzvppd | 1e-14 | 15936 | 15569.90 | 4690.95 | 3.3 |
| def2-qzvppd | 1e-12 | 15936 | 15569.90 | 4324.47 | 3.6 |
| cc-pvdz | 1e-14 | 3184 | 940.98 | 579.41 | 1.6 |
| cc-pvdz | 1e-12 | 3184 | 940.98 | 538.86 | 1.7 |
| cc-pvtz | 1e-14 | 7296 | 3050.26 | 1226.12 | 2.5 |
| cc-pvtz | 1e-12 | 7296 | 3050.26 | 1121.17 | 2.7 |
| cc-pvqz | 1e-14 | 14000 | 13546.61 | 2815.16 | 4.8 |
| cc-pvqz | 1e-12 | 14000 | 13546.61 | 2570.02 | 5.3 |
| aug-cc-pvdz | 1e-14 | 5344 | 1641.70 | 1618.09 | 1.0 |
| aug-cc-pvdz | 1e-12 | 5344 | 1641.70 | 1510.26 | 1.1 |
| aug-cc-pvtz | 1e-14 | 11408 | 6620.42 | 3724.20 | 1.8 |
| aug-cc-pvtz | 1e-12 | 11408 | 6620.42 | 3465.30 | 1.9 |
| aug-cc-pvqz | 1e-14 | 20704 | 31557.92 | 8738.65 | 3.6 |
| aug-cc-pvqz | 1e-12 | 20704 | 31557.92 | 8093.69 | 3.9 |

#### crambin, 642 atoms

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 6177 | 3649.58 | 1304.02 | 2.8 |
| def2-svp | 1e-12 | 6177 | 3649.58 | 1156.76 | 3.2 |
| def2-svpd | 1e-14 | 9294 | 6706.62 | 3763.13 | 1.8 |
| def2-svpd | 1e-12 | 9294 | 6706.62 | 3349.46 | 2.0 |
| def2-tzvp | 1e-14 | 12063 | 15350.67 | 3924.63 | 3.9 |
| def2-tzvp | 1e-12 | 12063 | 15350.67 | 3502.92 | 4.4 |
| def2-tzvpp | 1e-14 | 14613 | 20186.36 | 4609.76 | 4.4 |
| def2-tzvpp | 1e-12 | 14613 | 20186.36 | 4085.00 | 4.9 |
| def2-tzvpd | 1e-14 | 15180 | 21420.68 | 8687.95 | 2.5 |
| def2-tzvpd | 1e-12 | 15180 | 21420.68 | 7721.24 | 2.8 |
| def2-tzvppd | 1e-14 | 17730 | 26492.34 | 9748.20 | 2.7 |
| def2-tzvppd | 1e-12 | 17730 | 26492.34 | 8928.92 | 3.0 |
| def2-qzvp | 1e-14 | 28167 | 102744.48 | 12745.29 | 8.1 |
| def2-qzvp | 1e-12 | 28167 | 102744.48 | 11388.54 | 9.0 |
| def2-qzvpp | 1e-14 | 28167 | 102655.91 | 12742.22 | 8.1 |
| def2-qzvpp | 1e-12 | 28167 | 102655.91 | 11367.09 | 9.0 |
| def2-qzvpd | 1e-14 | 31284 | -- | 22167.71 | -- |
| def2-qzvpd | 1e-12 | 31284 | -- | 19981.53 | -- |
| def2-qzvppd | 1e-14 | 31284 | -- | 22159.78 | -- |
| def2-qzvppd | 1e-12 | 31284 | -- | 19956.65 | -- |
| cc-pvdz | 1e-14 | 6177 | 6853.00 | 2308.52 | 3.0 |
| cc-pvdz | 1e-12 | 6177 | 6853.00 | 2102.67 | 3.3 |
| cc-pvtz | 1e-14 | 14244 | 22937.01 | 5148.58 | 4.5 |
| cc-pvtz | 1e-12 | 14244 | 22937.01 | 4624.55 | 5.0 |
| cc-pvqz | 1e-14 | 27459 | 103033.38 | 12029.83 | 8.6 |
| cc-pvqz | 1e-12 | 27459 | 103033.38 | 10779.14 | 9.6 |
| aug-cc-pvdz | 1e-14 | 10380 | 12116.90 | 8156.40 | 1.5 |
| aug-cc-pvdz | 1e-12 | 10380 | 12116.90 | 7433.24 | 1.6 |
| aug-cc-pvtz | 1e-14 | 22311 | 49477.94 | 18930.62 | 2.6 |
| aug-cc-pvtz | 1e-12 | 22311 | 49477.94 | 17228.01 | 2.9 |
| aug-cc-pvqz | 1e-14 | 40674 | -- | 44352.28 | -- |
| aug-cc-pvqz | 1e-12 | 40674 | -- | 40141.93 | -- |

#### ubiquitin, 1231 atoms

| basis | threshold | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- |
| def2-svp | 1e-14 | 11577 | 23723.41 | 5004.68 | 4.7 |
| def2-svp | 1e-12 | 11577 | 23723.41 | 4401.56 | 5.4 |
| def2-svpd | 1e-14 | 17433 | 43953.05 | 15609.01 | 2.8 |
| def2-svpd | 1e-12 | 17433 | 43953.05 | 13711.35 | 3.2 |
| def2-tzvp | 1e-14 | 22442 | 98955.07 | 15024.17 | 6.6 |
| def2-tzvp | 1e-12 | 22442 | 98955.07 | 13180.53 | 7.5 |
| def2-tzvpp | 1e-14 | 27479 | 130750.69 | 18010.52 | 7.3 |
| def2-tzvpp | 1e-12 | 27479 | 130750.69 | 15690.41 | 8.3 |
| def2-tzvpd | 1e-14 | 28298 | 137340.08 | 35865.00 | 3.8 |
| def2-tzvpd | 1e-12 | 28298 | 137340.08 | 31627.84 | 4.3 |
| def2-tzvppd | 1e-14 | 33335 | -- | 41359.97 | -- |
| def2-tzvppd | 1e-12 | 33335 | -- | 35865.71 | -- |
| def2-qzvp | 1e-14 | 53197 | -- | 51076.91 | -- |
| def2-qzvp | 1e-12 | 53197 | -- | 44789.30 | -- |
| def2-qzvpp | 1e-14 | 53197 | -- | 50766.14 | -- |
| def2-qzvpp | 1e-12 | 53197 | -- | 44595.72 | -- |
| def2-qzvpd | 1e-14 | 59053 | -- | 95265.42 | -- |
| def2-qzvpd | 1e-12 | 59053 | -- | 84079.52 | -- |
| def2-qzvppd | 1e-14 | 59053 | -- | 95887.50 | -- |
| def2-qzvppd | 1e-12 | 59053 | -- | 83948.12 | -- |
| cc-pvdz | 1e-14 | 11577 | 43887.49 | 8254.35 | 5.3 |
| cc-pvdz | 1e-12 | 11577 | 43887.49 | 7321.89 | 6.0 |
| cc-pvtz | 1e-14 | 26870 | 148457.05 | 19260.87 | 7.7 |
| cc-pvtz | 1e-12 | 26870 | 148457.05 | 16926.39 | 8.8 |
| cc-pvqz | 1e-14 | 51984 | -- | 46050.47 | -- |
| cc-pvqz | 1e-12 | 51984 | -- | 40661.67 | -- |
| aug-cc-pvdz | 1e-14 | 19511 | 78025.18 | 37041.03 | 2.1 |
| aug-cc-pvdz | 1e-12 | 19511 | 78025.18 | 32973.26 | 2.4 |
| aug-cc-pvtz | 1e-14 | 42163 | -- | 88288.37 | -- |
| aug-cc-pvtz | 1e-12 | 42163 | -- | 77356.21 | -- |
| aug-cc-pvqz | 1e-14 | 77098 | -- | 203911.70 | -- |
| aug-cc-pvqz | 1e-12 | 77098 | -- | 181524.88 | -- |

#### The block size, and a floor of its own

`make_block_size` gives `npairs / (4 x nthreads)` and floors it at
`sparsity::min_block_size`, which is 256. Swept at fourteen threads, every case has the
same shape -- steep below 512 atom pairs a block, flat from there to 32768:

| case | default | 64 | 256 | 512 | 1024 | 4096 | 8192 | 32768 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| c60 def2-svpd | 25.74 | 41.50 | 25.64 | **24.58** | 24.61 | 25.02 | 24.99 | 25.56 |
| taxol def2-svpd | 58.41 | 85.24 | 57.76 | **55.17** | 55.22 | 57.14 | 56.42 | 56.84 |
| tagrisso def2-tzvpd | 44.08 | 56.76 | 43.62 | 41.89 | 42.37 | **41.83** | 42.45 | 43.53 |
| crambin def2-tzvp | 3429 | 9152 | 5342 | 4223 | 3651 | **3427** | 3475 | 3525 |
| ubiquitin def2-svp | 4360 | 14022 | 9812 | 7142 | 5520 | 4400 | **4310** | 4449 |

**It is a fixed cost per block and not the buffer.** The arena is the largest block's
pairs times the buffer rows of its momenta, and over this sweep it runs from 0.1 to 56
MB a thread with no structure in the timings at all -- crambin's best point is the 23 MB
one. It is not a scheduling effect either: c60 at four threads has the same shape, 53.60
at 256 against 52.29 at 512, and at one thread the same ordering with the gain down to
0.6 per cent.

**So the floor binds, and only for the small molecules.** crambin chooses 3685 and
ubiquitin 13541, both on the plateau and within one per cent of the best point in the
sweep -- there is nothing to win there. Everything under about 150 atoms falls to the
floor instead, and 256 costs it 4 to 6 per cent.

The nuclear potential driver therefore carries a floor of its own, 512, passed as
`min_pairs` through `make_pattern`. `sparsity::min_block_size` stays at 256: it was
*lowered* to that from 2048 on the evidence of the overlap driver, whose buffer is ten
rows against this one's hundreds, and reaching into it would change drivers this was not
measured on. What it bought, as a ratio of the committed numbers over the same grid:

| molecule | chosen block size | speedup |
| --- | ---: | ---: |
| tagrisso | 256 -> 512 | 1.057 |
| c60 | 256 -> 512 | 1.063 |
| taxol | 256 -> 512 | 1.058 |
| paracetamol cluster | 917, unchanged | 0.996 |

The paracetamol cluster is the control: its chosen size already exceeds both floors, and
it does not move. Three of the sixteen rows which lost to the reference now win.

**This is a laptop number and the floor is the kind of constant which inverts.** At 512
c60 holds four blocks where 256 gave it seven, and the small molecules are the ones with
the fewest tasks to spread over a node. The plateau is broad enough that 512 is not a
risky point, but it has not been measured above fourteen threads.

#### What these numbers say

**The advantage is screening, so it grows with the molecule and nothing else.** The
geometric mean by molecule: tagrisso 1.45, c60 1.64, taxol 1.70, paracetamol cluster
2.45, crambin 3.80, ubiquitin 4.93. The reference computes every atom pair; the driver
computes the pairs which survive, and on seventy atoms almost all of them do. The best
case in the grid is crambin at cc-pvqz, 9.6. At the other end **thirteen rows, seven
basis and molecule combinations, are slower than the reference**: tagrisso in def2-svp,
def2-svpd and def2-tzvppd, c60 in def2-svpd and aug-cc-pvdz, and taxol in def2-svpd and
aug-cc-pvdz, from 0.8 to 0.9. Every one is on the three smallest molecules and every one
is double or triple zeta. Little to screen, and the fixed cost of a block is what is
left.

**Diffuse functions are what defeats it, and they cost more the larger the molecule.**
At 1e-12, going from cc-pvdz to aug-cc-pvdz:

| molecule | nao | time | x ref |
| --- | ---: | ---: | --- |
| tagrisso | x1.68 | x1.79 | 1.6 -> 1.2 |
| crambin | x1.68 | x3.54 | 3.3 -> 1.6 |
| ubiquitin | x1.69 | x4.50 | 6.0 -> 2.4 |

The function count grows by the same 1.7 in all three; the time grows by 1.8 on seventy
atoms and by 4.5 on twelve hundred. The reference is indifferent -- it had no screening
to lose. So the ratio falling is the driver giving back exactly what it had gained, and
the larger the molecule the more there was to give back. def2-svp to def2-svpd is the
same story at 1.5 functions and 3.1 time, def2-tzvp to def2-tzvpd at 1.26 and 2.40.

**The threshold is worth between three and fourteen per cent, in that order.** 1e-12
against 1e-14: tagrisso 1.035, c60 1.038, taxol 1.048, paracetamol cluster 1.088,
crambin 1.113, ubiquitin 1.137. It buys nothing where there is nothing to screen and
the most where the screening is already doing the work -- the same axis as everything
else here. Since the reference has no threshold at all, the whole of it shows up in the
ratio, which is why the 1e-12 mean is 2.27 against 2.13.

### Where the time of the nuclear attraction driver goes

Profiled on 2026-09-15, fourteen threads, best of five. Nothing here is a share read
off a sampler: `sample` over-reports the allocator badly enough that two rounds of work
were once planned off its percentages and both measured as noise. Every number below is
a wall clock A/B, an experiment which removes a cost and re-times, taken against the
same baseline in one session.

#### The charge loop is the whole call

This operator has a knob the others do not: the number of charges. Timing the same
molecule and basis against k of them separates what is paid per charge from what is
paid once, with no instrumentation at all.

| case | fit | fixed, at k = N |
| --- | --- | ---: |
| tagrisso def2-svp, 70 atoms | 0.54 ms + 138 us x k | 5.3% |
| taxol def2-tzvp, 110 atoms | 1.07 ms + 715 us x k | 1.3% |
| crambin def2-tzvp, 642 atoms | 3.58 ms + 5419 us x k | 0.1% |

T(k) is linear over three decades of k. So the sparsity pattern, the coordinates, the
task list, the HRR transfer, the harmonic transform, the distributor and the dense fill
are together between five per cent and one part in a thousand. **There is nothing to
optimise outside the charge loop**, which is why no driver-level phase timers appear
below: they would be measuring noise.

#### Inside it

The body evaluated per charge is `compute_pc`, `compute_full_npot_boys_function`, the
VRR ladder, and `contract_primitives`. Each of the first three was ablated in turn --
the Boys ladder replaced by a fill of the same rows, the pair exponential by the factor
it multiplies, the accumulation by one column instead of ncols -- so that the writes and
the dependencies stay and nothing is eliminated as dead:

| case | baseline | Boys ladder | pair exp | contraction | remainder |
| --- | ---: | ---: | ---: | ---: | ---: |
| tagrisso def2-svp | 10.25 ms | 47.6% | 17.1% | 1.2% | 34.1% |
| taxol def2-tzvp | 79.53 ms | 39.3% | 21.9% | 3.6% | 35.2% |
| crambin def2-tzvp | 3470 ms | 34.5% | 26.2% | 5.3% | 34.0% |

Removing the Boys ladder and the pair exponential together gives 67.5, 60.9 and 61.1
per cent against 64.7, 61.2 and 60.7 for the two measured separately -- **additive
within three points**, which is the check that the attributions mean anything. The
remainder is by subtraction: the VRR ladder, `compute_pc`, `_make_scaled_arguments`,
the transform, and the once-paid part above.

The Boys share falls as the angular momentum rises, 48 to 35 per cent, because the VRR
ladder grows faster than the order of the Boys function does. The contraction is
vectorised already and is small.

#### The pair exponential is recomputed thousands of times

`_scale_pair_values` carries `std::exp(-mu * ab_2[k])` in its inner loop:

```
    for (size_t j = 0; j < nrows; j++)
    {
        auto *row = buffer.data(target + 1 + j);

        for (size_t k = 0; k < ncols; k++)
        {
            row[k] *= fj * std::exp(-mu * ab_2[k]);
        }
    }
```

`mu` is the pair of primitives and `ab_2` the atom pair. **Neither depends on the row or
on the charge**, and the call sits inside the loop over charges. For crambin at
def2-tzvp that is of order seven rows times six hundred and forty two charges -- some
four thousand evaluations of an exponential where one would do -- and it is 26 per cent
of the call. The share rises with the molecule for exactly that reason: the charge loop
is the molecule.

This is the same shape as the three-center kernel's `e_ab`, which recomputes its
exponential for every atom on the c side and was 17 per cent there.

#### Hoisting it out of the row loop

Formed once for the call into a scratch held per thread, every row then scaled by it.
The A/B is the same binary either way, both measured in one session:

| case | before | after | x |
| --- | ---: | ---: | ---: |
| tagrisso def2-svp | 10.24 ms | 9.26 ms | 1.106 |
| taxol def2-tzvp | 79.86 ms | 70.38 ms | 1.135 |
| crambin def2-tzvp | 3465.40 ms | 2966.55 ms | **1.168** |

The gain rises with the molecule because the charge loop is the molecule. What is left
of the 26 per cent is one `ncols` of exponentials per charge, where the operator wants
one per pair of primitives.

#### And then out of the charge loop

Removing that last part means keeping the values across the charge loop, which needs a
row of the buffer to keep them in. The generator now gives one: `BufferLayout` carries a
`pair_exp` section of a single row for an anchored operator, `compute_pair_exponent`
fills it once for the pair of primitives above the loop over the charges, and the Boys
wrapper takes that row where it used to take the reduced exponent. Every buffer row
behind it moves by one and the table of buffer rows counts one more.

| case | before | row loop | charge loop | whole |
| --- | ---: | ---: | ---: | ---: |
| tagrisso def2-svp | 10.24 ms | 9.26 ms | 8.07 ms | **1.269** |
| taxol def2-tzvp | 79.86 ms | 70.38 ms | 61.61 ms | **1.296** |
| crambin def2-tzvp | 3465.40 ms | 2966.55 ms | 2536.17 ms | **1.366** |

**The ablation predicted the ceiling and the change reached it.** Removing the
exponential entirely measured 26.2 per cent of crambin at def2-tzvp; forming it once a
pair of primitives instead of once a row of every charge took 26.8 per cent off. The
exponential is now evaluated `ncols` times for a pair of primitives where it was
evaluated `nrows` times `ncharges` times `ncols`, which on that case is four thousand
evaluations down to one.

#### The three-center kernels, once they were measured properly

The earlier reading of these was worthless: two cases, one run each, one of them a
fourteen second call which wanders by six per cent between runs of code nothing changed.
Measured again on four cases small enough to repeat, best of five with spreads under
1.12, ablating the exponential the same way:

| case | baseline | no exponential | its share |
| --- | ---: | ---: | ---: |
| tagrisso def2-svp | 286.62 ms | 247.52 ms | 13.6% |
| c60 def2-svp | 784.92 ms | 678.78 ms | 13.5% |
| taxol def2-svp | 999.54 ms | 861.77 ms | 13.8% |
| tagrisso def2-tzvp | 1094.17 ms | 980.87 ms | 10.4% |

The same change was tried first and reached only two thirds of that, because the loop
nest is not the nuclear attraction's. There the charges stand inside the pair of
primitives, so one row above them kills the whole repetition. Here the nest is atoms on
the c side, then the bra's two primitives, then the ket's: `mu` is the bra pair's and is
fixed only *inside* the loop over the atoms, so a row is refilled for every one of them.
That version measured 1.078 to 1.100.

**Reordering the nest is not the way to fix that.** With the atoms innermost every one
of them needs its partial sum live at once, and the contracted rows are reused per atom
today -- zeroed at the top, accumulated, transformed into that atom's slice. The section
would become `contracted x natoms`: 6 rows to 3.9 thousand for `(ss|d)`, 609 to 391
thousand for `(dd|g)`, 15176 to 9.7 million for `(ii|i)`. At a few thousand columns the
middle one alone is hundreds of gigabytes.

**The screening does not stand in the way of it, though**, which was worth checking
before ruling it out: the primitive bound neglects the position of the atom on the ket
side, so `dimensions` is indexed by the three primitives alone and a reorder would leave
it untouched. The kernels say so in a NOTE.

What works instead is to keep the nest and hold *every* pair's exponential at once, in a
scratch of `nprim_a * nprim_b` runs of atom pairs, filled once for the call above the
loop over the atoms and indexed by the pair. It costs no buffer row, no reorder and no
change to the accumulation:

| case | before | after | x | of the ceiling |
| --- | ---: | ---: | ---: | ---: |
| tagrisso def2-svp | 286.62 ms | 247.03 ms | 1.160 | 101% |
| c60 def2-svp | 784.92 ms | 676.12 ms | 1.161 | 102% |
| taxol def2-svp | 999.54 ms | 858.34 ms | 1.165 | 102% |
| tagrisso def2-tzvp | 1094.17 ms | 983.24 ms | 1.113 | 98% |

**It reaches the ablation**, which is the whole of what was there to take: removing the
exponential outright gave 247.52, 678.78, 861.77 and 980.87 ms, and forming every pair's
once gives 247.03, 676.12, 858.34 and 983.24.

The price is the scratch, and it is small: at most `nprim_a * nprim_b` runs of atom
pairs a thread, which is 25 for these def2-svp cases and 36 for def2-tzvp, against
buffers that are 22 rows for `(ss|d)` and 4998 for `(dd|g)`. Peak resident size is
unmoved -- 13.8 GB for taxol at def2-svp, which is the tensor it returns.

**The three-center kernels deliberately keep the old form.** They call the same
function and the same argument applies to them, but measured once each way the two
cases disagreed -- tagrisso/def2-svp 290 to 264 ms, taxol/def2-tzvp 13.49 to 15.59 s --
and at thirteen seconds a run both numbers are single samples. The change trades an
exponential for a stream of `ncols` doubles, which is a different trade for a kernel
whose blocks are larger, and 441 kernels is too many to move on a coin toss. With the
three-center path left alone it measures 287 ms and 13.79 s against 290 ms and 13.49 s,
which is where it was.

### The two-center Coulomb driver against the reference

The bases are the fitting sets, which is what this operator is used with. A dash in the
ref column marks a case the reference dispatcher cannot reach, which is everything
above angular momentum six. This grid merges the def2 and the correlation consistent
fitting sets, which the sections above keep apart, so its mean is not comparable with
the 1.07 and 0.96 recorded for them separately. The largest crambin and ubiquitin cases
above g were left out for runtime.

The geometric mean over the 44 cases with a reference is **1.21**. This driver stays the
weak one, and for a structural reason: nothing screens here, so both sides compute
every atom pair and the comparison is kernel against kernel rather than work against
work.

| molecule | basis | lmax | nao | ref ms | simd ms | x ref |
| --- | --- | --- | --- | --- | --- | --- |
| tagrisso | jfit | 4 | 2176 | 4.93 | 2.25 | 2.2 |
| tagrisso | jkfit | 4 | 3387 | 5.84 | 3.53 | 1.7 |
| tagrisso | cc-pvdz-rifit | 3 | 2534 | 2.07 | 1.72 | 1.2 |
| tagrisso | cc-pvtz-rifit | 4 | 3987 | 9.25 | 3.86 | 2.4 |
| tagrisso | cc-pvqz-rifit | 5 | 6699 | 15.41 | 13.88 | 1.1 |
| tagrisso | cc-pv5z-rifit | 6 | 10144 | 44.26 | 49.03 | 0.9 |
| tagrisso | aug-cc-pvdz-rifit | 3 | 3423 | 2.91 | 2.71 | 1.1 |
| tagrisso | aug-cc-pvtz-rifit | 4 | 5440 | 11.49 | 7.20 | 1.6 |
| tagrisso | aug-cc-pvqz-rifit | 5 | 8856 | 28.75 | 26.51 | 1.1 |
| c60 | jfit | 4 | 2940 | 7.59 | 3.32 | 2.3 |
| c60 | jkfit | 4 | 4500 | 11.17 | 6.06 | 1.8 |
| c60 | cc-pvdz-rifit | 3 | 3360 | 5.94 | 2.45 | 2.4 |
| c60 | cc-pvtz-rifit | 4 | 4860 | 6.27 | 5.50 | 1.1 |
| c60 | cc-pvqz-rifit | 5 | 7920 | 20.50 | 20.87 | 1.0 |
| c60 | cc-pv5z-rifit | 6 | 11580 | 58.95 | 79.14 | 0.7 |
| c60 | aug-cc-pvdz-rifit | 3 | 4320 | 4.42 | 3.86 | 1.1 |
| c60 | aug-cc-pvtz-rifit | 4 | 6360 | 10.86 | 10.22 | 1.1 |
| c60 | aug-cc-pvqz-rifit | 5 | 10080 | 35.33 | 39.67 | 0.9 |
| taxol | jfit | 4 | 3528 | 5.41 | 4.45 | 1.2 |
| taxol | jkfit | 4 | 5489 | 11.71 | 8.16 | 1.4 |
| taxol | cc-pvdz-rifit | 3 | 4102 | 3.87 | 3.54 | 1.1 |
| taxol | cc-pvtz-rifit | 4 | 6411 | 10.44 | 9.30 | 1.1 |
| taxol | cc-pvqz-rifit | 5 | 10747 | 35.88 | 34.87 | 1.0 |
| taxol | cc-pv5z-rifit | 6 | 16232 | 106.50 | 127.84 | 0.8 |
| taxol | aug-cc-pvdz-rifit | 3 | 5519 | 7.08 | 6.13 | 1.2 |
| taxol | aug-cc-pvtz-rifit | 4 | 8720 | 20.70 | 17.95 | 1.2 |
| taxol | aug-cc-pvqz-rifit | 5 | 14168 | 66.53 | 67.64 | 1.0 |
| Cu_PPh3_4_cation | jfit | 4 | 4483 | 9.24 | 6.93 | 1.3 |
| Cu_PPh3_4_cation | jkfit | 6 | 7256 | 19.09 | 14.74 | 1.3 |
| Cu_PPh3_4_cation | cc-pvtz-rifit | 6 | 8364 | 23.82 | 16.16 | 1.5 |
| Cu_PPh3_4_cation | cc-pvqz-rifit | 7 | 13798 | -- | 59.47 | -- |
| Cu_PPh3_4_cation | cc-pv5z-rifit | 8 | 20756 | -- | 218.97 | -- |
| Cu_PPh3_4_cation | aug-cc-pvtz-rifit | 6 | 11273 | 38.34 | 31.12 | 1.2 |
| Cu_PPh3_4_cation | aug-cc-pvqz-rifit | 7 | 18098 | -- | 116.00 | -- |
| paracetamol_cluster | jfit | 4 | 10208 | 38.87 | 37.66 | 1.0 |
| paracetamol_cluster | jkfit | 4 | 15888 | 82.02 | 71.45 | 1.1 |
| paracetamol_cluster | cc-pvdz-rifit | 3 | 11872 | 37.81 | 29.85 | 1.3 |
| paracetamol_cluster | cc-pvtz-rifit | 4 | 18576 | 88.97 | 83.13 | 1.1 |
| paracetamol_cluster | cc-pvqz-rifit | 5 | 31152 | 300.30 | 321.38 | 0.9 |
| paracetamol_cluster | cc-pv5z-rifit | 6 | 47072 | -- | 1355.72 | -- |
| paracetamol_cluster | aug-cc-pvdz-rifit | 3 | 15984 | 59.33 | 54.14 | 1.1 |
| paracetamol_cluster | aug-cc-pvtz-rifit | 4 | 25280 | 171.71 | 162.20 | 1.1 |
| paracetamol_cluster | aug-cc-pvqz-rifit | 5 | 41088 | -- | 730.90 | -- |
| crambin | jfit | 4 | 19500 | 136.51 | 135.97 | 1.0 |
| crambin | jkfit | 4 | 30751 | 274.47 | 263.42 | 1.0 |
| crambin | cc-pvdz-rifit | 3 | 22842 | 111.50 | 107.04 | 1.0 |
| crambin | cc-pvtz-rifit | 4 | 36183 | 319.40 | 381.11 | 0.8 |
| crambin | aug-cc-pvdz-rifit | 3 | 30909 | 219.37 | 196.04 | 1.1 |
| crambin | aug-cc-pvtz-rifit | 4 | 49398 | -- | 774.96 | -- |
| ubiquitin | jfit | 4 | 36419 | 1027.64 | 541.07 | 1.9 |
| ubiquitin | jkfit | 4 | 56971 | -- | 1124.64 | -- |
| ubiquitin | cc-pvdz-rifit | 3 | 42538 | -- | 471.45 | -- |
| ubiquitin | cc-pvtz-rifit | 4 | 67673 | -- | 1384.70 | -- |
| ubiquitin | aug-cc-pvdz-rifit | 3 | 57831 | -- | 891.39 | -- |

### What these numbers say

The overlap and the kinetic energy drivers moved from 3.11 and 2.51 to 3.89 and 3.94
against the reference on the same grids, and the gain is largest where it was needed:
tagrisso 1.65 to 2.41, c60 1.50 to 2.04, taxol 2.13 to 2.67. Those are the molecules
whose blocks are few and small, which is exactly what the fixed cost of a block was
punishing. The largest cases move four to six per cent, as their blocks exceed the
bound and keep the radix.

What has not been fixed is the floor. Fitting `T = S + P/N` from one and fourteen
threads, the part which does not scale is still one to four milliseconds, and it is 44
to 98 per cent of the wall time of a call: for everything under six hundred atoms in a
double or triple zeta basis it is 85 per cent or more. That number does not shrink with
cores, so on a large node it is what remains. sort_by_distance was 37 to 46 per cent of
it and is now much smaller; what is left is a tail of small phases -- the coordinates,
the pair groups, the task list, the diagonal blocks, the description -- of three
hundredths to a tenth of a millisecond each, several of which get *worse* as threads are
added.

## The three-center Coulomb driver

The driver takes a molecular basis and an auxiliary basis and returns a
`CSparseTensor`. Its kernels cover every combination to angular momentum six on the
two bra sides and eight on the auxiliary side, four hundred and forty one of them,
dispatched on the three momenta taken as one index.

### What it is checked against

`ThreeCenterElectronRepulsionDriver` dispatches to momentum four on the bra and six on
the auxiliary side, and **returns zeros without an error above that**. It can
therefore validate 175 of the 441 kernels, and the rest -- everything reaching h or i
on the bra, or k or l on the auxiliary side, `(ii|l)` among them -- has no independent
check. `tests/test_simd_three_center_electron_repulsion.py` sweeps the 175 element by
element and compares 1.5 million integrals; the worst disagreement is 3.1e-13, and the
error grows smoothly with angular momentum as accumulated round off should.

Two properties of that test are worth keeping. It proves its own index mapping on
`(ss|s)` and `(ps|s)` before the sweep, and it asserts on the number of elements it
compared, so a mapping which visits nothing fails rather than passes. **Summed values
must never be compared**: `T3FlatBuffer` keeps the upper triangle of an atom pair
while `CSparseTensor` keeps its own blocks, so their totals differ by tens of per cent
even when every integral agrees to 1e-15.

### The size of a block

A thread holds one buffer, spanning the largest combination any block carries, and its
size is the rows that combination needs times the atom pairs of a block. The rows
range over two orders of magnitude: three thousand for cc-pVDZ against its fitting
set, above half a million for `(ii|l)`. A fixed number of atom pairs therefore either
starves an ordinary combination or lets an extreme one run away, so the atom pairs are
chosen from a budget of 256 MB for the buffer, bounded to at most 256 and at least 8.

Measured by sweeping the size at fourteen threads, best of three:

| atom pairs | caffeine, def2-TZVP + jfit | tagrisso, def2-SVP + jkfit | c60, cc-pVDZ + RIFIT |
| --- | --- | --- | --- |
| 8 | 138.3 ms | 724.2 ms | 1803.9 ms |
| 32, the fixed size before | 99.4 | 417.3 | 1102.8 |
| 128 | 93.4 | 332.8 | 863.0 |
| 256 | 92.2 | 310.8 | 837.6 |
| 512 | 92.2 | 310.7 | 840.0 |

The curve is flat past 256, which is where the ceiling sits. Ordinary combinations
reach that ceiling far below the budget -- c60 holds 5.9 MB a thread, caffeine in
def2-TZVP 29.5 -- while `(ii|l)` is held to 57 atom pairs and the full 256 MB.

### Against the reference

`OMP_NUM_THREADS=14`. Every case runs in its own process, so the peak resident size is
its own; the driver is warmed up once and then timed best of two to five runs within a
four second budget. **ref** is `ThreeCenterElectronRepulsionDriver`, which carries no
threshold, so one number serves both threshold columns. The grid stops where the dense
tensor of the reference stops fitting, at twelve gigabytes.

| molecule | orbital | aux | nao | naux | ref ms | 1e-12 ms | 1e-14 ms | x ref | ref GB | simd GB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Caffeine | aug-cc-pvdz | rifit | 412 | 1238 | 91.2 | 81.3 | 82.0 | 1.11 | 0.78 | 0.81 |
| Caffeine | aug-cc-pvtz | rifit | 874 | 1944 | 781.8 | 490.3 | 492.3 | 1.59 | 5.54 | 5.56 |
| Caffeine | cc-pvdz | rifit | 246 | 924 | 44.5 | 32.0 | 32.8 | 1.36 | 0.21 | 0.21 |
| Caffeine | cc-pvqz | rifit | 1070 | 2398 | 1622.5 | 1099.7 | 1119.3 | 1.45 | 10.24 | 9.03 |
| Caffeine | cc-pvtz | rifit | 560 | 1434 | 194.4 | 155.7 | 158.8 | 1.22 | 1.68 | 1.59 |
| Caffeine | def2-qzvp | jfit | 1098 | 796 | 571.3 | 381.4 | 390.3 | 1.46 | 3.58 | 3.15 |
| Caffeine | def2-qzvpd | jfit | 1218 | 796 | 684.4 | 458.4 | 469.4 | 1.46 | 4.40 | 4.01 |
| Caffeine | def2-qzvpp | jfit | 1098 | 796 | 582.6 | 381.5 | 394.8 | 1.48 | 3.58 | 3.15 |
| Caffeine | def2-qzvppd | jfit | 1218 | 796 | 692.6 | 465.3 | 468.0 | 1.48 | 4.40 | 4.01 |
| Caffeine | def2-svp | jfit | 246 | 796 | 24.5 | 23.9 | 24.2 | 1.01 | 0.18 | 0.18 |
| Caffeine | def2-svpd | jfit | 366 | 796 | 43.8 | 44.5 | 45.3 | 0.97 | 0.40 | 0.40 |
| Caffeine | def2-tzvp | jfit | 494 | 796 | 90.0 | 90.4 | 91.9 | 0.98 | 0.73 | 0.71 |
| Caffeine | def2-tzvpd | jfit | 614 | 796 | 123.6 | 128.7 | 129.5 | 0.96 | 1.12 | 1.13 |
| Caffeine | def2-tzvpp | jfit | 574 | 796 | 115.2 | 104.0 | 106.6 | 1.08 | 0.98 | 0.92 |
| Caffeine | def2-tzvppd | jfit | 694 | 796 | 158.8 | 145.1 | 146.6 | 1.08 | 1.43 | 1.39 |
| Caffeine | def2-qzvp | jkfit | 1098 | 1242 | 803.6 | 503.7 | 524.9 | 1.53 | 5.58 | 4.90 |
| Caffeine | def2-qzvpd | jkfit | 1218 | 1242 | 969.2 | 614.5 | 641.5 | 1.51 | 6.87 | 6.24 |
| Caffeine | def2-qzvpp | jkfit | 1098 | 1242 | 832.8 | 505.7 | 511.8 | 1.63 | 5.58 | 4.90 |
| Caffeine | def2-qzvppd | jkfit | 1218 | 1242 | 981.4 | 640.6 | 617.1 | 1.59 | 6.87 | 6.24 |
| Caffeine | def2-svp | jkfit | 246 | 1242 | 33.4 | 31.3 | 31.7 | 1.05 | 0.28 | 0.27 |
| Caffeine | def2-svpd | jkfit | 366 | 1242 | 63.3 | 58.6 | 59.0 | 1.07 | 0.62 | 0.63 |
| Caffeine | def2-tzvp | jkfit | 494 | 1242 | 125.0 | 118.0 | 119.9 | 1.04 | 1.13 | 1.11 |
| Caffeine | def2-tzvpd | jkfit | 614 | 1242 | 173.8 | 170.2 | 178.7 | 0.97 | 1.75 | 1.76 |
| Caffeine | def2-tzvpp | jkfit | 574 | 1242 | 162.0 | 143.3 | 139.0 | 1.17 | 1.53 | 1.43 |
| Caffeine | def2-tzvppd | jkfit | 694 | 1242 | 217.4 | 191.9 | 194.5 | 1.12 | 2.23 | 2.17 |
| tagrisso | cc-pvdz | rifit | 683 | 2534 | 740.3 | 268.1 | 285.9 | 2.59 | 4.41 | 3.24 |
| tagrisso | def2-svp | jfit | 683 | 2176 | 470.3 | 222.2 | 235.0 | 2.00 | 3.79 | 2.62 |
| tagrisso | def2-svpd | jfit | 1010 | 2176 | 1177.9 | 536.5 | 557.9 | 2.11 | 8.28 | 7.08 |
| tagrisso | def2-svp | jkfit | 683 | 3387 | 628.5 | 288.6 | 306.2 | 2.05 | 5.89 | 4.05 |
| c60 | cc-pvdz | rifit | 840 | 3360 | 2250.3 | 797.3 | 845.3 | 2.66 | 8.84 | 8.94 |
| c60 | def2-svp | jfit | 840 | 2940 | 1420.7 | 605.4 | 642.6 | 2.21 | 7.74 | 7.14 |
| c60 | def2-svp | jkfit | 840 | 4500 | 1868.2 | 789.5 | 837.1 | 2.23 | 11.84 | 10.86 |

### What these numbers say

The geometric mean over the 32 combinations is **1.41**, and it divides sharply by
the extent of the molecule rather than by the basis:

| molecule | atoms | cases | x ref |
| --- | --- | --- | --- |
| caffeine | 24 | 25 | 1.23 |
| tagrisso | 70 | 4 | 2.18 |
| c60 | 60 | 3 | 2.36 |

**On caffeine the driver is close to parity**, and five of its twenty five cases are at
or just below it, the lowest at 0.96. A molecule of twenty four atoms has almost
nothing to screen away, so the sparse tensor holds very nearly the dense set and the
comparison becomes kernel against kernel. On the extended molecules the screening has
something to remove and the driver runs at better than twice the reference.

The same explains the thresholds. 1e-12 and 1e-14 differ by one to three per cent
throughout, because at either threshold nearly every atom pair of a compact molecule
survives.

**The memory is a wash: 1.08 times on the geometric mean, at best 1.46.** For
caffeine the two tensors are the same size to two digits. The sparse storage does not
pay here for the same reason the speed does not, and the regime where it should -- an
extended system -- is the one the twelve gigabyte cap on the reference excludes.

### A note on measuring this

The first two runs of this grid gave 1.50 and 1.71 for the geometric mean. Both were
wrong. The reference was timed once, with no warm up, and single shot timings of a
driver that allocates gigabytes moved by as much as 1.73 times between runs on
unchanged code -- six of thirty two cases beyond fifteen per cent. Warming up and
taking the best of several runs brought the reference times down by as much as 2.4
times and gave the 1.41 above. A timing which is not repeated is not a
measurement.


## Inverting the packed metric

The Coulomb metric of a fitting basis has to be inverted, and it is held in the
packed format, as its slow decay leaves no atom pair below the threshold. The
inversion is `packlin::invert`, which takes a packed matrix and returns the
inverted one.

### Why it does not invert in place

The packed layout of `CPackedMatrix` is the lower triangle in row major order, so
the element of row i and column j with j <= i sits at i (i + 1) / 2 + j. Read as
column major, which is what a math library expects, that is exactly the upper
triangle of LAPACK packed storage. The two agree bit for bit, which was checked
by handing our bytes straight to `dpptrf` and `dpptri` and getting the right
inverse back.

So the packed factorizations could be called on our values with no unpacking at
all, and would hold nothing beyond the matrix itself. They are still the wrong
choice. `dpptrf` and `dpptri` are level 2 BLAS and are not threaded, while the
dense `dpotrf` and `dpotri` are blocked and run on every core. The memory the
packed routines save is one dense matrix, which is a gigabyte at the sizes where
it matters and is not worth several times the run. The inversion therefore
expands the matrix, inverts it dense, and packs the result, at a peak of about
twice the dense matrix.

### The math library against Eigen

The math library is what is compiled when `Makefile.setup` sets one for the
platform, which is Accelerate on macOS and is MKL or OpenBLAS elsewhere. Eigen is
compiled in its place when none is set, so a machine without one needs no change
to build.

The two are not close. Times in seconds, on the M4 Max:

| molecule | fitting basis | n | math library | Eigen | ratio |
| --- | --- | ---: | ---: | ---: | ---: |
| Caffeine | def2-jfit | 796 | 0.003 | 0.023 | 7.7 |
| Caffeine | cc-pV5Z-rifit | 3612 | 0.105 | 2.324 | 22.1 |
| tagrisso | def2-jkfit | 3387 | 0.086 | 1.881 | 21.9 |
| tagrisso | cc-pV5Z-rifit | 10144 | 1.966 | 51.478 | 26.2 |
| c60 | def2-jkfit | 4500 | 0.189 | 4.625 | 24.5 |
| c60 | cc-pV5Z-rifit | 11580 | 3.103 | 76.609 | 24.7 |

Eigen sits at twenty to twenty two gigaflops at every size, which is one core.
The math library climbs to about five hundred, which is the machine. The gap is
threading and blocking, not the algorithm: both factorize and both invert. This
is why the math library is the default and Eigen is the fallback rather than the
other way round.

### The whole grid

With the math library. The build column is the two-center Coulomb driver making
the metric, and is there to show that it is not the cost. Counting n cubed for
the factorization and the inversion together:

| molecule | fitting basis | n | packed | build | invert | Gflop/s |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Caffeine | def2-jfit | 796 | 2.4 M | 0.00 | 0.003 | 191 |
| Caffeine | def2-jkfit | 1242 | 5.9 M | 0.00 | 0.007 | 292 |
| Caffeine | cc-pVDZ-rifit | 924 | 3.3 M | 0.00 | 0.004 | 211 |
| Caffeine | cc-pVTZ-rifit | 1434 | 7.8 M | 0.00 | 0.009 | 334 |
| Caffeine | cc-pVQZ-rifit | 2398 | 21.9 M | 0.00 | 0.032 | 428 |
| Caffeine | cc-pV5Z-rifit | 3612 | 49.8 M | 0.01 | 0.105 | 450 |
| tagrisso | def2-jfit | 2176 | 18.1 M | 0.00 | 0.023 | 446 |
| tagrisso | def2-jkfit | 3387 | 43.8 M | 0.00 | 0.086 | 454 |
| tagrisso | cc-pVDZ-rifit | 2534 | 24.5 M | 0.00 | 0.038 | 433 |
| tagrisso | cc-pVTZ-rifit | 3987 | 60.7 M | 0.00 | 0.136 | 465 |
| tagrisso | cc-pVQZ-rifit | 6699 | 171 M | 0.02 | 0.608 | 495 |
| tagrisso | cc-pV5Z-rifit | 10144 | 393 M | 0.06 | 1.966 | 531 |
| c60 | def2-jfit | 2940 | 33.0 M | 0.00 | 0.055 | 458 |
| c60 | def2-jkfit | 4500 | 77.3 M | 0.01 | 0.189 | 482 |
| c60 | cc-pVDZ-rifit | 3360 | 43.1 M | 0.00 | 0.079 | 478 |
| c60 | cc-pVTZ-rifit | 4860 | 90.1 M | 0.01 | 0.236 | 487 |
| c60 | cc-pVQZ-rifit | 7920 | 239 M | 0.03 | 0.941 | 528 |
| c60 | cc-pV5Z-rifit | 11580 | 512 M | 0.09 | 3.103 | 500 |
| taxol | def2-jfit | 3528 | 47.5 M | 0.01 | 0.095 | 464 |
| taxol | def2-jkfit | 5489 | 115 M | 0.01 | 0.350 | 472 |
| taxol | cc-pVDZ-rifit | 4102 | 64.2 M | 0.00 | 0.157 | 439 |
| taxol | cc-pVTZ-rifit | 6411 | 157 M | 0.01 | 0.559 | 471 |
| taxol | cc-pVQZ-rifit | 10747 | 441 M | 0.04 | 2.584 | 480 |
| taxol | cc-pV5Z-rifit | 16232 | 1005 M | 0.14 | 7.989 | 535 |
| crambin | def2-jfit | 19500 | 1451 M | 0.14 | 14.40 | 515 |
| crambin | def2-jkfit | 30751 | 3607 M | 0.39 | 61.54 | 473 |
| crambin | cc-pVDZ-rifit | 22842 | 1990 M | 0.19 | 23.97 | 497 |
| crambin | cc-pVTZ-rifit | 36183 | 4994 M | 0.37 | 99.23 | 477 |

**The throughput plateaus near five hundred gigaflops from about n of 2500 and
holds it to 36183.** Nothing degrades at scale. The small cases are below the
plateau because a factorization of a few hundred rows cannot fill the machine,
not because anything is wrong with them.

**Building the metric is not the cost.** For crambin with jkfit it is 0.39
seconds against 61.5 for the inversion, a factor of 160. The integrals are cheap
and the cubic step is what is paid for.

The two largest crambin cases were skipped, and on memory rather than on time.
cc-pVQZ-rifit is n of 60645, whose peak of two dense matrices is 54.8 gigabytes
against the 36 of the machine, and cc-pV5Z-rifit is 126 gigabytes.

### The residual, and what it says about the metric

The inversion was checked as the largest element of A A inverse minus the
identity, wherever the dense matrices of the check were small enough to be cheap:

| molecule | fitting basis | n | max abs residual |
| --- | --- | ---: | ---: |
| Caffeine | def2-jfit | 796 | 4.0e-10 |
| Caffeine | def2-jkfit | 1242 | 2.2e-07 |
| Caffeine | cc-pVTZ-rifit | 1434 | 2.3e-09 |
| Caffeine | cc-pV5Z-rifit | 3612 | 8.2e-07 |
| tagrisso | def2-jfit | 2176 | 6.1e-10 |
| tagrisso | def2-jkfit | 3387 | 9.3e-06 |
| tagrisso | cc-pVTZ-rifit | 3987 | 8.7e-09 |
| c60 | def2-jfit | 2940 | 2.6e-09 |
| taxol | def2-jfit | 3528 | 1.7e-09 |

**The residual tracks the basis, not the size.** tagrisso with jfit at n of 2176
gives 6e-10, and with jkfit at n of 3387 gives 9e-06, four orders apart at
comparable size. The residual of a solved system is about the condition number
times the machine epsilon, so 9e-06 says the jkfit metric of tagrisso has a
condition number near 4e10. That is the fitting basis being close to linearly
dependent, which is a property of the basis and not of the inversion. It is worth
knowing before an inverse of a jkfit or a high zeta RIFIT metric is fed to
anything that cares about its accuracy.

### What this residual cannot be used for

The number above is the largest element of A A inverse minus the identity. It is
a fair measure of one inversion against another **only when both results are
constrained the same way**, and it is worthless across that line. Two traps,
both of which were walked into while writing this section.

**Do not compare a symmetric inverse against an unconstrained one.** The packed
format stores one triangle, so what comes back is symmetric to the last bit. The
LU route of the reference returns a matrix which is not, and the residual rewards
it for that: the error is free to be antisymmetric, and antisymmetric error is
invisible to the identity it is measured against. Taking the same LU inverse and
symmetrizing it, which is all that storing it in the packed format would do,
moves its residual from 3.1e-08 to 4.5e-01 on caffeine with cc-pV5Z-RIFIT. Seven
orders, from the same numbers, for no other reason than being made symmetric. A
comparison run this way reported the packed inversion as ten to two hundred
times less accurate than numpy. On forward error, against a refined reference,
it is within one to five times of it and better on some cases.

**Do not read small differences as quality.** The condition number of these
metrics reaches 1e10, and the residual multiplies whatever perturbs the inverse
by roughly that. A Newton step written two ways which differ only in the order
the library accumulates -- the same matrix mathematically, and agreeing to
1.6e-11 where the two were compared element by element -- gives residuals of
3.0e-07 and 5.0e-01 on tagrisso with jkfit. A four order spread from an
accumulation order. Anything the residual says at that scale is the condition
number talking, not the algorithm.

What the residual is good for is what the table above uses it for: the order of
magnitude, against the conditioning of the same matrix, for results produced the
same way. For comparing two methods, use the forward error against a reference
refined until it stops moving, and check that the reference is converged further
than the difference being claimed.

Why the plain Cholesky inverse is hard to beat here is worth recording too. With
A about L L transposed and X about L inverse transposed times L inverse, the
product A X cancels structurally rather than numerically, which is why its
residual sits below the condition number times the epsilon. Refining it breaks
that cancellation. In working precision there is nothing to gain: iterative
refinement lowers the backward error only when the residual is formed in higher
precision than the working one, and without that the inversion is already at the
floor.

### The Cholesky and the fallback

The metrics are positive definite and are inverted through their Cholesky
factorization. Near linear dependence can still push one numerically indefinite,
and rather than fail there the inversion falls back to the Bunch-Kaufman
factorization and prints a warning. Nothing in the grid above reached the
fallback; the test suite reaches it with a constructed indefinite matrix, and
gets the right inverse from it.

## The Coulomb matrix of the resolution of the identity

The chain is four steps. The metric (p|J|q) of the fitting basis is factorized as
L L transposed and its factor inverted; the B vectors are formed as B(q)_ij = sum
over p of Linv_qp (ij|p); a density is contracted into Y(q) = sum over i and j of
B(q)_ij D_ij; and the Coulomb matrix is F_ij = sum over q of B(q)_ij Y(q).

### Why the factor and not the inverse

The closing step is what fixes which matrix builds the B vectors. Writing M for
that matrix,

    sum over q of B(q)_ij B(q)_kl = sum over p and t of (ij|p) [M transposed M]_pt (kl|t)

and the resolution of the identity needs the bracket to be the inverse of the
metric. With the metric equal to L L transposed its inverse is L inverted
transposed times L inverted, so M is L inverted. **The plain inverse of the metric
is not a substitute for it**, and not by a little:

| combination | closes with L inverted | closes with the inverse |
| --- | ---: | ---: |
| (ss\|s) | 7.4e-16 | 8.97e-01 |
| (ps\|p) | 4.8e-16 | 7.09e-01 |
| (pp\|d) | 8.0e-16 | 5.41e-01 |
| (ds\|p) | 7.4e-16 | 7.12e-01 |

The right matrix closes to the last bit and the wrong one is wrong by most of the
answer. The test suite asserts both halves of that table, so the two cannot be
confused again.

Inverting the factor is also the cheaper of the two. The factorization costs a
third of the cube of the dimensions and the inversion of the factor another
third, against a further third to form the whole inverse:

| molecule | fitting basis | n | invert | cholesky_inverse | ratio |
| --- | --- | ---: | ---: | ---: | ---: |
| Caffeine | def2-jfit | 796 | 0.003 | 0.002 | 1.35 |
| Caffeine | cc-pV5Z-RIFIT | 3612 | 0.104 | 0.072 | 1.44 |
| tagrisso | def2-jkfit | 3387 | 0.087 | 0.060 | 1.44 |
| c60 | cc-pVQZ-RIFIT | 7920 | 0.962 | 0.681 | 1.41 |
| taxol | cc-pV5Z-RIFIT | 16232 | 7.980 | 5.628 | 1.42 |
| crambin | def2-jfit | 19500 | 14.480 | 10.234 | 1.41 |
| crambin | def2-jkfit | 30751 | 59.980 | 42.291 | 1.42 |

The measured 1.42 is the ratio of the flop counts, which is 1.5, less what the
triangular inversion loses to the shape of its blocks.

### What the screening is worth

The B vectors are screened on the pair of atomic orbitals, as the three-center
integrals are, and are dense in q, as the metric is. Their memory is therefore the
surviving pairs times the auxiliary basis. Against the dense triangle of the
reference, with the overlap bound at 1e-12:

| molecule | basis | nao | naux | dense triangle | screened | kept |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Caffeine | def2-svp | 246 | 796 | 30381 | 26712 | 87.9% |
| Caffeine | def2-tzvp | 494 | 796 | 122265 | 107767 | 88.1% |
| tagrisso | def2-svp | 683 | 2176 | 233586 | 132082 | 56.5% |
| c60 | def2-svp | 840 | 2940 | 353220 | 283430 | 80.2% |
| taxol | def2-svp | 1099 | 3528 | 604450 | 281682 | 46.6% |
| crambin | def2-svp | 6177 | 19500 | 19080753 | 2934925 | 15.4% |

Screening is worth a tenth on a compact molecule and six and a half times on
crambin. It pays where it matters and nowhere else, which is the same shape the
timings below have.

### Forming the B vectors

Against the two routes of the reference which form the same thing,
compute_screened_bq_vectors and compute_bq_vectors of CRIJKFockDriver, at 1e-12:

| molecule | basis | nao | naux | simd | ref screened | ref full | vs screened | vs full |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Caffeine | def2-svp | 246 | 796 | 0.343 | 0.394 | 0.421 | 1.15 | 1.23 |
| Caffeine | def2-tzvp | 494 | 796 | 1.657 | 1.777 | 2.152 | 1.07 | 1.30 |
| tagrisso | def2-svp | 683 | 2176 | 4.548 | 15.917 | 41.580 | 3.50 | 9.14 |
| c60 | def2-svp | 840 | 2940 | 10.163 | 96.498 | 141.301 | 9.50 | 13.90 |
| taxol | def2-svp | 1099 | 3528 | 21.192 | 130.246 | 444.399 | 6.15 | 20.97 |

**The gap grows with the molecule**, from a tenth on caffeine to six and ten times
on taxol and c60. The reference walks the auxiliary functions one at a time and
adds one scaled vector per surviving pair of them, which is level 1 work and as
many calls as the square of the auxiliary basis. The new driver contracts a block
at a time with one matrix product per pair of angular components, which is level
3. As the fitting basis grows, level 1 against level 3 is the whole story.

### The Coulomb matrix

Against CRIFockDriver, whose route keeps the raw integrals and applies the metric
inside every build. Setup is the B vectors for the new driver and prepare_buffers
for the reference; the Coulomb column is one build for one density:

| molecule | basis | nao | naux | setup | ref setup | F | ref F | speedup | B | ref |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Caffeine | def2-svp | 246 | 796 | 0.346 | 0.029 | 0.0051 | 0.0045 | 0.89 | 0.17 GB | 0.18 GB |
| Caffeine | def2-tzvp | 494 | 796 | 1.494 | 0.116 | 0.0200 | 0.0178 | 0.89 | 0.70 GB | 0.73 GB |
| tagrisso | def2-svp | 683 | 2176 | 4.450 | 0.593 | 0.0306 | 0.0955 | 3.12 | 2.49 GB | 3.79 GB |
| c60 | def2-svp | 840 | 2940 | 10.328 | 1.705 | 0.0565 | 0.1975 | 3.49 | 6.90 GB | 7.74 GB |
| taxol | def2-svp | 1099 | 3528 | 20.992 | 3.274 | 0.0933 | 0.4114 | 4.41 | 8.76 GB | 15.89 GB |

**Three to four times faster per build on the real systems**, and level with the
reference on caffeine, where nothing is large enough for the fixed costs to
disappear. The memory is 1.8 times better on taxol, as the reference stores a
dense triangle of atomic orbitals for every auxiliary function and the new driver
only the pairs which survive.

The intermediate step, the Y vector alone, is close to level: 1.0, 1.1 and 1.5 on
tagrisso, c60 and taxol, and slower than the reference on caffeine at 0.3. Both
read a structure of the same size and are bound by the memory rather than by the
arithmetic. The new driver loses on the small molecules because it walks a blocked
sparse structure with a short inner loop per pair of angular components, while the
reference runs one long contiguous loop per auxiliary function.

### Where the two routes cross

The new driver spends its metric transform once, in the B vectors, and the
reference spends it inside every build. Their totals for n densities cross where
the difference of the setups equals n times the difference of the builds, which
for taxol is

    (20.99 - 3.27) / (0.4114 - 0.0933) = 56 densities

Above a field calculation and well inside a response one. Assembling the Coulomb
matrix is what moves that number: on the Y vector alone the same molecule crosses
at about 780, and adding the step which uses the B vectors a second time brings it
down by a factor of fourteen.

### A note on which setup is being compared

There are two reference routes and they do not do the same work, which makes
"setup" ambiguous unless it is said which one is meant. On taxol with def2-svp:

| what it produces | route | time |
| --- | --- | ---: |
| B vectors, metric applied | the new driver | 21.0 |
| B vectors, metric applied | CRIJKFockDriver::compute_screened_bq_vectors | 130.2 |
| raw three-center integrals, no metric | CRIFockDriver::prepare_buffers | 3.3 |

Against the route which produces the same thing the new driver is six times
faster. Against the route which only computes the integrals it is six times
slower, and that is not a comparison of like with like: prepare_buffers never
applies the metric, and CRIFockDriver pays for it again in every build. The two
ratios both land near six by coincidence, which is exactly how a comparison of
this kind gets reported the wrong way round.

### Accuracy against the reference

The Coulomb matrix agrees with the one CRIFockDriver builds to between 1.5e-10 and
2.9e-10 relative across the set, and the Y vector to between 7.8e-11 and 9.0e-10.
The two routes sum in different orders and take their three-center integrals from
different drivers, and the metrics of the fitting bases have condition numbers
which reach 1e10, so agreement of this order is what a correct implementation
looks like rather than a sign of one.

Turning the screening from 1e-12 down to nothing moves the Coulomb matrix by 1e-13
relative, so the screening is not what limits the agreement.

### On reproducibility

The Y vector and the Coulomb matrix are summed by the threads into vectors of
their own, which are added up in the order of the threads. The order the blocks
reach the threads is not fixed, though, as they are handed out dynamically, so
which of them a thread sums varies between runs and the last bit of the total
varies with it. The measured spread over five runs is one unit in the last place.
Making it exact would need a static schedule and the imbalance which comes with
it, and is not worth that.

## The W matrices, and the matrix unit of the machine

The transformation of one index of the B vectors into the molecular orbitals,
W(q)_is = sum over r of B(q)_ir C_rs, turned out to be where a calculation spends
most of its time, and the way it was written could not reach the machine. This is
how that was found and what was done about it.

### Where a calculation spends its time

Tagrisso in def2-svp against def2-universal-jkfit, the whole calculation:

| phase | time | share |
| --- | ---: | ---: |
| setup: the B vectors and the metric | 10.79 | 13.3% |
| of which the two-center integrals | 0.00 | 0.0% |
| of which inverting the Cholesky factor | 0.06 | 0.1% |
| the Fock builds, twenty three of them | 61.83 | 76.5% |
| everything else | 8.26 | 10.2% |

and one Fock build:

| | time | share of the build |
| --- | ---: | ---: |
| Coulomb and exchange | 2.671 | 100% |
| the Coulomb alone | 0.050 | 1.9% |
| the exchange, by difference | 2.620 | 98.1% |

**The exchange is the calculation.** It is 98 per cent of a build and three
quarters of the run, and the Coulomb matrix, which the earlier sections are about,
is under two per cent of a build. Within the exchange the W matrices are 84 per
cent and the rank k update which follows them is 16.

The metric is nothing at all: six hundredths of a second to invert the Cholesky
factor. The choice between that and the inverted square root, which the earlier
section measures as a factor of eight to twenty seven, cannot matter to a
calculation of this shape.

### What the sum was reaching, and what the machine has

The transformation as first written walks the values of the B vectors and adds a
scaled row of the coefficients for each of them. Over the threads:

| threads | time | gigaflops per second |
| ---: | ---: | ---: |
| 1 | 20.00 | 13.7 |
| 2 | 10.28 | 26.6 |
| 4 | 5.38 | 50.8 |
| 8 | 2.98 | 91.9 |
| 12 | 2.34 | 117.2 |
| 16 | 2.19 | 125.0 |

Nine times on sixteen cores, and not flat, so the sum is not held up by the memory.
Against that, what a matrix product reaches on the same machine:

| threads | shape | gigaflops per second |
| ---: | --- | ---: |
| 1 | 683 x 683 x 133 | 417.7 |
| 1 | 2000 cubed | 415.6 |
| 16 | 2000 cubed | 887.9 |
| 16 | 683 x 3387 x 133 | 747.2 |

**One thread of a matrix product is three times the whole of the sum on sixteen.**
That is the matrix unit of the processor, which the library reaches through its
matrix product and which a loop of the compiler does not reach at any number of
threads. Per core the two are 13.7 against 417.7, a factor of thirty.

So the question is not how to write a better sum. It is whether the work can be
put into a matrix product at all.

### The trade

The B vectors are not sparse. Of the square that one auxiliary function expands
into, better than half is filled, and for a compact molecule better than nine
tenths. Expanding them and handing the square to a matrix product is about one and
a half times the arithmetic of walking their values, and the matrix unit runs it
several times faster than the sum runs the smaller amount. Measured, against
def2-universal-jkfit throughout:

| molecule | basis | nao | naux | norb | sum | product | gain | sum Gf/s | product Gf/s |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Caffeine | def2-svp | 246 | 1242 | 51 | 0.06 | 0.03 | 2.32 | 118.8 | 286.2 |
| tagrisso | def2-svp | 683 | 3387 | 133 | 2.05 | 0.66 | 3.13 | 133.4 | 640.7 |
| tagrisso | def2-svpd | 1010 | 3387 | 133 | 6.16 | 1.60 | 3.86 | 123.7 | 575.7 |
| c60 | def2-svp | 840 | 4500 | 180 | 7.67 | 1.56 | 4.92 | 131.9 | 734.1 |
| taxol | def2-svp | 1099 | 5489 | 223 | 12.61 | 4.41 | 2.86 | 128.0 | 670.1 |

and over the basis sets of one molecule, where the orbitals stay at fifty one
however large the basis grows:

| basis | nao | sum | product | gain | product Gf/s |
| --- | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 246 | 0.06 | 0.03 | 2.36 | 284.2 |
| def2-svpd | 366 | 0.14 | 0.05 | 2.66 | 312.9 |
| def2-tzvp | 494 | 0.25 | 0.10 | 2.57 | 322.0 |
| def2-tzvpd | 614 | 0.39 | 0.15 | 2.56 | 313.2 |
| def2-qzvp | 1098 | 1.01 | 0.63 | 1.61 | 242.8 |
| def2-qzvpd | 1218 | 1.33 | 0.76 | 1.74 | 246.1 |

**The product won every one of the eleven.** The sum sits at 119 to 133 gigaflops
whatever it is given, which is the rate of the loop and not of the problem. The
product reaches 640 to 734 where the orbitals are many, and falls to 243 at
quadruple zeta, where fifty one orbitals make the third dimension of the product
too thin for the matrix unit. Caffeine has fifty one occupied orbitals whatever
basis it is given, so that is a property of the molecule and not of the basis.

### The density, and the threshold which is not set

Which form is taken follows from how dense the B vectors are, and the density is
the filled places of the square over its size. An off-diagonal pair of atoms is
held once and fills two places; a diagonal pair is held with the basis functions
of both sides and fills one. Counting the values twice over, which is the obvious
thing to write, counts the diagonal pairs twice and puts the answer above one:

| molecule | basis | counting twice | counted properly |
| --- | --- | ---: | ---: |
| Caffeine | def2-svp | 96.1% | 91.2% |
| Caffeine | def2-svpd | 100.9% | 96.0% |
| Caffeine | def2-tzvpd | 99.7% | 94.3% |
| Caffeine | def2-qzvp | 85.4% | 80.9% |
| tagrisso | def2-svp | 65.1% | 63.4% |

The density columns of the two tables above are of the first kind, as they were
taken before this was found. They are two to five points high and none of the
timings depends on them.

**The threshold is zero, which is to say that the product is always taken.** The
density at which the two forms change places was never reached: every case which
fits in the memory of this machine is denser than half, and the product won all of
them. A threshold naming a density would be a guess dressed as a measurement. A
caller which meets a sparse enough set of B vectors can raise it, and the two
forms are held to one another by a test so that the sum does not rot while it is
not the default.

## The field calculation through the resolution of the identity

The chain is reachable from the input of a closed shell calculation, as the
alternative path through the Fock build which ri_jk_simd selects. This is what it
does to a whole calculation rather than to one matrix.

### Caffeine and tagrisso, whole calculations

*Superseded. Both molecules at Hartree-Fock and B3LYP, four builds, with the B
vectors timed apart from the builds:
`benchmarks/data/scf/2026-09-15_m4max_caffeine.md` and
`benchmarks/data/scf/2026-09-15_m4max_tagrisso.md`. The tables which stood here
compared a full build against the two resolution of the identity routes over the
def2 sets; they crossed two runs and are replaced by measurements which do not.*

### What the memory check did

The diffuse row did not run at the default budget. The B vectors need 10.67
gigabytes and the driver was given 10.57, which is half of what was free, and it
refused:

    RIJKFockDriver.prepare: The B vectors need 10.674862 GB and the budget is 10.568748 GB

It refused before computing the integrals, from the sparsity pattern alone, rather
than at the end of the ten minutes it would have spent forming them. The row above
was taken with ri_memory_budget set to 24 gigabytes, which the machine has.

**The default is too tight for a case of this shape.** Half of what is free is a
reasonable guard when the B vectors are one term among several, and a poor one
when they are the term which dominates, as they are here: they alone are a third
of the memory of the machine. A default of what is free less a fixed reserve would
have taken this run. The check itself did what it exists to do, which is to say
which number is too small rather than to die in the allocator.

### On the metric of these bases

The Cholesky factorization succeeded for every basis of every molecule here,
including the diffuse ones, so the fallback to the inverted square root was never
taken and no warning was printed. The metrics of the universal fitting set are
well enough conditioned for the cheaper factorization at these sizes. The fallback
is therefore still covered by constructed matrices alone and not by a calculation.

## Forming the B vectors, and the depth of its products

Once the W matrices were formed by a matrix product, the setup became the largest
fixed cost of a calculation: thirty seconds of a hundred and fifty for tagrisso in
def2-svpd, a fifth of the run. This is what it was made of and what was done to it.

### What the setup is

| stage | time | share |
| --- | ---: | ---: |
| the two-center integrals | 0.007 | 0.0% |
| inverting the Cholesky factor | 0.070 | 0.2% |
| the three-center integrals | 0.75 | 2.6% |
| contracting them with the metric | 28.33 | 97.4% |

**The integrals are not the setup; the contraction is.** Forming every three-center
integral of the molecule takes three quarters of a second, and turning them into
the B vectors took twenty eight. Each of the 1.43 billion values of B is a sum over
all 3387 auxiliary functions, which is 9.7 teraflops, and it was running at 343
gigaflops per second where a large product reaches seven hundred.

### Why it was slow

The contraction was one matrix product for each block of atom pairs, each basis
function on the auxiliary side of the output, **each basis function on the
auxiliary side of the input**, each combination of basis functions and each angular
component. The auxiliary side of tagrisso against def2-universal-jkfit is four
groups of atoms carrying eighty one basis functions between them, and the products
were of this shape:

| | rows | depth | columns | flops |
| --- | ---: | ---: | ---: | ---: |
| one product | 42 | 42 | 31 | 107 thousand |

A hundred and seven thousand flops is a third of a microsecond of arithmetic, which
is the same order as the call itself. The depth of forty two is the functions of
one group, and the machine was given eighty one slivers where it wanted one
product.

### The gather

Every one of those eighty one products adds into the same rows of the output, so
they are one sum over the whole auxiliary basis broken into pieces. Gathering the
integrals of all the groups into one buffer first makes the depth the whole
auxiliary basis:

| | rows | depth | columns | flops |
| --- | ---: | ---: | ---: | ---: |
| gathered | 42 | 3387 | 31 | 9 million |

The buffer is the auxiliary basis by the atom pairs of the block, which is under a
megabyte, and it is written once and read eighty one times. Against the arithmetic
it enables that is some eight hundred and fifty flops for every byte written, which
is why it costs nothing worth measuring. A group whose block is absent, and the
atom pairs a group keeps fewer of than the widest, leave zeros in the buffer, and a
zero adds what the sum of that group would have added.

| | before | after | gain |
| --- | ---: | ---: | ---: |
| the contraction | 343 Gflop/s | 554 Gflop/s | 1.62 |
| forming the B vectors | 29.08 | 18.26 | 1.59 |

**The products are still thin.** Thirty one columns is the atom pairs of a block,
and no gathering of the auxiliary side changes that, which is why 554 and not the
seven hundred and fifty a squarer product reaches.

### What it is worth

| mode | before | after | gain |
| --- | ---: | ---: | ---: |
| Hartree-Fock | 150.22 | 138.75 | 1.08 |
| B3LYP | 212.07 | 199.91 | 1.06 |

**A sixth off the setup is a fifteenth off the calculation**, because the setup was
a fifth of it. This is the smallest of the three changes these benchmarks record:
the matrix product form of the W matrices was worth 1.8 times on a run and this is
worth 1.08. It is worth more to a calculation which forms the B vectors once and
uses them many times, as a geometry optimization or a response calculation does,
than to a single field calculation.

The energies are unchanged to all ten digits in both modes, and the Fock builds and
the rest of the calculation are unchanged, so the change is where it was meant to
be and nowhere else.

## The two modes of a hybrid

Tagrisso in def2-svpd, Hartree-Fock against B3LYP, with the setup above:

| mode | total | setup | Fock builds | share | per build | the rest | share | iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Hartree-Fock | 138.75 | 18.60 | 60.21 | 43.4% | 2.509 | 59.93 | 43.2% | 24 |
| B3LYP | 199.91 | 18.60 | 58.14 | 29.1% | 2.528 | 123.17 | 61.6% | 23 |

and one Fock build of each:

| | time |
| --- | ---: |
| Hartree-Fock, exchange scaled by one | 2.510 |
| B3LYP, exchange scaled by 0.20 | 2.498 |
| the Coulomb alone, exchange scaled by zero | 0.115 |
| the exchange, by difference | 2.395 |

**Scaling the exchange costs nothing.** A build of B3LYP and a build of
Hartree-Fock are the same time to a twentieth of a per cent, because the whole of
the W matrices and the whole of the exchange are formed and the fraction is applied
to the result. A hybrid pays the full price of the exact exchange whatever fraction
of it the functional asks for.

**What a functional adds is its quadrature, and for B3LYP that is the calculation.**
The Fock builds of the two modes are the same, and everything else doubles from
sixty seconds to a hundred and twenty three, which is sixty two per cent of the
run. Nothing in the resolution of the identity is the bottleneck of a hybrid
calculation of this size: the exchange is 2.4 seconds a build against some 2.8
seconds a step of quadrature.

## The exchange build, once both of its halves are products

The matrix product form of the W matrices changed which half of the exchange build
costs what. Tagrisso in def2-svpd, with 3387 auxiliary functions and 133 occupied
orbitals:

| stage | time | share | teraflops | gigaflops per second |
| --- | ---: | ---: | ---: | ---: |
| W, expanding and multiplying | 1.63 | 66.8% | 0.92 | 563.0 |
| the exchange, a rank k update | 0.81 | 33.2% | 0.46 | 565.5 |
| the build | 2.44 | 100% | 1.38 | |

against what it was when the W matrices were formed by walking the values of the B
vectors, where W was 84 per cent of the build at 119 gigaflops per second and the
rank k update 16 per cent at 480.

**The two halves now run at the same rate**, 563 against 565, within half a per
cent of one another. **The two to one split of the time is the arithmetic and not
an inefficiency.** W is the product of a square by the orbitals, which is two
flops for every element of the square and every orbital; the rank k update writes
one triangle and is half of that. Twice the arithmetic at the same rate is twice
the time, which is what the table shows.

### What is left in it

Both halves are at 563 against the 750 a product of these dimensions reaches and
the 888 of a square one. Two things account for the difference and neither is worth
much.

The occupied orbitals are 133, which is a thin third dimension for the matrix unit
and is a property of the molecule rather than of anything that can be written
differently.

The squares the B vectors are expanded into are zeroed and scattered into once for
every auxiliary function, which is 25.7 gigabytes of writing. At a hundred
gigabytes a second that is 0.28 of the 1.63 seconds above, a sixth of the W matrices
and an eighth of the build. It cannot be removed, as the squares have to be built
somehow, but it could be trimmed: the auxiliary functions of one group write the
same places, so clearing those places rather than the whole square would do.

**The whole of that is worth about a twentieth of a calculation**, and the exchange
build is left as it is. What is left in a Hartree-Fock run of this size is not the
exchange: the setup is an eighth of it, the Fock builds are a little over two
fifths, and everything else, which is not the resolution of the identity at all, is
the other two fifths. For B3LYP the quadrature is three fifths of the run on its
own.

## Caffeine, the routes side by side and against the build which makes none

*Superseded. Caffeine at Hartree-Fock and B3LYP, four builds, eight orbital and fitting pairs, with the B vectors timed apart from the builds:
`benchmarks/data/scf/2026-09-15_m4max_caffeine.md`.*

## The larger molecules, the two routes side by side

c60 and taxol at Hartree-Fock in def2-svp against def2-universal-jkfit, both routes
one after the other in the same process. Not measured again since; the caffeine and
tagrisso tables they were taken beside have been.

| molecule | nao | naux | occupied | conventional | simd | gain | iterations | energy, conventional | energy, simd |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| c60 | 840 | 4500 | 180 | 391.53 | 91.28 | 4.29 | 24 | -2269.9104513062 | -2269.9104513050 |
| taxol | 1099 | 5489 | 223 | 623.92 | 207.34 | 3.01 | 23 | -2907.5505875335 | -2907.5505875339 |

Both converged, in the same number of iterations by either route, to energies which
agree to the ninth decimal. The last digit or two differ, which is the order the
arithmetic is summed in over twenty odd iterations.

**These are the largest gains of any calculation in this file.** Taken with the
tables above, and with the caffeine and tagrisso rows as they stood when this was
written rather than as `benchmarks/data/scf` now has them:

| molecule | auxiliary functions | occupied orbitals | simd against conventional |
| --- | ---: | ---: | ---: |
| caffeine | 1242 | 51 | 1.89 to 8.17 |
| tagrisso | 3387 | 133 | 2.73 to 4.23 |
| taxol | 5489 | 223 | 3.01 |
| c60 | 4500 | 180 | 4.29 |

The advantage follows the fitting set and the occupied orbitals together, which is
what the benchmarks of the exchange said it should: those two are the rows and the
third dimension of the products the exchange is built from, and a product which is
thin in either of them leaves the matrix unit idle. Caffeine has fifty one occupied
orbitals whatever basis it is given, and that is what holds it to the bottom of the
table however large its orbital basis grows.

c60 gains more than taxol although its fitting set is smaller. The likely reason is
that its sixty identical atoms in a cage make the B vectors denser, 88.5 per cent
against 54.6 for taxol by the counting of the section on the W matrices, and the
denser they are the more the expanded form of the transformation pays. That is not
separated by these runs and is offered as the likely reason rather than a measured
one.

## A transition metal complex

A dinuclear copper guanidinate, C18H40Cu2N6, built from the SMILES

    CC(C)N1C(N(C)C)=[N+](C(C)C)[Cu-]N(C(C)C)C(N(C)C)=[N+](C(C)C)[Cu-]1

which is two copper centers bridged by two guanidinate ligands in an eight
membered ring, each ligand carrying a dimethylamino group on its central carbon
and an isopropyl on each of its two ring nitrogens. Sixty six atoms, two hundred
and forty eight electrons, a neutral closed shell singlet: the SMILES writes the
ring nitrogens as positive and the metals as negative, which nets to nothing and
reads chemically as copper in its first oxidation state, with a filled d shell.

Hartree-Fock against def2-universal-jkfit with 3060 auxiliary functions, both
routes one after the other in the same process.

| basis | nao | conventional | simd | gain | iterations | energy, conventional | energy, simd |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 598 | 77.76 | 26.48 | 2.94 | 23 | -4307.9982685253 | -4307.9982685252 |
| def2-tzvp | 1074 | 374.85 | 97.20 | 3.86 | 24 | -4309.5477190075 | -4309.5477190075 |

Both converged, in the same number of iterations by either route, to energies which
agree to the tenth decimal.

**The gain barely moves as the orbital basis nearly doubles**, 2.13 to 2.04, which
is the same thing caffeine showed over its own basis sets. What differs is where
it sits: caffeine fell from 1.57 to 1.15 over a comparable range and this holds
near two. The fitting set is 3060 functions in both rows and the occupied orbitals
are 124 in both, and 124 is enough to keep the products of the exchange from being
thin, where the fifty one of caffeine is not.

Taken with everything else in this file, the two axes separate: **the occupied
orbitals set the level, the orbital basis hardly matters, and the fitting set moves
it further.**

| molecule | auxiliary functions | occupied orbitals | simd against conventional |
| --- | ---: | ---: | ---: |
| caffeine | 1242 | 51 | 1.89 to 8.17 |
| copper guanidinate | 3060 | 124 | 2.94 to 3.86 |
| tagrisso | 3387 | 133 | 2.73 to 4.23 |
| taxol | 5489 | 223 | 3.01 |
| c60 | 4500 | 180 | 4.29 |

### On the geometry, and on what would not fit

The structure is a distance geometry embedding and is **not optimized at any level
of theory**. Its copper to nitrogen distances are 1.885 to 1.908 angstrom and its
shortest contact of any kind is a carbon to hydrogen bond, so it is a reasonable
structure to time a calculation on, and it is not a structure to draw any chemistry
from. The copper to copper separation in particular, which is the interesting
quantity in these compounds, is whatever the embedding produced.

Building it through the force field which the SMILES reader applies gave a copper
to nitrogen distance of 0.776 angstrom, shorter than a carbon to hydrogen bond,
because the force field has no parameters for copper and left the metals
unconstrained while it pulled everything else into place. The embedding without
that step is what the table above was run on. A geometry from a structure builder
is worth looking at before it is used.

A larger complex was tried first and did not fit: tetrakis(triphenylphosphine)
copper as its cation, a hundred and thirty seven atoms, needs 1411 orbital
functions and 7256 auxiliary ones in def2-svp, which is 34.4 gigabytes of B vectors
against the 36 of this machine, and 53.9 gigabytes for the dense triangle the
conventional route keeps. Neither route can take it here. That brackets what this
machine holds at around a hundred atoms of that composition.

### The metric of a fitting set which carries a metal

The Cholesky factorization succeeded for both bases, so the fallback to the
inverted square root was not taken. The compact functions a first row transition
metal puts into a fitting set did not make its metric indefinite at either zeta
level. The fallback is still exercised by constructed matrices alone.

## A metal oxide cluster, where the occupied orbitals are many

A Ti15O30 cluster of forty five atoms, cut from the bulk and left unrelaxed: one
titanium at its center keeps the six neighbours of the bulk, the rest carry four or
five, and most of its oxygens bridge two titaniums where the bulk gives them three.
Five hundred and seventy electrons and a neutral closed shell singlet.

This is the system which pushes the occupied orbitals hardest. Titanium puts many
electrons behind few basis functions, so the occupied orbitals are 285 of 885,
against about a fifth for every organic molecule in this file.

| basis | nao | naux | occupied over nao | B vectors | conventional triangle |
| --- | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 885 | 6270 | 0.32 | 15.17 GB | 18.31 GB |
| def2-tzvp | 1605 | 6270 | 0.18 | 46.05 GB | 60.21 GB |

**Only the single zeta basis fits.** The triple zeta one needs forty six gigabytes
of B vectors against the thirty six of this machine, and sixty for the triangle the
conventional route keeps.

| route | time | iterations | energy |
| --- | ---: | ---: | ---: |
| RI-JK conventional | 792.01 | 23 | -14971.0615474580 |
| RI-JK simd | 193.43 | 23 | -14971.0615474553 |

**Three times the conventional route**, and both converged in twenty three
iterations, which a cut cluster with undercoordinated surface atoms was not certain
to do.

### What it says about which axis matters

| molecule | auxiliary functions | occupied orbitals | occupied over nao | simd against conventional |
| --- | ---: | ---: | ---: | ---: |
| caffeine | 1242 | 51 | 0.21 | 1.89 to 8.17 |
| copper guanidinate | 3060 | 124 | 0.21 | 2.94 to 3.86 |
| tagrisso | 3387 | 133 | 0.19 | 2.73 to 4.23 |
| taxol | 5489 | 223 | 0.20 | 3.01 |
| c60 | 4500 | 180 | 0.21 | 4.29 |
| Ti15O30 | 6270 | 285 | 0.32 | 4.09 |

The oxide has the most occupied orbitals and much the largest fitting set of
anything here, and it sits at the top of the table beside c60. Both quantities are
dimensions of the products the exchange is built from, and this system is large in
both while its orbital basis stays small. It is the clearest case in the file of
what the exchange wants: **many occupied orbitals and a large fitting set behind a
compact orbital basis.** Caffeine is the opposite of it in every one of those and
sits at the bottom.

### Two things about the input

The coordinates as given were in bohr, and the reader of an xyz file takes
angstrom. Read as angstrom the nearest titanium to oxygen distance is 3.43
angstrom and the cluster has no bonds at all: forty five atoms sitting apart from
one another, which would have run and returned a number rather than failing.
Divided by the bohr radius it is 1.82 angstrom and the coordination comes out as
above. **A geometry whose bonds are not checked can be computed successfully and
mean nothing.**

The two routes agree to 2.7e-09 here, which is the largest disagreement in this
file, against about 1e-11 elsewhere. Against a total energy of fifteen thousand
hartree that is two parts in ten to the thirteenth, so it is the arithmetic being
summed in different orders over twenty three iterations rather than a difference of
the approximation. It is worth recording that the absolute agreement of two routes
follows the size of the number they are computing.

## The four ways of building a Fock matrix

*Superseded. Tagrisso at Hartree-Fock and B3LYP, four builds, def2-svp and def2-svpd, with the B vectors timed apart from the builds:
`benchmarks/data/scf/2026-09-15_m4max_tagrisso.md`.*

### The direct way costs less than counting its passes suggests

The direct way holds no B vectors. It forms the three-center integrals again for
every batch of occupied orbitals and once more for the Coulomb matrix, which for
this molecule is several sweeps of them for every Fock matrix, where the way which
holds them sweeps them once for the whole calculation.

| | def2-svp | def2-svpd |
| --- | ---: | ---: |
| against the way which holds them | 2.57 slower | 2.06 slower |
| against the route VeloxChem had | **1.26 faster** | **2.37 faster** |
| against the four center build | **1.97 faster** | **7.84 faster** |

Those are the Hartree-Fock rows of
`benchmarks/data/scf/2026-09-15_m4max_tagrisso.md`, recomputed from it rather than
carried over from the table which used to stand above.

A factor of two for holding nothing, not the five or ten a count of the passes
would suggest. The reason is in the section on the setup: **the three-center
integrals are 2.6 per cent of forming the B vectors** and the contraction with the
metric is the rest, so repeating the integrals costs far less than repeating the
work around them.

That makes the direct way more than a fallback. It is faster than the route
VeloxChem had in both basis sets while holding a small fraction of the memory, on a
molecule which fits either way. The gap to the way which holds the B vectors also
narrows as the basis grows, 2.57 to 2.06, which is the direction that suits it:
the calculations which need it are the large ones.

These are the numbers as the driver stands, measured again after all of the work of
the sections below. Before the two corrections of the next section the direct way
took 100.54 and 283.05 seconds, and was 1.04 slower than the route VeloxChem had on
the smaller basis rather than 1.20 faster.

The work which came after those, and which was driven entirely by what a node with
many cores showed, was measured here again to see what it cost a machine with few.
**It cost nothing and gained a little**: the direct way fell 7.5 per cent at def2-svp
and 13.4 at def2-svpd, the way which holds the B vectors fell 4.8 and 1.8, and the
two builds which were not touched moved by under one per cent, which is the noise of
this machine. The exchange in particular now takes a triangle for every thread where
it took one matrix between them, which is more memory in exchange for work the
threads can divide -- a trade made for a hundred and twenty eight cores which turns
out to pay at sixteen as well, only by less.

## Where the direct mode spent its time

The direct way was two to four times slower than the way which holds the B
vectors, which sounded like the price of holding nothing. It was not. Most of the
gap was work being done more than once, and a profile of one Fock build said so
plainly.

Tagrisso in def2-svpd, direct mode, the whole SCF, the twenty four Fock builds
averaged. The driver was instrumented for the measurement and the instrumentation
thrown away afterwards.

| phase | before | after | share after |
| --- | ---: | ---: | ---: |
| build and zero the half transformed integrals | 0.056 | 0.049 | 0.7% |
| three-center integrals, first pass | 1.430 | 1.434 | 19.4% |
| **the half transform** | **3.709** | **1.703** | 23.1% |
| the closure onto the fitting coefficients | 0.066 | 0.068 | 0.9% |
| copies into and out of the stacked array | 0.201 | 0.262 | 3.5% |
| **the triangular solve** | 2.334 | 2.342 | **31.7%** |
| the exchange square | 0.709 | 0.696 | 9.4% |
| three-center integrals, second pass | 0.709 | 0.716 | 9.7% |
| the Coulomb matrix and its accumulate | 0.083 | 0.071 | 1.0% |
| rest | 0.022 | 0.045 | 0.6% |
| **total, seconds per Fock build** | **9.318** | **7.387** | |

**Every arithmetic phase already ran at 650 to 745 Gflop/s.** The copies and the
allocations, which is where one looks first, were 2.8 per cent between them. There
was no slow code to speed up.

### The half transform did three times the arithmetic it needed to

The sparsity pattern describes its blocks pair block by pair block and every
auxiliary group within one:

```cpp
const auto npatterns = static_cast<int>(blocks.size() * naux);
patterns[index].emplace(blocks[index / naux], aux_groups[index % naux], ...);
```

so the list runs pair block major, auxiliary group minor. The direct mode swept it
in runs sized to fit the integrals in memory, and **a contiguous run of that list
touches every auxiliary function**. Counted rather than assumed: 3387 auxiliary
functions carried blocks in each of the three sweeps, all of them, every time.

The dense path of the W build scatters one auxiliary function into a square of the
whole basis and hands it to a matrix product, and it did so whether or not the
sweep carried anything for that function. Three sweeps, three products for every
function, two thirds of them multiplying zeros — real arithmetic at a good rate,
for nothing.

Restricting the product to the rows a run touches would not have helped: a run
spanning every auxiliary group touches nearly every atom, so the square it needs is
nearly the whole square. **The fix is to sweep parts which own disjoint auxiliary
functions**, so a function belongs to one part and is transformed once. The
auxiliary basis has only four groups here, the largest of them 6.62 GB of a 10.675
GB whole, so the parts are cut by auxiliary atom rather than by group, each atom's
share measured from the blocks which carry it.

With the parts disjoint, a guard for a function no block of the call carries turns
from dead code into the thing which collects the saving. It is inert for the way
which holds the B vectors, which sweeps once and carries every function, and that
way is unchanged at 137 seconds.

### The integrals were formed once per batch of orbitals, and the batch was a constant

The batch of occupied orbitals came from a hardcoded four gigabytes rather than
from the memory budget the driver is given. Tagrisso in def2-svpd ran two batches
where the whole 133 orbitals fit in one, and so formed the integrals twice per Fock
build instead of once. Taking the batch from the budget removed the second pass.

The budget feeds two things which are live together — the half transformed
integrals held twice over, and the integrals of the part being swept — so each is
given half of it and the whole stays inside the number asked for.

### What is left is the floor

| | before | after |
| --- | ---: | ---: |
| tagrisso def2-svp | 100.54 | **87.72** |
| tagrisso def2-svpd | 283.05 | **237.38** |
| [Cu(PPh3)4]+ def2-svp | 10164.16 | **2466.47** |

The triangular solve is now the largest phase at 31.7 per cent, and it is the one
thing the direct way cannot avoid: **it solves against the factor on every
iteration, where the way which holds the B vectors solves once at the setup.** At
2.33 seconds a build over twenty four builds that is 56 seconds. The integrals are
the other thing paid every iteration rather than once, about 1.44 seconds a build
now that the orbitals fit in one batch, another 35 seconds. Together 91 of the 100
seconds by which the direct way trails the resident one, which is as close an
account as these phases allow. There is no third thing hiding in the gap.

### The molecule which fits no other way

Tetrakis(triphenylphosphine)copper(I), 137 atoms, closed shell d10 singlet at
charge +1, in def2-svp against def2-universal-jkfit: 1411 basis functions, 7256
auxiliary functions, 290 occupied orbitals. Both modes in one process.

| mode | time | iterations | energy |
| --- | ---: | ---: | ---: |
| RI-JK simd, direct | 2466.47 | 31 | -5755.2979742886 |
| full four-center | 2860.53 | 30 | -5755.2998281850 |

The B vectors would want **34.44 GB** and the conventional dense triangle **53.85
GB**, against 36 GB of machine. Neither the way which holds them nor the route
VeloxChem had can start this calculation. The direct way is the only one of the
four which runs, which is what it was written for.

It also wins. **1.16 times faster than making no approximation at all**, on a
molecule where the four center build is at its best: extended, mostly light atoms,
four phenyl-laden arms, which is the shape screening likes most. Before the two
corrections the same calculation took 10164 seconds and lost to the four center
build by 3.56, so the whole of that reversal is the duplicated work, not the
formulation.

The gain is larger here than on tagrisso — 4.12 against 1.19 — because the
duplication was paid more times over: twelve batches of orbitals rather than two,
and a larger integral set cut into more parts.

### A note on the budget

The table of the four ways pins the budget at 24 GB. Run instead with the default
this machine computes, 19.90 GB, the same def2-svpd calculation took 216 seconds
rather than 237. The two were measured in separate processes so the difference is
not established, but the direction is worth knowing: **more memory is not
automatically faster.** A larger budget makes one part where two would do, and the
one part plus the half transformed integrals held twice over comes to about 25 GB
on a 36 GB machine. Whether that is memory pressure or something else has not been
measured.

## The direct mode across the cores of a node

Everything above was measured on a laptop with sixteen cores. The work has to run
on a node with hundreds, and the first measurement there was disappointing in a way
the laptop could never have shown.

Tagrisso through the direct way on an AMD EPYC 9755, Turin, 128 cores at 2.7 GHz.
The driver is timed on its own rather than through an SCF -- the setup once, and one
Fock build, which is what every iteration repeats -- so that the diagonalisation and
the DIIS, which thread on their own, do not blur the curve.

### What the first curve said

| threads | def2-svp | speedup | def2-svpd | speedup |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 26.009 | 1.00 | 49.677 | 1.00 |
| 2 | 14.155 | 1.84 | 27.472 | 1.81 |
| 4 | 8.157 | 3.19 | 14.344 | 3.46 |
| 8 | 5.046 | 5.15 | 8.115 | 6.12 |
| 16 | 3.368 | 7.72 | 4.972 | 9.99 |
| 32 | 2.568 | 10.13 | 3.791 | 13.11 |
| 48 | 2.158 | **12.05** | 3.381 | **14.69** |

Twelve times on forty eight cores, a quarter of them doing anything. But the shape
of it says more than the number. Fit `T(n) = T(1) (s + (1 - s) / n)` and the whole
table follows from one figure: **a serial fraction of 6.3 per cent**, which predicts
7.74 seconds at four threads against 8.16 measured, 3.17 against 3.37, and 2.16
against 2.16. Nothing pathological, no threading storm. One serial section, and it caps the
speedup at `1/0.063`, about **sixteen times, however many cores are added**.

The larger basis scales better, which is the direction that suits the method: the
serial traffic grows as `naux nao nocc` while the arithmetic grows as `naux nao^2
nocc` and `naux^2 nao nocc`, so a bigger calculation dilutes the defect.

### Four things which did not widen with the threads

They were found by timing the phases of a build, behind VLX_RIJK_PROFILE, at one
thread and at forty eight. The phase whose time does not fall between the two is the
one to work on. None of the four had been visible on the laptop.

| what | why it was serial |
| --- | --- |
| the zeroing of the half transformed integrals | the matrix of one auxiliary function is smaller than the chunk the packed matrix divides its zeroing by, so its constructor zeroed it on the calling thread, several thousand times over |
| the copies into and out of the stacked array | plain loops, and the array was value initialised besides, for values every one of which is written before it is read |
| the closure onto the fitting coefficients | the vectorised inner loop was inside a serial outer one |
| the staging of the exchange update | a plain loop, gigabytes of it per build |

and a fifth, which was not serial at all but would not spread:

| what | why it would not spread |
| --- | --- |
| the rank k update of the exchange | the library divides an update over the blocks of the triangle it writes, and the triangle of a thousand functions holds about ten of them, which is nothing for forty eight cores. The depth is thousands, but the depth is the sum, which it cannot divide |

The last one is the interesting case, because it is not a defect in a loop of ours
but a mismatch between the shape of the work and the way a library divides it. The
answer is to divide it ourselves: **a triangle for every thread, a share of the
auxiliary functions each, the library left to run one update on one thread, and the
triangles summed at the end.** The functions are thousands and divide perfectly where
the triangle does not.

### What it came to

| at 48 threads, def2-svpd | build | speedup | serial |
| --- | ---: | ---: | ---: |
| as it was | 3.381 | 14.69 | 4.8% |
| the copies and the zeroing divided | 2.724 | 18.29 | 3.4% |
| a thread to each matrix it zeroes | 2.572 | 19.67 | 3.1% |
| the staging of the exchange divided | 1.857 | 27.0 | 1.65% |
| a triangle to each thread | **1.474** | **34.0** | **0.88%** |
| a thread to each Coulomb matrix, and the threads bound | **1.392** | **36.0** | **0.71%** |

**2.43 times faster on the same cores, and the ceiling from sixteen to a hundred
and forty.** Perfect scaling at forty eight threads would be 1.044 seconds, so what
is left is seventy five per cent of ideal, against thirty per cent at the start.

The last row is two changes at once, and they are separated in the sections below:
the Coulomb matrix was the one phase which took longer the more threads it was
given, and binding the threads to the cores is worth a third at 128 threads though
almost nothing at 48, which is why it appears here as a rounding rather than a
result.

### Where the time goes now

def2-svpd, one thread against forty eight, the fastest of three builds.

| phase | 1 thread | 48 threads | scaling |
| --- | ---: | ---: | ---: |
| integrals, second pass | 7.853 | 0.183 | 42.9 |
| integrals, first pass | 7.839 | 0.185 | 42.4 |
| the half transform | 12.726 | 0.317 | 40.1 |
| the triangular solve | 15.511 | 0.399 | 38.9 |
| the exchange | 4.753 | 0.142 | 33.5 |
| the zeroing | 0.490 | 0.046 | 10.7 |
| the copies | 0.415 | 0.051 | 8.1 |
| the closure | 0.099 | 0.014 | 7.1 |
| the Coulomb matrix | 0.340 | 0.066 | 5.2 |
| the rest | 0.115 | 0.071 | 1.6 |
| **whole build** | **50.142** | **1.474** | **34.0** |

The five which carry the work now run between 33 and 43 times, and are 1.23 seconds
of the 1.474. The five which lag are 0.248 seconds together, and most of that is not
a threading defect at all: **the copies move about 14.6 GB in 0.051 s, which is 286
GB/s, and that is the memory of the machine rather than its cores.** No number of
threads improves it; only not copying would. The rest, the part the phases do not
account for, is most likely the allocation and destruction of the sparse tensor of
integrals, which is gigabytes through mmap and munmap twice a build.

So about 0.2 seconds is what remains to be had from the code, which would be 1.26
seconds and a speedup near forty.

### Lifting the thread limit, and binding the threads

The limit was lifted by building OpenBLAS 0.3.30 again, in a directory of our own,
with `NUM_THREADS=256` and the same `COOPERLAKE` target. **The target was kept on
purpose.** OpenBLAS has a `ZEN` target, but it is an AVX2 one inherited from
Haswell, where COOPERLAKE is AVX-512, and Zen 5 carries AVX-512. The kernels of the
module were already reaching 70 to 93 per cent of what the cores can do, which is
not a library needing to be replaced. Only the thread count was wrong.

Two things about loading it are worth recording, as neither is obvious and the
first cost an afternoon.

**`LD_LIBRARY_PATH` did not work.** The library is named by an RPATH written into
the module at link time, and an RPATH is read before `LD_LIBRARY_PATH` is. The
`ldd` of the built module still named the one from the software tree, with the path
expanded as `.../lib/../lib64/`, which is the shape of a recorded search path rather
than of an environment one. `LD_PRELOAD` is read before either and is what worked.

**A process may hold more than one of them.** numpy carries its own, under a name of
its own, loaded by an absolute path, and a tool which asks the process for its
configuration will answer with whichever copy it meets first. The thread limit of a
library is a property of the file, not of the machine, so an answer from the wrong
file is worse than no answer -- it reports a limit which does not apply, or misses
one which does. The way to ask is to find which file answers for dgemm, by dladdr,
and read the configuration from that one.

### Binding the threads is worth a third

The node is two EPYC 9755, 128 cores each, 256 in all with no threading. Left
unbound, the threads wander between the sockets and read memory placed on the other
one. Bound, with `OMP_PROC_BIND=spread` and `OMP_PLACES=cores`, the first touch of
each buffer holds where it was made.

| threads | unbound | bound |
| ---: | ---: | ---: |
| 32 | 1.973 | 1.912 |
| 48 | 1.412 | 1.392 |
| 64 | 1.160 | 1.157 |
| 96 | 1.095 | 1.015 |
| 128 | 1.172 | **0.878** |
| 256 | 1.166 | 0.978 |

**Binding is worth 1.33 at 128 threads**, and it moves the best point from 96 to
128. Below 64 it is worth nothing at all, which is the tell: one socket holds 128
cores, so a run which fits inside one has nothing to wander across.

It also corrected a diagnosis. Unbound, the half transform grew with the threads --
0.256 at 64, 0.285 at 96, 0.377 at 128, 0.433 at 256 -- and the cause looked like
the memory it moves: it clears a square of the whole basis for every auxiliary
function, 27.6 GB a call at this size, which no number of cores can hurry. Bound,
the same phase reads 0.248, 0.200, 0.191, 0.284. **It was the placement, not the
bandwidth.** The zeroing is real but it was never the limit, and the work which
would have removed it would have bought little.

### The curve on two sockets

Tagrisso in def2-svpd, the direct way, bound, the fastest of three builds.

| threads | build | speedup | efficiency | serial |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 50.155 | 1.00 | 100% | |
| 8 | 6.727 | 7.46 | 93.2% | 1.04% |
| 16 | 3.467 | 14.47 | 90.4% | 0.71% |
| 32 | 1.912 | 26.23 | 82.0% | 0.71% |
| 48 | 1.392 | 36.03 | 75.1% | 0.71% |
| 64 | 1.157 | 43.35 | 67.7% | 0.76% |
| 96 | 1.015 | 49.41 | 51.5% | 0.99% |
| **128** | **0.878** | **57.12** | 44.6% | 0.98% |
| 256 | 0.978 | 51.28 | 20.0% | 1.57% |

**The best point is 128 threads, one socket's worth, and the second socket makes it
worse.** The serial fraction holds at about 0.7 per cent to 64 threads, rises to 1.0
at 128, and to 1.6 at 256, which is the interconnect appearing in the arithmetic.
At this size the calculation does not have enough work per unit of traffic to pay
for crossing between the sockets. A larger one might: the copper complex has twice
the basis and twice the orbitals, so four times the arithmetic against twice the
traffic.

**From where this started -- 3.381 seconds on 48 threads -- to 0.878 is 3.85 times**,
of which 2.4 is the code and 1.6 the threads and the library it may now use.

### What is left, and why it is not worth much

| phase | 1 thread | 128 threads | scaling |
| --- | ---: | ---: | ---: |
| integrals, first pass | 7.729 | 0.103 | 75.0 |
| integrals, second pass | 7.717 | 0.103 | 74.9 |
| the half transform | 13.018 | 0.191 | 68.2 |
| the exchange | 4.813 | 0.077 | 62.5 |
| the triangular solve | 15.570 | 0.257 | 60.6 |
| the zeroing | 0.366 | 0.009 | 40.7 |
| the copies | 0.412 | 0.029 | 14.2 |
| the closure | 0.096 | 0.007 | 13.7 |
| the Coulomb matrix | 0.356 | 0.059 | 6.0 |
| the rest | 0.077 | 0.043 | 1.8 |

The five which carry the work run between 60 and 75 times. The four which lag are
**0.138 seconds together, 15.7 per cent of the build**, and taking all four to the
rate of the others would give 0.75 seconds. **A ceiling of 1.17 for the hardest
work left**, which is where this stops being worth doing.

The Coulomb matrix is the largest of the four and the one which cannot simply be
divided harder: it holds a triangle for every thread and sums them at the end, so
both the memory it takes and the sum at the end grow with the cores. It no longer
grows in time, which was the defect worth fixing, but it will not fall much either.

That 1.17 is the ceiling for the Fock build alone, which is what every section above
this one measures. The section after it measures a whole calculation and finds most
of an iteration is no longer the Fock build at all, so the ceiling on the build is
not the ceiling on the calculation.

### The whole calculation, on the node

Everything above times a Fock build. This times what an SCF does with it: tagrisso
at restricted Hartree-Fock, the four builds one after another in a single process,
128 threads, bound, against the rebuilt OpenBLAS.

| basis | build | time | per iteration | against the exact one |
| --- | --- | ---: | ---: | ---: |
| def2-svp | four center | 33.21 | 1.581 | 1.00 |
| | RI-JK veloxchem | 44.92 | 1.953 | 0.74 |
| | RI-JK simd, in memory | 48.71 | 2.118 | 0.68 |
| | RI-JK simd, direct | **18.77** | **0.816** | **1.77** |
| def2-svpd | four center | 155.85 | 7.421 | 1.00 |
| | RI-JK veloxchem | 146.86 | 6.119 | 1.06 |
| | RI-JK simd, in memory | 104.54 | 4.356 | 1.49 |
| | RI-JK simd, direct | **51.80** | **2.158** | **3.01** |

**These were taken before the first exchange was fixed**, two sections below, and
every row of them carries one build of the four center integrals which the simd
rows no longer pay. On this node that build is about 7.4 seconds, so the two simd
rows should each fall by about six once measured again.

The three ways of the approximation agree to a ten thousand millionth of a hartree,
as they do everywhere else.

**The way which holds the B vectors is now the slower of the two, by two times.**
On the laptop it is 1.52 times the faster. On the node it is beaten by the direct
way at both basis sets and, in the smaller one, by the build which makes no
approximation at all. The two differ in one thing which matters here: the way which
holds them forms the W matrices in batches of sixty four auxiliary functions, so a
call has sixty four pieces of work to divide however many threads are waiting. The
direct way has no such batch.

That paragraph was written of this measurement and **neither half of it holds now**.
The batch was the first of four things wrong with the way which holds them, and with
all four mended it leads on both machines: 2.0 to 2.5 times on the laptop and 14 to
30 per cent on the node. The sections on the four builds on the node carry the
current figures and what each of the four was worth.

**The build which makes no approximation is also much better than the laptop said.**
It is within six per cent of the RI-JK route VeloxChem had, where on the laptop it
was 3.3 times behind. Screening divides over cores well; the old route does not.

### What an iteration spends outside the Fock matrix

The direct way at def2-svpd took **2.158 seconds an iteration** on the node, against
a Fock build of **0.878** for the same molecule, basis, mode, threads and binding.
That looked like 1.28 seconds an iteration of something else, and a first reading of
it blamed the diagonalisation.

Profiling a whole calculation said otherwise. The labels the driver keeps -- the
error vectors, the effective Fock matrix, the new orbitals, the new density -- come
to 0.19 seconds an iteration between them, and the diagonalisation is 0.06 of that.
**The rest was not an iteration cost at all.** It was one call, before the first
iteration, to the build of the four center integrals: on this machine 53.99 seconds
of a 207 second calculation, spread over twenty four iterations by the arithmetic
and made to look like overhead.

The reason is that the exchange is formed from the occupied orbitals, and at the
first iteration there are none, so **both** resolution of the identity routes fell
back to the build which needs none. Every RI-JK calculation paid one exact build
before it could start. The section below removes it from the simd route.

What is left after that, on the node, is about 0.97 seconds an iteration against
0.23 here -- real, and larger where the threads are many, which is the shape of
starting a hundred and twenty eight threads for work of a thousand rows. It has not
been profiled on the node and is not explained by anything measured here.

### The first exchange, from the density which has the orbitals in it

Any C whose product with its own transpose is the density gives that density's
exchange. The eigenvectors of the density, scaled by the roots of its eigenvalues,
are such a C, and the initial guess is a sum of atomic densities whose rank is the
occupied orbitals of the atoms: **218 for tagrisso, against 133 occupied orbitals of
the molecule and 683 or 1010 basis functions, and the same 218 in either basis.**
The eigenvalues below it are at the level of the arithmetic, 5e-16, so the rank is
sharp and no threshold has to be chosen carefully.

So the first exchange costs about 1.64 times an ordinary one, in place of a build of
the four center integrals. What that is worth, on this machine, with the energies
and the iteration counts unchanged in every case:

| calculation | before | after | |
| --- | ---: | ---: | ---: |
| caffeine def2-tzvp, in memory | 14.77 | 6.11 | 2.42 |
| caffeine def2-tzvp, direct | 20.10 | 11.20 | 1.79 |
| tagrisso def2-svpd, in memory | 137.25 | 84.36 | 1.63 |
| tagrisso def2-svpd, direct | 210.36 | 161.16 | 1.31 |

The savings match what the arithmetic says they should to within a second: direct
saved 44.3 seconds where one build of 54 replaced by 1.64 of 6 predicts 44, and the
way which holds the B vectors saved 50.6 against 48.6 predicted. **The way which
holds them gains more**, as its own builds are cheaper and the fixed cost was a
larger share of them.

It is worth most where the four center build is dearest beside the rest of the
calculation, which is the small molecule in the large basis: caffeine in def2-qzvpd
went from 311 seconds to 44.

**The conventional route still pays it.** It takes the orbitals as an object rather
than as a matrix of coefficients, so giving it the same treatment is more than the
one change made here. Every table in this file which sets the two side by side
therefore compares a route with this fixed against a route without it, and part of
each gain is this rather than the driver.

### The four builds on the node, and which mode to choose

Everything in the sections above was measured on a laptop with sixteen cores. The
same two molecules on the node, 128 cores, bound, with numpy given the machine.

Caffeine, Hartree-Fock, def2-universal-jkfit with 1242 auxiliary functions in every
row:

| basis | nao | four center | veloxchem | in memory | direct | best against exact |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 246 | 7.66 | 2.24 | **1.34** | 1.91 | 5.71 |
| def2-svpd | 366 | 12.29 | 4.47 | **2.25** | 3.01 | 5.46 |
| def2-tzvp | 494 | 35.48 | 8.25 | **3.21** | 4.01 | 11.05 |
| def2-tzvpd | 614 | 61.12 | 13.62 | **4.49** | 5.27 | 13.60 |
| def2-qzvp | 1098 | 487.83 | 59.38 | **12.41** | 14.42 | **39.30** |

Tagrisso, the same fitting set with 3387 functions:

| basis | nao | four center | veloxchem | in memory | direct | best against exact |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 683 | 29.51 | 40.85 | **12.75** | 12.81 | 2.31 |
| def2-svpd | 1010 | 142.00 | 136.08 | 30.76 | **29.68** | 4.78 |
| def2-tzvp | 1345 | 330.77 | 215.74 | **48.10** | 51.65 | 6.88 |

**The way which holds the B vectors is the one to use on a machine of this width.**
It leads five rows of caffeine by 14 to 30 per cent, leads tagrisso at def2-tzvp by
7, ties it at def2-svp and loses def2-svpd by 4. Every earlier version of this
section said the opposite, and the one which follows says why.

### What the way which holds the B vectors was paying, and for what

Against the tables this section has replaced twice, the way which holds them is
**2.8 to 3.8 times quicker on caffeine and 2.9 to 3.4 on tagrisso**, and the direct
way is unchanged to within a per cent -- which is what it should be, as nothing here
touched its path. The ratio between the two went from a flat 2.2 to 3.4 to between
0.70 and 1.04.

| in memory against direct | caffeine | tagrisso |
| --- | ---: | ---: |
| def2-svp | 0.70 | 1.00 |
| def2-svpd | 0.75 | 1.04 |
| def2-tzvp | 0.80 | 0.93 |
| def2-tzvpd | 0.85 | |
| def2-qzvp | 0.86 | |

Four things were in the way. None of them was where it looked, and one of the fixes
made the calculation slower until the next one landed.

**The range of the W matrices was fixed at sixty four.** The ratio above used to
hardly move across basis sets spanning 246 to 1345 functions, which is what a fixed
quantity predicts and not what anything growing with the calculation would. The
transformation divides the range it is given over the threads, so a call had sixty
four pieces of work however many threads were waiting -- sixty four being four times
the sixteen cores it was tuned on. Taken from the threads instead, tagrisso at
def2-tzvp fell from 150.41 to 103.51 seconds.

**Forming the B vectors cost forty four seconds of a hundred.** What was left after
the batch was not the Fock build at all: the builds of the two modes were 1.496 and
1.350 seconds, three and a half seconds apart over the whole calculation, while
everything outside them differed by 42.75. The mode which holds the B vectors forms
them before the first iteration and the direct mode does not, which is the whole of
what the two do differently outside a build. Timing the setup put 43.68 seconds on
forming them against 0.13 for the direct setup, and timing that put **42.51 of the
43.68 on the contraction of the metric with the integrals, and 0.27 on the integrals
themselves**. The front-loading of the integrals was doing exactly what it should
and costing nothing; the multiply afterwards was the price.

That multiply is the same operation the direct mode performs on every build, where
it costs 0.265 seconds and runs at 7.7 teraflops. Forming B carries all 1345 columns
where the direct solve carries only the 133 occupied ones, so six and a half times
the arithmetic is intrinsic. Running it at 318 gigaflops, a twenty fourth of the
rate, was not.

**The contraction issued one product per output group and threaded over the blocks.**
The comment above that loop already explained that the integrals of every atom basis
group are gathered into one buffer so that the depth of the product is the whole
auxiliary basis rather than one group's forty functions, "and one product replaces
eighty one". That was done for the depth and never for the rows: the output was
still written one group at a time, so the gathered buffer was read eighty one times
over and every product was fifteen rows tall. Counted: 0.29 GB gathered and 23.80 GB
read back from it on caffeine at def2-svp, the ratio 81.1 in every case measured.

Forming the whole output in one product and scattering it afterwards **made it
worse**, from 42.51 to 49.62 seconds. The rate rose by a third, to 427 gigaflops,
and the dense intermediate then formed every auxiliary row of a block including
those whose sparse block had been screened away -- 23.29 GB against 14.88 GB of B
vectors, 57 per cent more arithmetic, which swallowed the gain.

427 gigaflops on 128 cores was the tell. The contraction divided over the blocks of
atom pairs, and **a molecule has a dozen or two of those whatever its size**: the
distinct pairs of atom bases, fourteen for caffeine at every basis set and 81 for
tagrisso, against 128 threads. The columns of a block are its atom pairs, a hundred
or so, which do not divide usefully either. What is plentiful is the combinations of
angular momenta -- every one gathers and multiplies on its own and writes where no
other one does -- and a basis of triple zeta quality makes tens of thousands of
them: 31 116 for tagrisso at def2-tzvp where the blocks gave 81.

| tagrisso def2-tzvp, contraction | seconds | rate |
| --- | ---: | ---: |
| one product per output group, over blocks | 42.51 | 318 GF/s |
| one product per block, over blocks | 49.62 | 427 GF/s |
| one product per block, over the angular combinations | **3.33** | **6.35 TF/s** |

**Twelve and a half times, and the two changes only make sense together**: the
shape is what lets a product reach the matrix unit, and the threading is what puts
128 cores on it. Either alone reads as a regression. The rate now sits beside the
7.7 teraflops of the direct solve, so the 57 per cent of rows formed and dropped
costs about 1.2 seconds and is worth the shape that earns it.

**The B vectors were first touched by one thread.** With the setup down from 43.68
seconds to 4.48, what was left between the two modes was the half transformation,
which both call and which is 68 per cent of a build of the way which holds them.
Timing its phases put the answer beyond argument -- the two modes do the same work
there, and one phase carries all of the difference:

| a build of tagrisso def2-tzvp | in memory | direct |
| --- | ---: | ---: |
| entries | 0.007 | 0.001 |
| fill | 0.051 | 0.029 |
| scatter | **0.697** | **0.196** |
| product | 0.218 | 0.219 |
| rest, and what falls outside the phases | 0.306 | 0.011 |

**The product is identical to three digits**, which is the control: the same shapes,
the same count, the same library. The scatter moves the same 14.88 gigabytes in
both and takes 3.6 times as long, 21.3 gigabytes a second against 75.9.

The cause is not in the transformation. `CSparseTensor::zero` walked its blocks
serially, and since a large allocation is handed back untouched, that walk was the
first touch of all 14.88 GB -- placing every page on the memory of whichever socket
the master thread was on, for 128 threads to read afterwards, half of them across
the link. The setup's own report gave it away: 0.761 seconds to zero 14.88 GB is
19.6 gigabytes a second, which is one core and not a hundred and twenty eight.

**The direct mode never had the problem.** Its integrals tensor is allocated and not
zeroed; the pages are first touched by the integral kernel writing into them, which
is parallel, so they are spread over both sockets by the threads which made them.

Dividing the zeroing over the threads took it from 0.761 seconds to 0.054, and the
scatter from 21.3 to 59.6 gigabytes a second. Direct reaches 85.9, so about a third
of the gap remains: the zeroing places by block and the transformation reads by
auxiliary function, which is better than one socket and short of perfect.

**Seven calls a build where the direct mode makes one.** What a call of the
transformation costs beside its tasks -- the parallel region, and a square of the
basis allocated and zeroed for every thread -- is paid once for the call however
long its range is, and it was 0.3 seconds a build. The range is now the longest the
memory allows rather than the shortest the threads need, which is three calls at
def2-tzvp and one at every basis set of caffeine. It is worth most where the basis
is smallest, as those were the rows making twenty calls a build: caffeine at
def2-svp gains 30 per cent from it and def2-qzvp 14.

| tagrisso def2-tzvp | before | after |
| --- | ---: | ---: |
| zeroing the B vectors | 0.762 | **0.054** |
| scatter, a build | 0.697 | **0.268** |
| the half transformation, a build | 1.279 | 0.763 |
| the Fock build | 1.870 | **1.061** |
| the whole calculation | 61.29 | **48.10** |

**The Fock build of the way which holds the B vectors is now quicker than the direct
one**, 1.061 seconds against 1.320, having been 1.496 against 1.350 when this
started.

### What the approximation is worth depends on the machine

| caffeine, best route against the exact build | sixteen cores | 128 cores |
| --- | ---: | ---: |
| def2-svp | 7.9 | 5.7 |
| def2-svpd | 15.4 | 5.5 |
| def2-tzvp | 29.0 | 11.1 |
| def2-tzvpd | 46.0 | 13.6 |

**The approximation is worth one and a half to three times less on the larger
machine**, because the build which makes none divides over cores well and caffeine
is too small to keep a hundred and twenty eight of them busy: at def2-svp the whole
calculation is 1.34 seconds. The denominator falls faster than the numerator. The
gap has closed as the approximating routes were mended -- it was two to four times
when this column read 4.6 and 3.7 -- and what remains of it is the molecule.

**Caffeine at def2-svp is quicker on the laptop than on the node**: 1.26 seconds
against 1.34, with a ninth of the cores, both of them the way which holds the B
vectors. The exact build does go 1.6 times quicker on the node, having work enough
to divide, but the approximation has so little left to do that a hundred and twenty
eight cores cannot be given any of it. The molecule is the limit there, not the
machine.

It returns as the basis grows. At def2-qzvp the way which holds the B vectors is
thirty nine times the exact build on the node, the largest figure in this file, and
there the exact build has work enough to spread. **Its cost also steepens**: from
def2-tzvpd to def2-qzvp it grows as the 3.45 power of the basis where the earlier
steps gave 2.5, as screening loses its grip on diffuse and high angular momentum
functions.

**At def2-svp the conventional route loses to making no approximation at all** --
40.85 seconds against 29.51 for tagrisso, and at def2-svpd it ties the exact build
to within four per cent. On a node of this size VeloxChem's RI-JK is not worth using
below a thousand basis functions. The way which holds the B vectors used to lose
there too, at 43.93 seconds; it now takes 12.75 and is more than twice the exact
build.

### A caveat on the tables above

The four builds of a row are measured one after another in a single process, and the
direct way is measured last, after three others have taken and released gigabytes.
An earlier version of this file charged it 8 and 18 per cent for that, from a
comparison against runs made alone on another day -- 27.23 seconds at def2-svpd and
44.41 at def2-tzvp against 29.48 and 52.25 in the table. **That penalty is not there
now.** The direct way at tagrisso def2-tzvp measured 52.38, 52.61, 52.80 and 51.65
seconds with three other builds running before it, and 53.13 alone, so the order
costs nothing and the spread of the machine, under a per cent and a half, covers the
difference. The older figures predate several changes to the setup and should not be
read against the present table.

**Profiling costs little at the present call counts.** The same calculation with
VLX_RIJK_PROFILE set and without gave 48.71 and 48.10 seconds. It was worth 15 per
cent when the half transformation made seven calls a build, as the phases are summed
over the threads of every call; at three calls it is not worth correcting for.


### A note on measuring this at all

Two things made the laptop useless for this question, and both are worth knowing
before trusting a scaling number from one.

Accelerate exposes no way to set its thread count, so a run with one OpenMP thread
still has a fully threaded BLAS underneath it. Its one thread column is not one
thread, every speedup computed against it is wrong, and the phases which are mostly
BLAS -- the solve, the exchange -- appear not to scale at all, because they were
already parallel at the first point. **Five predictions were made from this side and
all five were wrong**, three of them about which phase was even the problem, and one
blaming the memory of the machine for what turned out to be where its threads were
standing. The node settled each of them in a single run, and the lesson is the dull
one: measure on the machine the answer is for.

The first build of a run is slower than the rest, by about eight per cent here,
which is first touch settling. Three builds and the fastest kept; one build
overstates.

## Binding the threads costs numpy a factor of thirty

This one is not about the driver at all, and it is the most broadly useful thing in
this file: **a VeloxChem calculation which pins its threads runs numpy's dense
algebra on a single core**, silently, and the cost grows with the basis.

### What it looks like

Tagrisso at def2-qzvp, 3099 basis functions, on 128 cores of the node, timing the
same calculation twice:

| per iteration | pinned | pinned, with numpy given the cores back |
| --- | ---: | ---: |
| the Fock build | 7.10 | 7.00 |
| the new orbitals | 3.83 | **1.05** |
| the error vectors | 2.34 | **0.14** |
| the new density | 0.32 | 0.06 |
| the effective Fock matrix | 0.25 | 0.12 |
| **whole iteration** | **14.83** | **9.41** |
| **whole calculation** | **355.95** | **225.72** |

**One and a half times, from two lines of setup.** The share of an iteration spent
in the Fock build goes from 48 back to 74 per cent.

### What is happening

A product of three thousand square, timed in the same process:

| | Gflop/s |
| --- | ---: |
| numpy alone, threads pinned | 3418 |
| numpy after veloxchem is imported, threads pinned | **105** |
| numpy after veloxchem is imported, threads not pinned | 3519 |
| numpy pinned, with the mask widened and its pool rebuilt | 3485 |

Neither pinning nor importing the code is enough on its own. Together they are, and
the reason is in what the libraries report of themselves:

```
libopenblas ......... 128 threads, threading_layer openmp     (the driver's)
libgomp ............. 128 threads
libscipy_openblas ...   1 thread,  threading_layer pthreads   (numpy's)
```

**numpy carries its own BLAS, and it is a pthreads build which sizes its pool from
the affinity mask the first time it is used.** Importing veloxchem loads libgomp,
which honours `OMP_PROC_BIND` by pinning the process to one core. numpy then builds
a pool of one thread and keeps it for the life of the process. The driver's own
library is untouched because it threads through OpenMP, which binds per region from
its own list of places rather than from the mask.

Three things which do not fix it, each tried:

| | |
| --- | --- |
| `OMP_WAIT_POLICY=passive`, `GOMP_SPINCOUNT=0` | 105 Gflop/s. It is not spin waiting |
| widening the mask after the first product | 105 Gflop/s. The pool already exists |
| forcing 64 threads without widening the mask | **16** Gflop/s. Sixty four threads on one core is worse than one |

### The fix

Widen the mask, then rebuild the pool, in that order, before numpy touches a
matrix:

```python
import os
import veloxchem as vlx                     # this is what pins the process

os.sched_setaffinity(0, range(os.cpu_count()))

from threadpoolctl import ThreadpoolController
ThreadpoolController().select(prefix='libscipy_openblas').limit(limits=64)
```

Only numpy's library is retargeted. Limiting every BLAS in the process would cap
the driver's own at sixty four, which it does not want.

### What is left after it, and what it says about the basis

The parts which are matrix products gain by the factor the measurement above
predicts -- the residual of the gradient by twenty, the transformation of the
orbitals by twenty one. **The eigen decomposition gains only 2.1**, from 0.322 to
0.152 seconds, because an eigensolver divides over cores far worse than a product
does. At 268 gigaflops in 0.152 seconds it is reaching 1.76 teraflops, which is as
much as it is going to give. That one is real work.

### It also undoes a conclusion this file nearly reached

Pinned, the part outside the Fock build appeared to grow as the cube of the basis
while the build grew as its square -- the fitting set does not grow with the orbital
set -- and the two were seen to cross at def2-qzvp: 48 per cent Fock build against
52 everything else. The reading was that the calculation around the driver would
overtake the driver as the basis grew.

The three basis sets measured again, with numpy given the machine:

| basis | nao | the iteration | the Fock build | the rest | the rest, as a share |
| --- | ---: | ---: | ---: | ---: | ---: |
| def2-svpd | 1010 | 1.134 | 0.90 | 0.23 | 21% |
| def2-tzvp | 1345 | 1.850 | 1.47 | 0.38 | 21% |
| def2-qzvp | 3099 | 9.405 | 7.00 | 2.41 | 26% |

and the whole calculations, 32.88 to 27.23, 52.74 to 44.41, 355.95 to 225.72. The
Fock build is unmoved in every row, which is what makes the three comparisons
clean.

**The share is flat at about a fifth, and there is no crossing in sight.** Taking
the two steps as they stand, the part outside the build grows with exponents of 1.68
and 2.21 against the build's 1.70 and 1.87. Near enough the same. **The cube was the
single core, not the arithmetic**: an eigen decomposition divides over the cores
better the larger it is, so its cost climbs far more slowly than its flop count, and
it was only the pinned measurement which made it look otherwise.

**The numbers of the sections above this one were taken pinned**, so the part of
them outside the Fock build is overstated -- by a fifth at def2-svpd and by a factor
of three at def2-qzvp.

## Density functional theory, where the integration is the other half

Everything above is Hartree-Fock. A hybrid functional adds the exchange correlation
integration to every iteration and asks for only a fifth of the exact exchange, so
it changes both sides of the balance at once. Caffeine and tagrisso at B3LYP against
def2-universal-jkfit, the grid at the level the functional asks for.

Caffeine on the laptop, sixteen cores, beside the same rows at Hartree-Fock:

| basis | nao | HF, in memory | B3LYP, in memory | B3LYP against the exact build |
| --- | ---: | ---: | ---: | ---: |
| def2-svp | 246 | 1.27 | 4.13 | 3.6 |
| def2-svpd | 366 | 2.69 | 8.86 | 6.3 |
| def2-tzvp | 494 | 4.81 | 13.73 | 13.1 |
| def2-tzvpd | 614 | 7.80 | 22.92 | 18.8 |

**B3LYP costs caffeine three and a half times what Hartree-Fock does here, and the
approximation is worth a third of what it was**: the advantage over the exact build
falls from 9.7 to 3.6 and from 54.2 to 18.8. The integration is a fixed addition
which lands on the quick routes and the slow one alike, so it hurts whichever is
quickest. It adds 2.7 seconds to a calculation of 1.27 and 8 seconds to one of 423.

Tagrisso on the node, 128 threads, all four builds:

| basis | nao | four center | veloxchem | in memory | direct | best against exact |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| def2-svp | 683 | 26.14 | 41.85 | **13.62** | 14.92 | 1.92 |
| def2-svpd | 1010 | 139.81 | 140.74 | 34.20 | **33.11** | 4.22 |
| def2-tzvp | 1345 | 302.22 | 217.49 | **45.29** | 48.07 | 6.67 |

**On tagrisso B3LYP costs what Hartree-Fock costs, and at def2-tzvp it costs less**
-- 45.29 against 48.10, a factor of 0.94, with 1.07 and 1.12 at the two smaller
bases. The integration adds 0.275 seconds an iteration and the fifth of the exchange
gives about as much back. The contrast with caffeine is the balance: the grid grows
with the atoms and the fitting set does not, so the integration is 1.5 times the
Fock build for caffeine and the build is 3.9 times the integration for tagrisso.

**The two ways of the driver are level under a hybrid** -- 0.91, 1.03 and 0.94 --
where at Hartree-Fock the way which holds the B vectors led by seven per cent at
def2-tzvp. What separates them is a smaller share of a larger iteration.

**The conventional route loses to making no approximation at all at two of the three
bases**, 41.85 against 26.14 and 140.74 against 139.81. A hybrid makes the exact
build cheaper, as only a fifth of the exchange is wanted, and the older route does
not gain from that. The crossing moves against it.

### What the integration is made of

The phases the integrator already timed, reported under VLX_XC_PROFILE. One call,
128 threads, def2-tzvp:

| phase | caffeine | | tagrisso | |
| --- | ---: | ---: | ---: | ---: |
| gtoeval | 0.013 | 26.5% | 0.063 | 22.9% |
| Generate density grid | 0.011 | 22.4% | 0.065 | 23.6% |
| Vxc matmul and symm. | 0.008 | 16.3% | 0.055 | 20.0% |
| Vxc matrix G | 0.004 | 8.2% | 0.025 | 9.1% |
| Density matrix slicing | 0.001 | 2.0% | 0.015 | 5.5% |
| XC functional eval. | 0.001 | 2.0% | 0.003 | 1.1% |
| Vxc dist. | 0.001 | 2.0% | 0.012 | 4.4% |
| rest | 0.010 | 20.4% | 0.036 | 13.1% |
| **total** | **0.049** | | **0.275** | |

**The functional itself is one per cent.** Whatever B3LYP costs, it is not the
formula: it is evaluating the basis functions on the grid points and multiplying by
them. Three phases -- the density on the grid, the Vxc product and the values
themselves -- are two thirds of the integration on both molecules, and the screening
and the distribution are nothing.

Measured on the laptop the three scale as the algorithm says they should: the two
products as the square of the basis, 2.09 and 2.15 against the 2.0 they must be, and
exactly linear in the points; the evaluation of the values as the first power of the
basis, 1.08. Nothing is being computed twice.

### One critical section was two thirds of it

The phases above are of the current build. Before it, on the same machine and the
same calculation, the distribution was **0.104 seconds of a 0.159 second call, sixty
five per cent**, where on the laptop it was 0.005 and one per cent.

| caffeine def2-tzvp, one call | laptop, 14 threads | node, 128 threads |
| --- | ---: | ---: |
| everything but the distribution | 0.386 | 0.028 |
| the distribution | 0.005 | **0.104** |

**The arithmetic scaled better than the cores did** -- 13.8 times on 9.1 times the
threads -- and the distribution went twenty times the other way. Every box added its
partial matrix into the shared Kohn-Sham matrix inside one `omp critical`, so the
threads queued. At fourteen it does not show; at a hundred and twenty eight it is
the calculation.

Giving each thread a matrix of its own and adding them at the end, as the exchange
of the RI-JK driver does with its triangles:

| caffeine def2-tzvp, node | before | after | |
| --- | ---: | ---: | ---: |
| Vxc dist. | 0.104 | **0.001** | 100 |
| rest | 0.027 | 0.010 | 2.7 |
| the whole call | 0.159 | **0.049** | 3.2 |
| an iteration | 0.327 | 0.212 | 1.54 |
| the calculation | 6.87 | **4.45** | 1.54 |

The `rest` line fell with it, which is what says most of that was tasks waiting on
the same lock rather than work.

**The energy moved by one unit in the tenth decimal**, -680.6196125485 against
-680.6196125486. A hundred and twenty eight partial sums are added where there was
one running total; at 1, 4 and 14 threads it is unchanged to every digit. It is the
order of the arithmetic.

**The memory is the square of the basis for every thread**, so it is bounded, and a
basis too large for the bound keeps the critical section rather than half of the
new way. Sixteen gigabytes carries four thousand functions at 128 threads and
twenty nine hundred at 256. It was two to begin with, which stopped at fourteen
hundred -- and tagrisso in def2-tzvp is 1345, which is closer than a bound should
ever be to the case it was measured on.

**Only the closed shell GGA path has this.** LDA, meta-GGA and both open shell paths
carry the same critical section, and will show the same two thirds on a machine of
this width.

### What is left

| def2-tzvp, node, one iteration | caffeine | tagrisso |
| --- | ---: | ---: |
| the Fock build | 0.073 | 1.060 |
| the integration | 0.049 | 0.275 |
| everything else | 0.090 | 0.691 |
| the Fock build against the integration | 1.49 | 3.85 |
| everything else, as a share | **42%** | **34%** |

**Neither of the two is the largest piece any more.** What is left is the
diagonalisation and the DIIS, a third to a half of an iteration in both molecules,
and nothing in this file has yet looked at it.

The integration also still divides its work by boxes of the grid, of which caffeine
at the level B3LYP asks for has 384 and tagrisso 1036 -- three and eight for every
thread at 128. The `rest` of the two reports, 20.4 and 13.1 per cent, is what that
granularity costs, and it eases as the molecule grows. On 256 threads caffeine would
have 1.5 boxes a thread, which cannot work.

## Three hundred and twenty atoms, where the integrals are three per cent

A cluster of paracetamol, 320 atoms and 640 occupied orbitals, at def2-svp against
def2-universal-jkfit: 3184 basis functions and 15888 auxiliary ones. The B vectors
would be 196.94 GB, which is most of the memory of the node, so this is the direct
way. Three iterations on 128 threads, bound.

| a build | seconds | |
| --- | ---: | ---: |
| solve | 69.5 | **44.7%** |
| transform | 66.6 | **42.8%** |
| exchange | 10.2 | 6.6% |
| integrals a + b | 4.6 | **3.0%** |
| allocate, copies, closure, coulomb, rest | 4.7 | 3.0% |
| **the whole build** | **155.6** | |

**The integrals are three per cent.** The direct way forms them again on every
build, twice, and at this size that costs 4.6 seconds of 155.6. Screening is what
does it: three hundred and twenty atoms in a cluster leave most pairs of them with
nothing to compute, while the two dense phases grow with the square of the basis
times the occupied orbitals. On tagrisso at def2-tzvp the same two sweeps were
19.7 per cent of a build. **The larger the system, the less of the direct way is
integrals** -- which is the opposite of what the name suggests, and it is worth
knowing before choosing what to make quicker.

**Eighty eight per cent of the build is two phases**, and neither is an integral.
The solve applies the Cholesky factor of the metric to the half transformed
integrals, 15888 squared by 3184 by 640, about 514 teraflops a build at 7.4
teraflops a second, which is what this machine gives for a triangular solve. The
transform is the half transformation, of which seven tenths is its product and a
fifth its scatter.

The setup is nothing: 0.63 seconds, five sixths of it the Cholesky factor of a
15888 by 15888 metric.

### A third cap which was chosen on a laptop

The exchange gathers the W matrices into triangles of the basis, one for each
thread, and `_syrk_triangles` bounded them together at two gigabytes. A triangle
here is 3184 squared by eight, 81.1 megabytes, so the bound allowed **twenty six of
them on a hundred and twenty eight threads**. A fifth of the build ran on a fifth of
the cores.

| | 2 GB | 16 GB |
| --- | ---: | ---: |
| taxol, 1099 functions | 128 of 128 | 128 |
| tagrisso def2-tzvp, 1345 | 128 of 128 | 128 |
| tagrisso def2-tzvpd, 1673 | 95 | 128 |
| paracetamol, 3184 | **26** | 128 |
| 4096 | 16 | 128 |

Raised to sixteen, the exchange of this build fell from 38.7 seconds to 10.2, and
three iterations from 880.2 to 731.9. Everything which should not have moved did
not -- the solve to a tenth of a second, the transform by one per cent, the
integrals not at all -- which is what makes the one number trustworthy.

**The rank k update itself went from 28.71 seconds to 7.26**, a factor of 3.95 where
twenty six threads to a hundred and twenty eight predicts 4.92. The missing fifth is
the price of adding a hundred and twenty eight triangles of 81 megabytes instead of
twenty six: 10.4 gigabytes of summation which was 2.1 before. It is worth paying
four times over, and it is a real cost rather than a rounding error.

**The energy moved by one unit in the tenth decimal**, -8191.8161913019 against
-8191.8161913018, from the same change of the order of the arithmetic.

**Three bounds were found too small on one day**, and all three were two gigabytes:
the copies of the Kohn-Sham matrix, the range of the W matrices, and the triangles
of the exchange. All three are the square of something by the threads, so all three
bind at about fourteen hundred and fifty basis functions on a hundred and twenty
eight threads, and all three were chosen on a machine with sixteen cores where that
is four thousand. The lesson is not about any of the three numbers.

## Dividing the resolution of the identity over the ranks

MPI in VeloxChem lives in the Python layer alone: the C++ takes a share of the work
and answers a partial quantity, and Python reduces. The SIMD RI-JK driver was
refused on more than one rank until now, with an assertion saying so. Three changes
lifted it, and each way of building is divided over a different index.

**The way which holds the B vectors is divided over the atoms of the auxiliary
basis.** Every term of both the Coulomb and the exchange is a sum over the auxiliary
basis, so a rank given a share of its atoms forms a share of every Fock matrix and
the existing reduction at the end of the build adds the shares.

**The direct way is divided over the orbitals, and over the parts it sweeps.** Its
exchange pass is a sum over the orbitals; the auxiliary basis is not an index it can
be divided over, because the triangular solve of that pass reaches across the whole
of it. Dividing it that way was tried first and gives a Fock matrix wrong in the
third figure -- a two way split of a matrix of order 1e3 was out by 7.2e+03, which is
the forward substitution missing the rows above its own. The Coulomb pass is a sum
over the parts and is divided that way, with the fitting between the two passes
gathered by the one communication a build makes, of one value per auxiliary function.

**The metric is inverted once on the master and broadcast.** Which fallback the
inversion takes is decided from the matrix itself -- the Cholesky factor, or the
inverted square root when a fitting basis is close to linear dependence -- so two
ranks could decide differently and build with metrics which are not the same. The
packed matrix learned to cross the ranks without a copy of itself for this, through
an array which writes into the matrix rather than into a buffer beside it.

### The energies do not move with the number of ranks

A water dimer and carbon monoxide, 8 atoms, def2-svp against def2-universal-jkfit,
against the conventional RI-JK driver of the same approximation.

| ranks | the B vectors held | the direct way |
| --- | ---: | ---: |
| 1 | -264.536932129285 | -264.536932129281 |
| 2 | -264.536932129285 | -264.536932129282 |
| 3 | -264.536932129285 | -264.536932129281 |
| 5 | -264.536932129285 | -264.536932129282 |

The conventional driver gives -264.536932129281, so the largest disagreement is
4.1e-12 and it does not grow with the ranks. B3LYP, which scales the exchange rather
than taking all of it, agrees to the same figure.

**Bit for bit was the wrong thing to ask for.** Calling the driver twice on one rank
with the same input already differs by 8.9e-16: the threaded reductions are not order
stable. Agreement to rounding is the most any division can be held to, and it is what
the shares meet.

### What the division was for: the memory divides

Tagrisso at def2-tzvp, 70 atoms, 1345 basis functions and 3387 auxiliary ones. The
memory is answered from the sparsity pattern before any integral is computed, which
is what the driver itself asks before choosing a way of building.

| ranks | the B vectors together | the largest rank | the smallest | spread |
| --- | ---: | ---: | ---: | ---: |
| 1 | 15.28 GB | 15.28 GB | 15.28 GB | 1.000 |
| 2 | 15.28 GB | 7.77 GB | 7.50 GB | 1.036 |
| 4 | 15.28 GB | 4.06 GB | 3.71 GB | 1.094 |
| 8 | 15.28 GB | 2.04 GB | 1.69 GB | 1.207 |

**This is the point of the exercise.** The atoms are dealt out sorted by their
distance from the centre of mass and by element, so the shares are not equal, but at
eight ranks the largest is a fifth above the smallest and nothing is lost: the sum
is the same 15.28 gigabytes at every count. A molecule whose B vectors do not fit on
one rank fits on enough of them, and the choice of way is made from the largest
share against the smallest budget rather than from the whole against one machine.

### The exchange does not divide, and gets worse

Caffeine at def2-tzvp, the way which holds the B vectors, **four threads on every
rank** so that what changes is the division and not the cores. One build, its phases.

| phase | 1 rank | 2 ranks | 4 ranks |
| --- | ---: | ---: | ---: |
| transform | 0.217 s | 0.116 s | 0.063 s |
| coulomb | 0.074 s | 0.067 s | 0.048 s |
| **exchange** | **0.034 s** | **0.058 s** | **0.092 s** |

**The transform divides, 3.4 times over four ranks. The exchange grows, 2.7 times
the wrong way.** The correctness tests could not see this, and did not: every rank
answers its share and the shares add to the right matrix. What they do not measure is
whether the work was divided at all.

The cause is that a rank sweeps the whole auxiliary basis whichever share of it it
holds. `_compute_in_memory` steps `first` from zero to `naux`, `compute_w_vectors`
builds its table of entries over the whole range asked for, and the functions this
rank carries nothing for are zeroed and then multiplied as zeros. So the rank
allocates the full range of W matrices, zeroes seven eighths of them at eight ranks,
and the rank k update of the exchange runs over all of them. **The transform divides
because only a function which carries something is transformed; the exchange does not
because a zero is multiplied like anything else.**

### A cap which is per rank on a machine which is not

Tagrisso at def2-svp on fourteen ranks of a laptop asked for **68 gigabytes on a
machine with 36**, and was killed.

| | |
| --- | ---: |
| basis functions | 683 |
| auxiliary functions | 3387 |
| occupied orbitals | 133 |
| one W matrix | 0.69 MB |
| the range `_w_batch_memory` allows | 23640 functions |
| the range after `min(naux, ...)` | 3387 |
| W matrices held by one rank | **2.29 GB** |
| by fourteen ranks | **32.09 GB** |

`_w_batch_memory` is sixteen gigabytes, and it is a constant of the process. Fourteen
processes on one node are therefore allowed two hundred and twenty four gigabytes
between them, and here the auxiliary basis cut the range down to 2.29 gigabytes
each -- which is 32.09 together, before the B vectors, the dense Fock matrices and
fourteen Python interpreters.

**The other two sixteen gigabyte bounds do not multiply, and it is worth being clear
why.** `_syrk_triangles` and the copies of the Kohn-Sham matrix are `min(nthreads,
cap / square)`: they are bounded by the threads as well as by the constant, and the
threads of a rank shrink as the ranks of a node grow, so their total over a node is
what one rank with all the cores would have taken. `_w_batch_memory` is bounded by
the auxiliary basis instead, and the auxiliary basis does not shrink when a node is
divided. A bound which is a constant is safe where something else already scales with
the share; this one had nothing else.

**`_get_ri_memory_budget` was divided by the ranks sharing a host and the caps were
not.** The budget counts the host names of the communicator and gives each rank its
share, which is right; the three sixteen gigabyte constants beneath it know nothing
about the communicator. Dividing the budget alone is worth very little when what the
driver allocates is not bounded by the budget.

### What a laptop can and cannot say

Fourteen cores, the total held constant, so that R ranks get 14/R threads each. This
measures what the division costs, not what it buys: no core is added, and a rank is
never cheaper than a thread doing the same work.

| caffeine | 1 rank | 2 ranks | 7 ranks | 14 ranks |
| --- | ---: | ---: | ---: | ---: |
| def2-svp, held | 1.02 s | 1.21 s | 2.15 s | 4.62 s |
| def2-svp, direct | 2.51 s | 3.68 s | 10.03 s | 20.03 s |
| def2-tzvp, held | 4.02 s | 4.65 s | 8.16 s | 17.16 s |
| def2-tzvp, direct | 8.67 s | 13.20 s | 38.69 s | 75.37 s |

Every row is slower, and that alone would say nothing -- one rank on fourteen threads
ought to beat fourteen ranks on one, since the OpenMP division of the same work has
no communication and no duplicated sweep. **What the rows do say is how much slower.**
Four and a half times on caffeine def2-tzvp held, eight and a half times direct, for
a division into fourteen. If the work divided cleanly these would be near one.

The direct way is the worse of the two here for a reason which is not a defect: each
rank sweeps the integrals for its own range of orbitals, so fourteen ranks make
fourteen sweeps where one rank made one, and the sweep is a fixed cost per rank. At
caffeine's size the integrals are a large share of a direct build. At three hundred
and twenty atoms they were three per cent, which is the size the division is for.

**Nothing here is a measurement of scaling.** Scaling is more cores, not the same
cores cut up, and it belongs on the node. What the laptop settles is that the answers
are right at every rank count, that the memory divides as it should, and that two of
the phases do not.

### Sweeping the share instead of the whole

The share of a rank is a set of auxiliary functions and not a range of them. The
dense index runs over the angular momenta of the whole molecule before it runs over
the atoms, so the functions of one atom are scattered through it, and a rank given
every fourth atom holds no interval of anything. The transformation could only be
asked for a range, `qfirst` to `qlast`, so the build asked for the range which covers
the share -- which is the whole auxiliary basis on every rank.

`compute_w_vectors` gained a form which takes the functions by name. It looks up
where a function belongs in the call, in a table of one entry per auxiliary function
of the molecule, instead of subtracting the first of a range. The build keeps the
dense indices its atoms give it and sweeps those.

**The functions a rank sweeps, caffeine at def2-tzvp, 1242 auxiliary functions:**

| ranks | swept by each rank | together |
| --- | ---: | ---: |
| 1 | 1242 | 1242 |
| 2 | 621, 621 | 1242 |
| 4 | 340, 340, 281, 281 | 1242 |

Before this they read 1242 on every rank at every count, and the sum was 1242 times
the ranks. The values of the B vectors always divided exactly -- 147 980 640 of them,
at one rank and at four -- which is why the memory divided and nothing else did.

**One build, one thread on every rank.** Read this as one rank's share against one
rank's resource: the cores grow with the ranks here, so it is not a comparison at a
fixed machine, and the section below is. It says whether the share of a rank shrank,
and nothing about what a division costs.

| phase | 1 rank | 2 ranks | 4 ranks | divides by |
| --- | ---: | ---: | ---: | ---: |
| transform | 0.625 s | 0.354 s | 0.174-0.208 s | **3.3** |
| exchange | 0.047 s | 0.026 s | 0.014-0.041 s | **2.5** |
| coulomb | 0.241 s | 0.202 s | 0.135-0.189 s | 1.5 |

The exchange divided by 2.5 where it had **grown** by 2.7, and the transform kept the
division it already had. Measured at four threads on every rank instead of one, the
Coulomb looks flat rather than 1.5; four ranks of four threads is sixteen threads on
fourteen cores and the phase is bandwidth bound, so that reading is of the laptop and
not of the code. **At one thread the Coulomb still divides by 1.5 where its values
divide by 4**, which is the code: the walk keeps its per block and per combination
overhead while only the values inside divide, and a rank holds the same 56 blocks
however few auxiliary atoms it was given. That one is not fixed.

**The whole calculation, caffeine def2-tzvp, fourteen cores divided R ways:**

| | 1 rank | 2 ranks | 7 ranks | 14 ranks |
| --- | ---: | ---: | ---: | ---: |
| held, before | 4.02 s | 4.65 s | 8.16 s | 17.16 s |
| held, after | 4.00 s | **4.31 s** | **6.07 s** | **8.01 s** |
| direct | 8.60 s | 13.08 s | 38.76 s | 71.32 s |

Fourteen ranks of one thread went from four and a quarter times the cost of one rank
of fourteen threads to twice it. It is still a cost, and it should be: the same work
on the same cores, divided by processes which do not share memory rather than by
threads which do. The direct way does not move, and should not -- it is divided over
the orbitals, and its sweep of the auxiliary basis is the whole of it by construction.

**And the memory a node is asked for stops growing with the ranks.** Tagrisso at
def2-svp, a budget of 28 GB divided by the ranks sharing the host:

| ranks | budget a rank | functions held | W matrices a rank | over the node |
| --- | ---: | ---: | ---: | ---: |
| 1 | 28.00 GB | 3387 | 2.29 GB | 2.29 GB |
| 2 | 14.00 GB | 1719 | 1.16 GB | 2.33 GB |
| 7 | 4.00 GB | 526 | 0.36 GB | 2.49 GB |
| 14 | 2.00 GB | 263 | **0.18 GB** | **2.49 GB** |

Two and a half gigabytes over the node at every rank count, against 32.09 at fourteen
ranks before. Two changes do it and both were needed: the range is now the share, and
the bound on it is the smaller of the constant and a quarter of the budget, which is
itself already divided by the ranks of the host.

### The test which would have caught it

Every correctness test passed throughout. They had to: a rank which sweeps a function
it holds nothing of adds a matrix of zeros, and the shares add to the right Fock
matrix whether the zeros were multiplied or skipped. **An energy cannot tell a
division which divides from one which only looks like it.**

`number_of_aux_functions` is on the driver for that reason, and the test asserts of a
partition that the shares sum to the auxiliary basis and that no share is the whole
of it. It is two numbers rather than a time, so it does not flicker with the machine,
and it fails on the code as it was written the first time.

## Two changes which measured nothing, and one which was not the plan

The Coulomb of the way which holds the B vectors divided by 1.5 where its values
divided by 4. A two term model of a probe which times `compute_y_vector` and
`compute_fock_matrix` on their own -- a cost following the blocks and a cost
following the values -- put the block bound term at 44 per cent of the phase at one
rank and 74 at eight, which pointed at the one thing in those loops which is paid per
block and does not depend on the auxiliary side:

```cpp
for (size_t k = 0; k < npairs_max; k++) {
    row = starts[a_atoms[k]*nmoms + lval_a] + ia + ma*strides[lval_a];
    col = starts[b_atoms[k]*nmoms + lval_b] + jb + mb*strides[lval_b];
    weights[k] = ... density.at(row, col) ...;
}
```

Two changes followed from that reading, and the redundancy was real: the blocks of
the B vectors which carry the same atom pairs differ only in their group on the
auxiliary side, and there are **four of them to a group for caffeine and 3.4 for
tagrisso**, so this gather was being done three or four times over.

| caffeine def2-tzvp | blocks | distinct lists of atom pairs | repeated |
| --- | ---: | ---: | ---: |
| 1 rank | 56 | 14 | 4.00 |
| 2 ranks | 56 | 14 | 4.00 |
| 4 ranks | 42 | 14 | 3.00 |
| 8 ranks | 28 | 14 | 2.00 |

**Both changes measured nothing.** Hoisting the places of the pairs out of the loops
over the angular components: 0.1490 seconds against 0.1491, which is noise. Gathering
the density once for a group of blocks instead of once for each of them, which is two
hundred lines and a restructured loop: **0.1602 against 0.1491, seven per cent
slower** -- the indirection through the group and the jumping between blocks for every
pair of components cost more than the gather it saved. Both were reverted.

**The cost was `CSparseTensor::_check_values`.**

```cpp
errors::assertMsgCritical(_values_state == valstat::allocated,
                          std::string("SparseTensor.") + label + std::string(": ..."));
```

The message is an argument, so it is built whether or not the check fails: two
concatenations of a string, and so two turns of the allocator, on a call which
otherwise reads one pointer out of a vector. `values(block, la, ia, lb, jb, lc, kc)`
is called once for every combination of basis functions of every block, which is
where the cost that followed the blocks was. `CSparseMatrix::_check_values` had been
fixed for exactly this, with a comment saying so; the tensor had been missed.

| caffeine def2-tzvp, one thread | before | after | |
| --- | ---: | ---: | ---: |
| 1 rank | 0.2397 s | **0.1494 s** | 1.60 |
| 2 ranks | 0.1997 s | **0.1092 s** | 1.83 |
| 4 ranks | 0.1310 s | **0.0696 s** | 1.88 |
| 8 ranks | 0.0715 s | **0.0375 s** | 1.91 |

The per value cost falls from 0.81 to 0.50 nanoseconds at one rank, and the phase now
divides by 4.0 over eight ranks where it divided by 3.4. It is not only the Coulomb:
the transformation calls the same accessor, and a build of caffeine def2-tzvp at one
thread moved from 0.625 to 0.534 seconds there and 0.241 to 0.150 in the Coulomb.

**The whole calculation, fourteen cores divided R ways:**

| caffeine def2-tzvp, held | 1 rank | 2 ranks | 7 ranks | 14 ranks |
| --- | ---: | ---: | ---: | ---: |
| before the sweep was the share | 4.02 s | 4.65 s | 8.16 s | 17.16 s |
| after it | 4.00 s | 4.31 s | 6.07 s | 8.01 s |
| after this | **3.73 s** | **3.76 s** | **5.14 s** | **8.08 s** |

**What to take from it.** The model was right that the cost followed the blocks and
wrong about which cost it was, and the two changes it recommended were built and
measured before that was known. Reading the loop for what is obviously expensive
found the density gather; what was actually expensive was a string built inside an
assertion which never fires. A cost model says where to look and not what to change,
and the only way to tell the difference is to measure each change on its own -- which
is what said that two hundred lines of grouping were seven per cent slower than the
code they replaced.

## The same cores divided into ranks, and the rank which finishes last

Every phase breakdown above is either one thread to a rank, where the cores grow with
the ranks, or four threads to a rank on fourteen cores, which oversubscribes. Neither
is what anybody runs. A hybrid calculation has `ranks x threads` equal to the cores
of the machine, and this is caffeine def2-tzvp at fourteen cores divided four ways.

**The number which matters is the slowest rank**, because the Fock matrices are
reduced and every rank waits for it. The range is given beside it.

| ranks x threads | transform | exchange | coulomb | the three together |
| --- | ---: | ---: | ---: | ---: |
| 1 x 14 | 0.075 | 0.022 | 0.024 | 0.121 |
| 2 x 7 | 0.082 (0.065) | 0.020 (0.012) | 0.027 (0.026) | 0.129 |
| 7 x 2 | 0.093 (0.054) | 0.017 (0.006) | 0.065 (0.030) | 0.175 |
| 14 x 1 | **0.152** (0.046) | 0.028 (0.004) | 0.075 (0.035) | 0.255 |

**The transform of the slowest rank is 0.152 seconds where the quickest is 0.046, a
spread of 3.3.** The phases of the quickest rank fall with the ranks exactly as the
sections above say they do. The build does not, because it ends when the last rank
ends, and taking a single rank's report -- which is what every measurement of this
session did -- shows the share and hides the wait.

### The spread is the machine, not the partition

The obvious reading of the table above is that the atoms are badly divided, and it is
wrong. The control is to give every rank provably identical work -- a square matrix
product of six hundred, twelve times over, timed as the best of seven -- and ask the
same machine what it does to that.

| ranks x 1 thread | identical work, spread | the transform, spread |
| --- | ---: | ---: |
| 2 | 1.00 | 1.01 |
| 4 | 1.96 | 2.16 |
| 8 | 3.01 | 2.35 |
| 10 | 3.14 | 3.94 |
| 14 | **4.93** | **4.00** |

**The same curve, and the transform sits below it.** This laptop is an Apple M4 Max:
`hw.perflevel0.logicalcpu` is 10 and `hw.perflevel1.logicalcpu` is 4, ten performance
cores and four efficiency cores. A rank which lands on an efficiency core takes two
to three times as long whatever it is doing, and the shape of the measured spread --
a tight cluster and a tail of two or three slow ranks at every count -- is that and
not a distribution of work.

**So the partition was not what the table showed.** Measured on what it does divide,
the round robin is uneven but nothing like three fold:

| predicted work of the heaviest rank over the lightest | by count | by work |
| --- | ---: | ---: |
| caffeine def2-tzvp, 4 ranks | 1.211 | 1.056 |
| caffeine def2-tzvp, 14 ranks | 1.268 | 1.233 |
| tagrisso def2-svp, 4 ranks | 1.094 | 1.021 |
| tagrisso def2-svp, 14 ranks | **1.293** | **1.078** |

`Molecule.partition_atoms` sorts the atoms by distance from the centre of mass and by
element and deals them round robin, which balances their count and so roughly the
count of auxiliary functions -- 93 and 95 at fourteen ranks -- but not the work. An
auxiliary function on an atom in the middle of a molecule survives screening against
far more atom pairs than one at the edge, so equal counts of functions are unequal
counts of values, and it is values which cost.

`Molecule.partition_atoms_by_weight` deals them heaviest first, each to the rank
carrying least so far, by weights the driver measures out of the sparsity pattern and
now answers through `aux_atom_weights`. It takes the ratio at fourteen ranks from
1.293 to 1.078 on tagrisso. **That is the whole of what it is worth**, and it is worth
having; it is not a cure for a three fold spread, because there was no three fold
spread to cure. Caffeine gains least because twenty four atoms over fourteen ranks is
one or two atoms a rank and there is nothing left to balance.

**What this does not say.** Fourteen cores of a laptop with two kinds of core is not a
node, and a molecule of twenty four atoms is not what any of this is for. What it
settles is the shape of the measurement: at a fixed machine, report the rank which
finishes last, never a single rank's profile -- and run the control before blaming the
code for what the machine did.

## The ranks of a node, on the node

Everything above this was measured on a laptop with two kinds of core, which is not
what a hybrid code is for. This is the node: two EPYC 9755 Turin, 256 cores, eight
NUMA domains of thirty two. A rank is meant to hold a domain, so the sweep is
1 x 256, 2 x 128, 4 x 64 and 8 x 32, the cores held constant and only the division
changing. Tagrisso at def2-tzvp against def2-universal-jkfit, 1345 basis functions
and 3387 auxiliary, B vectors of 15.28 GB, eighteen iterations.

| ranks x threads | Fock build | spread | control spread | orbitals | outside | whole |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 x 256 | 0.978 s | 1.00 | 1.00 | 0.199 s | 11.08 s | 31.34 s |
| 2 x 128 | 0.733 s | 1.07 | 1.04 | 0.212 s | 11.77 s | 26.69 s |
| 4 x 64 | 0.591 s | 1.06 | 2.77 | 0.178 s | 11.88 s | **23.53 s** |
| 8 x 32 | **0.565 s** | 1.09 | 1.76 | **0.157 s** | 13.44 s | 23.68 s |

**The Fock build is 1.73 times quicker on eight ranks than on one, with the same two
hundred and fifty six cores.** Ranks beat threads, which is the opposite of what a
division usually costs, and the reason is where the memory sits: one rank is one
allocation of fifteen gigabytes touched by two hundred and fifty six threads across
eight domains, and eight ranks is each rank's share sitting in the domain which
reads it. The whole calculation gains 1.33, best at four ranks with eight tied inside
the noise.

**Nothing here is load imbalance.** The slowest rank's build is within 1.09 of the
quickest, while the machine's own spread on identical work is 1.76 at eight ranks and
2.77 at four. The ranks differ less than the hardware does, so the partition by work
has nothing left to give and the next thing to fix is not it.

### What the division costs, which is not what was expected

| an iteration | 1 x 256 | 8 x 32 | |
| --- | ---: | ---: | --- |
| the Fock build | 0.978 s | **0.565 s** | 1.73 quicker |
| the orbitals, on the master | 0.199 s | **0.157 s** | 1.27 quicker |
| everything else | 0.153 s | **0.326 s** | **2.13 slower** |

**The eigen decomposition was the thing to worry about and is not.** It is behind a
rank guard, so it runs on the master's threads while every other rank waits, and
cutting the node into eight ought to have cost it eight fold. It got quicker: an
eigensolver divides over cores so poorly that thirty two threads beat two hundred and
fifty six, and the library numpy carries is built with a ceiling of sixty four in any
case, so the rows above it were never using what they asked for.

**What does cost is everything else**, which doubles: the extrapolation, the density,
the energy and the error vector, all of them master only or replicated, all of them
on threads which shrink as the ranks grow. It gives back about two fifths of what the
build gains, and at eight ranks the part outside the builds is 57 per cent of the
calculation. That is where the next work is, and it is not in the driver.

### Four ways the measurement lied before it told the truth

The numbers above are the fourth run of this sweep. The first three were wrong, and
each was wrong in a way worth writing down, because none of them looked wrong.

**numpy on one thread, every rank, every row.** The published fix in this file widens
the affinity mask to the whole machine before numpy builds its pool. Generalising it
to a communicator, the mask of a rank looked like the right thing to widen to instead
of the machine -- but a launcher binding a rank to a NUMA domain pins the process to
one core of it and lets OpenMP spread from its own list of places, so the mask reads
one core where the rank has thirty two. Sizing the pool from it put every rank's
linear algebra on one core. The control said so and was not read: 0.115 seconds on
one rank, on two, on four and on eight, when the same work on more cores cannot take
the same time.

**A driver on one core, in the row which had no launcher.** Told not to bind, the
launcher leaves the inherited mask, and libgomp builds its places from what is
available -- one core, two hundred and fifty six threads on it. The Fock build read
28.3 seconds against 0.73 on two ranks with the same cores, while numpy in that row
was the quickest of the four. A single rank is now run with no launcher at all, which
is how every pure OpenMP measurement in this file was taken.

**A printout which announced `4 ranks x 1 threads = 4 cores`** for a run on two
hundred and fifty six. The thread count was derived from the mask and presented as a
fact about the machine, when the driver and numpy take their threads by different
mechanisms and had different answers. Both are now read from their own sources -- the
driver's from the library, numpy's from threadpoolctl -- and the header says what was
asked for rather than asserting what was got.

**A label of 256 numpy threads on a library built with a ceiling of 64.** Asked for
more than it can give, it clamps and says nothing. The pool is now read back after it
is set, and where it cannot be read the line says `not verified` rather than
repeating the request.

The common thread is one mistake in four costumes: **printing a number that was
assumed in the place where a measured one belongs.** A control of identical work,
reported as a rate rather than as seconds, catches all four of them in the first row
of output -- about a hundred gigaflops a rank is one core whatever the launcher was
told. It was in the script from the start and it was reported in units which needed a
number nobody had.

### The direct way is not divided over the ranks of a node

The same sweep, the same molecule, the way which forms the integrals again on every
build instead of holding the B vectors.

| ranks x threads | Fock build | spread | outside | setup | whole |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 x 256 | **1.420 s** | 1.00 | 7.47 s | 1.58 s | **35.45 s** |
| 2 x 128 | 1.596 s | 1.22 | 8.03 s | 1.64 s | 38.97 s |
| 4 x 64 | 2.442 s | 1.23 | 8.31 s | 1.63 s | 55.23 s |
| 8 x 32 | 4.267 s | 1.22 | 8.99 s | 1.97 s | 90.73 s |

**Three times slower on eight ranks than on one, with the same cores.** Where the way
which holds the B vectors gained 1.73, this loses 3.0, and the whole calculation goes
from 35 seconds to 91.

**It is not the imbalance it was predicted to be.** The Coulomb pass of the direct way
is divided over the parts the auxiliary basis is swept in, and the parts are cut to
fit a memory budget rather than to fit the ranks, so a budget of hundreds of
gigabytes gives one part and one rank takes the whole of that pass while the others
take none. That is real -- seven of eight ranks take no part at eight -- and it is
not what costs: the spread of the slowest rank over the quickest is 1.22 at two
ranks, 1.23 at four and 1.22 at eight. Flat. An imbalance which one rank carries
alone would grow with the ranks, and this does not.

**What costs is that every rank sweeps every integral.** The exchange pass is divided
over the orbitals: a rank takes a range of them, forms the half transformed integrals
of that range, and for each batch of the range it sweeps the whole of the three
center integrals. The number of batches is set by the memory a rank may hold:

| | |
| --- | ---: |
| one orbital of the half transform, twice over | 2 x 3387 x 1345 x 8 = 72.9 MB |
| orbitals a batch may hold, at a budget of hundreds of gigabytes | about 5000 |
| occupied orbitals of the molecule | 133 |
| batches a rank makes, at any rank count | **one** |

One rank makes one sweep of the integrals. Eight ranks make eight, one each, on a
thirty second of the cores apiece. **The integral work multiplies by the number of
ranks**, and it does so precisely when the memory is plentiful, which is the case
this machine is.

The plan for this division said the opposite -- that the sweeps a rank makes fall as
its orbitals do, so the integral work is conserved -- and that is true only where
`nbatch` is smaller than the occupied orbitals, so that there are several batches to
divide. It was written for a molecule of three hundred and twenty atoms where the
integrals are three per cent of a build and the batches are several; it was measured
on one of seventy where they are a fifth and the batch is one. Integrals at 19.7 per
cent of a build, multiplied by eight, is 158 per cent of the original build in
integrals alone. Measured: 0.28 seconds of integrals a build at one rank, eight
sweeps on an eighth of the cores is 2.24, plus 1.14 for the phases which do divide,
against 4.27 measured. The mechanism accounts for it.

### What to run, then

| tagrisso def2-tzvp, 256 cores | best | |
| --- | --- | ---: |
| the B vectors held | **8 x 32** | 23.68 s |
| " | 4 x 64 | 23.53 s |
| the direct way | **1 x 256** | 35.45 s |

**The way which holds the B vectors wants the ranks; the direct way wants one rank
and all of the threads.** Within a node the direct way should not be divided at all,
and across nodes each node will duplicate the sweep, which is the price of not
sending several gigabytes of half transformed integrals between them.

Two ways out, neither taken yet. The parts of the exchange pass could be divided over
the ranks as well as the orbitals, with the half transformed integrals summed between
them -- 4.85 gigabytes a batch for this molecule, against eight sweeps of the
integrals, and which of the two is dearer is a question with a number rather than an
opinion. Or the orbitals could be divided only as far as the batches go, leaving the
ranks beyond that to divide the parts instead, which costs nothing to decide and
needs the same sum.

### Why the duplication is not worth removing on one node

The obvious fix for the direct way is to divide the parts of its exchange pass rather
than the orbitals, summing the half transformed integrals between the ranks so that
each integral is formed once. It was planned and then not built, because of an
arithmetic which applies to any division at a fixed machine.

Take `C` cores split into `R` ranks of `C/R`. A phase whose work divides by `R` runs
on `C/R` cores, so its wall time is `(W/R)/(C/R) = W/C` -- **the same as one rank**. A
phase which is duplicated runs its whole work on `C/R` cores and takes `R` times as
long. **No phase can be quicker than it is on one rank**, so no division of the direct
way can beat one rank on one node; removing the duplication would take the eight rank
build from 4.27 seconds to about 1.42 plus the communication, and 1.42 is the one rank
build.

**The way which holds the B vectors gained 1.73 from the ranks, and that is not a
counter-example.** Its gain is not from dividing work but from where the memory sits:
fifteen gigabytes of B vectors touched by two hundred and fifty six threads over eight
NUMA domains, against each rank's share sitting in the domain which reads it. The
direct way allocates its half transformed integrals afresh for every batch and first
touches them in parallel, so they are already placed, and there is no such gain to
collect. **A division at a fixed machine wins only where the threads were not using
the machine properly.**

**Across nodes it is a different question and the answer is not yet known.** With one
rank to a node the integral sweep is formed on every node while everything else
divides, which is an Amdahl term equal to the integral fraction of a build: a fifth
for tagrisso at def2-tzvp, which caps the speedup at five however many nodes are
given, and three per cent for the cluster of three hundred and twenty atoms, which
caps it at thirty three. That ceiling is real and the parts division would lift it.
Whether it binds before the nodes run out is a measurement rather than an argument.

### Two nodes, and what the division was built for

Sixteen ranks of thirty two threads over two nodes, five hundred and twelve cores,
against the same molecule on one node. The launcher is srun with PMIx; the site's
mpirun has no transport between nodes.

| | nodes | cores | Fock build | spread | outside | setup | whole |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 x 256 | 1 | 256 | 0.978 s | 1.00 | 11.08 s | 4.73 s | 31.34 s |
| 8 x 32 | 1 | 256 | 0.565 s | 1.06 | 13.44 s | 4.87 s | 23.68 s |
| 16 x 32 | 2 | 512 | **0.305 s** | 1.16 | 10.87 s | 3.79 s | **16.58 s** |

**The build is 1.85 times quicker on two nodes than on one, against a doubling of the
cores** -- ninety three per cent of what the cores could give, across an interconnect,
for a phase which sends one Fock matrix and one metric between the ranks and nothing
else. The B vectors of a rank are 0.95 gigabytes where one rank held 15.28, which is
the whole point: the molecule which does not fit fits on enough of them.

**The part outside the builds fell as well**, 13.44 seconds to 10.87, because the
setup divides too -- the B vectors are formed by sixteen ranks rather than eight. It
is still 66 per cent of the calculation, and it is still where the next work is.

| an iteration | 1 x 256, one node | 16 x 32, two nodes |
| --- | ---: | ---: |
| the Fock build | 0.978 s | **0.305 s** |
| the orbitals, on the master | 0.199 s | 0.134 s |
| everything else | 0.153 s | 0.259 s |

**Two cautions about these rows.** The single node rows were taken where the module's
BLAS has a high thread ceiling and the two node rows where it has one of forty eight,
so the rows asking for sixty four threads a rank and more are not on the same library;
the rows at thirty two threads are unaffected, which is the comparison above. And the
control of equal work now reads 110 to 1277 gigaflops a rank, a spread of 11.6, which
is the measurement and not the machine: numpy's mask is widened to the whole node so
that its pool can be built at all, and sixteen pools of thirty two threads then roam
two nodes unplaced. The driver's own spread is 1.16 and is placed by OpenMP, which is
why the build is trustworthy where the control beside it is not.

### The direct way across two nodes, where the duplication finally shows its price

The rows above which asked for more than forty eight threads a rank were on the
module's BLAS and warned; with the rebuilt OpenBLAS preloaded they run clean.

| direct, tagrisso def2-tzvp | nodes | cores | build, slowest | quickest | spread | whole |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 x 256 | 1 | 256 | 1.420 s | 1.420 s | 1.00 | 35.45 s |
| 2 x 256 | 2 | 512 | **1.004 s** | 0.772 s | **1.30** | 26.81 s |

**The model of the last section predicted 0.85 seconds and the quickest rank reached
0.772.** Taking `build(N) = T_int + T_other/N` with the integral sweep at 0.28 seconds
and everything else at 1.14, two nodes should give 0.85; the rank which is not carrying
the Coulomb pass gives 0.772, near enough that the model is the right one.

**The 0.232 seconds between the two ranks is the single part.** The Coulomb pass is
divided over the parts the auxiliary basis is swept in, the parts are cut to fit a
memory budget, and a budget of hundreds of gigabytes gives one part -- so one rank
sweeps the integrals a second time for the Coulomb matrix and the other waits. On one
node this was invisible, because the duplicated exchange sweep was eight times larger
than it and swamped it. On two ranks it is a fifth of the build.

**So the direct way has two defects and they are now separable**, which they were not
before:

| | what it costs | |
| --- | ---: | --- |
| the integral sweep formed on every node | 0.28 s, fixed | caps the speedup at about five nodes' worth |
| the Coulomb pass on one rank | 0.232 s, on one rank | a part count which knows how many ranks there are |

The second is a small change -- the parts are cut by memory alone and could be cut by
memory or by the rank count, whichever gives more -- and it is worth making before the
first, which is a restructure. Neither is worth anything on a single node, where no
division of this way can beat one rank.

### Cutting the two passes of the direct way apart

The Coulomb pass of the direct way is divided over the parts the auxiliary basis is
swept in, and with one part it lands on one rank. Cutting more parts balances it --
and the exchange pass sweeps every part on every rank, so cutting more parts also
gives every rank another call of the transformation, which costs a square of the basis
allocated and zeroed for every thread whatever the part holds. At 1345 functions on
256 threads that square is 3.7 gigabytes a call.

Tagrisso def2-tzvp, two nodes, one rank each, 512 cores:

| | parts | build, slowest | quickest | spread | |
| --- | --- | ---: | ---: | ---: | --- |
| cut by memory alone | 1 | 1.004 s | 0.772 s | 1.30 | one rank carries the Coulomb pass |
| cut for the ranks, both passes | 3 and 3 | 1.034 s | 1.017 s | 1.02 | balanced, and two extra sweeps |
| **cut apart** | **3 and 1** | **0.947 s** | 0.915 s | 1.03 | balanced, one sweep |

**Balancing the Coulomb pass by cutting the sweep cost more than it saved**, which is
why the second row is the slowest of the three despite being the best balanced. The
parts now come in two lists: the sweep is cut by memory, as it was, and the Coulomb
pass is cut again and finer when more parts are asked for than the memory gave.

The average work of a rank is 0.931 seconds against the 0.888 the model predicts for a
perfect division, and the difference is the Coulomb pass paying for three patterns and
three contractions where it paid for one. That is the right trade: its parts are swept
once between the ranks, not once each.

**A caution on the whole-calculation column.** Across four runs of this configuration
the part outside the builds measured 7.15, 8.05, 8.22 and 8.76 seconds, a fifth of
spread, while the build -- measured on every rank and averaged over seventeen
iterations -- moved by a few per cent. The build is the number to read here; the totals
are too noisy at this size to rank three configurations by.

### Naming every BLAS in the process before timing anything

Two of today's failures were a library quietly doing something other than what was
asked: numpy's pool built on one core, and the module's OpenBLAS given 256 threads
when it was compiled for 48, which warned once per thread and then died in its
allocator. Both are visible in one place, and it now prints before the run:

```
# blas libopenblas_cooperlakep-r0.3.30.so: 256 threads, openblas
# blas libgomp.so.1.0.0: 256 threads, openmp
# blas libscipy_openblas64_-f48b354e.so: 64 threads, openblas  <-- BELOW THE THREADS THIS RANK WAS GIVEN
```

**A process holds several, they are given their threads by different mechanisms, and
they do not agree.** The driver's takes them from OpenMP and is untouched by a launcher
pinning the process; numpy's sizes a pool from the affinity mask and is ruined by it.
Reporting one of them as though it described the process is how a run on two hundred
and fifty six cores came to be labelled as four.

## The math library, and what the Eigen fallback costs

The dense linear algebra has two implementations. `Makefile.setup` selects a hardware
math library per platform and defines `VLX_USE_MATHLIB`; where no library is named the
define is absent and an implementation in Eigen is compiled in its place, so a machine
without one needs no further change. Three files carry the branch -- the packed linear
algebra, which inverts the metric, and the two resolution of the identity drivers,
which hold the half transformation, the triangular solve and the rank k update of the
exchange.

The two built on this laptop, an M4 Max of fourteen cores, Accelerate against Eigen.
The Fock build of an iteration:

| molecule, basis, way | Accelerate | Eigen | Eigen is |
| --- | ---: | ---: | ---: |
| caffeine def2-svp, held | 0.030 s | 0.046 s | 1.53 |
| caffeine def2-svp, direct | 0.131 s | 0.499 s | **3.81** |
| caffeine def2-tzvp, held | 0.120 s | 0.177 s | 1.48 |
| caffeine def2-tzvp, direct | 0.430 s | 1.208 s | **2.81** |
| tagrisso def2-svp, held | 0.845 s | 1.491 s | 1.76 |
| tagrisso def2-svp, direct | 3.155 s | 23.761 s | **7.53** |
| tagrisso def2-tzvp, held | 3.937 s | 6.244 s | 1.59 |
| tagrisso def2-tzvp, direct | 9.801 s | 51.578 s | **5.26** |

The energies agree to the twelfth decimal in every row, which is what makes the
comparison a comparison.

**The two ways are not affected alike, and the two figures do not mean the same
thing.** Both builds were run at fourteen threads, so neither library was held back;
what differs is where the threads come from.

| operation | where it is called | Eigen | the math library |
| --- | --- | --- | --- |
| the half transformation | inside an OpenMP region | one thread a call, fourteen calls at once | the same |
| the rank k update of the exchange | inside an OpenMP region | one thread a call | the same |
| **the triangular solve** | **at the top level** | **one thread** | **fourteen** |
| the Cholesky of the metric | at the top level, in the setup | one thread | fourteen |

**The 1.5 to 1.8 of the way which holds the B vectors is a comparison of kernels.**
Its build is a half transformation and a rank k update, both called from inside a
parallel region with one small product to a thread, so neither library threads and what
is measured is Eigen's arithmetic against the matrix unit of the machine.

**The 2.8 to 7.5 of the direct way is mostly a comparison of thread counts.** It solves
the Cholesky factor of the metric against the half transformed integrals on every
build, at the top level, and Eigen parallelises general matrix products and little
else -- a triangular solve and a Cholesky are not among them. So that phase ran on one
core of fourteen where the library used all of them. It is a fair account of what each
build gives, and it is not a fair account of the two kernels.

**The setup divides the same way.** Forming the B vectors of tagrisso def2-tzvp took
31.73 seconds against 58.01, and at def2-svp 8.11 against 15.21 -- 1.8 in both, which
is the matrix product figure, as it should be.

**The whole calculation of tagrisso def2-tzvp by the direct way was 197.53 seconds
against 1020.19**, a factor of 5.2, from a build flag.

**What to take from it.** The Eigen path is a fallback for correctness and not for
speed, and the difference is large enough to change which way of building is the
quicker one: with Accelerate the direct way costs 2.5 times the held way at tagrisso
def2-tzvp, and with Eigen it costs 8.3. A machine built without a math library named
would reach conclusions about this driver which do not hold on a machine with one.

## c60 in def2-tzvp, which was not the driver

Open since 12 September: c60 at restricted Hartree-Fock in def2-tzvp did not converge
in fifty iterations through the SIMD RI-JK path. The first ten iterations were healthy
and monotonic, reaching a gradient of 1.1e-4, so it was a late stall and not an
unstable start. It was left open because each attempt was believed to cost 1.9 hours
on this machine, and the comparison against the four-center build another five or six.

It took two eigendecompositions and no SCF at all.

**The fitting metric is not the difference, and cannot be.** `def2-universal-jkfit`
does not depend on the orbital basis, so c60 has the same 4500 auxiliary functions and
the same metric at def2-svp, which converges, as at def2-tzvp, which does not:

| molecule | basis | naux | smallest | condition |
| --- | --- | ---: | ---: | ---: |
| caffeine | def2-tzvp | 1242 | 1.812e-06 | 9.877e+08 |
| tagrisso | def2-tzvp | 3387 | 4.259e-07 | 7.139e+09 |
| c60 | def2-svp | 4500 | 3.805e-08 | 1.322e+11 |
| c60 | def2-tzvp | 4500 | 3.805e-08 | 1.322e+11 |

c60's metric is the worst conditioned of the four by twenty times, and it is the same
1.322e+11 in the run which converges. A number shared by both sides of a comparison
explains neither.

**The orbital basis is the difference.** The twenty smallest eigenvalues of the overlap:

```
c60 def2-svp   3.202e-05  3.241e-05  3.279e-05  5.825e-05 ...
c60 def2-tzvp  9.043e-07  9.088e-07  9.151e-07 | 1.441e-06  1.463e-06  1.475e-06 ...
```

The default `ovl_thresh` is 1e-6, and it falls inside a cluster: **three directions are
dropped at 9.0e-07 and five are kept at 1.44e-06**, within a factor of 1.6 of the ones
just discarded. Near-null directions retained and allowed into the orbital rotations
are the textbook late stall. At 1e-5 the whole cluster goes -- 37 of 1860 dropped, the
smallest kept 1.564e-05, a clean gap.

| | smallest overlap eigenvalue | below 1e-6 |
| --- | ---: | ---: |
| caffeine def2-tzvp | 4.518e-05 | 0 |
| tagrisso def2-tzvp | 4.321e-06 | 0 |
| c60 def2-svp | 3.202e-05 | 0 |
| **c60 def2-tzvp** | **9.043e-07** | **3** |

**With `ovl_thresh = 1e-5` it converges in twenty five iterations** to
-2272.3448501388, monotonic throughout, passing through the region it used to stall in
at iterations ten and eleven without pausing. 1330 seconds on the M4 Max by the direct
way.

**A late stall is worth checking against the overlap spectrum before the Fock driver is
suspected.** One eigendecomposition, four seconds, against 1.9 hours an attempt -- and
the driver was never in question: c60 at def2-svp agrees between the two ways to
9.9e-10 and caffeine at def2-tzvp to 1.4e-11.

**What this does not settle.** The direct way solves with the Cholesky factor of the
metric and has no other option, so that run used it throughout. Whether the eigenvalue
route -- which the conventional RI-JK driver uses, and which `ri_metric_route` now
selects for the way which holds the B vectors -- also cures it at the default
`ovl_thresh` is untested. The two are independent and both may be real.

### And it was not a stall either

The section above is right that the driver is not at fault and wrong about what is.
It was written from one controlled comparison -- `ovl_thresh` varied within the direct
way on a laptop -- against a failure recorded on 12 September with older code. The
control was never run. Run now, with the same binary, and on a node:

| machine | way | ovl_thresh | verdict | iterations | energy |
| --- | --- | ---: | --- | ---: | ---: |
| laptop | direct | 1e-6 | **not converged** | 50 | -2272.3450690097 |
| node, 8 x 32 | direct | 1e-6 | converged | 25 | -2272.3450690097 |
| node, 8 x 32 | held | 1e-6 | converged | 25 | -2272.3450690105 |
| laptop | direct | 1e-5 | converged | 25 | -2272.3448501388 |
| node, 8 x 32 | held | 1e-5 | converged | 25 | -2272.3448501395 |

**The run which does not converge reaches the same energy as the runs which do**, to
ten decimals. Its last ten iterations:

```
 41  -2272.345069009752  grad 1e-08      46  -2272.345069009743  grad 1e-08
 42  -2272.345069009749  grad 1e-08      47  -2272.345069009752  grad 2e-08
 43  -2272.345069009758  grad 2e-08      48  -2272.345069009749  grad 1e-08
 44  -2272.345069009758  grad 1e-08      49  -2272.345069009734  grad 2e-08
 45  -2272.345069009749  grad 1e-08      50  -2272.345069009745  grad 1e-08
```

The energy is settled to 2.4e-11 and the gradient oscillates between 1e-8 and 2e-8.
**It is not a stall: it is a noise floor which sits exactly on `conv_thresh`.** At
1e-5 the floor lands just under the threshold and the test passes at iteration 25; at
1e-6 it lands just over and fifty iterations report failure. The node, with OpenBLAS
instead of Accelerate and eight ranks of thirty two threads instead of fourteen, lands
just under at 1e-6 and converges in twenty five.

**So c60 at def2-tzvp is a calculation whose achievable gradient is `conv_thresh`**,
and whether it is called converged is settled by the arithmetic: the math library, the
thread count, the rank count, the linear dependence threshold. None of them change the
answer. `conv_thresh = 1e-7` is met at iteration 21 in every configuration.

**Raising `ovl_thresh` is not a free fix.** 1e-5 drops 37 directions where 1e-6 drops
three, and the energies differ by 2.19e-04 Eh accordingly -- it is a smaller
variational space, not a better converged one. What it also does is lower the gradient
noise slightly, which is why it crossed the threshold; that was the effect measured
and mistaken for the cause.

**And the metric route changes none of it.** On the node, cholesky against
eigenvalues: converged in twenty five iterations either way, in both ways of building,
differing by 3.4e-07 Eh -- the same difference in both, so a systematic property of the
two inversions rather than noise. It costs 5.5 per cent of a direct calculation and
6.8 of a held one, which answers the question of what multiplying by the root costs
against solving the factor: not much.

## Against pyscf, where the core Hamiltonian was the thing being wrong

Everything above compares the SIMD RI-JK driver against VeloxChem's own conventional
build. That comparison was worth very little above g, and this section is how that was
found out.

Water at the geometry below, `conv_thresh = 1e-8` here and `conv_tol = 1e-12` in pyscf.
The AO map is built from the quantum numbers of the two labellings and **verified on
the overlap before anything else is compared** -- a wrong map gives a confident wrong
answer. It agrees to 3e-11 at every basis, so the map is not in question anywhere
below.

```
O   0.000000   0.000000   0.117790
H   0.000000   0.755453  -0.471161
H   0.000000  -0.755453  -0.471161
```

### The one-electron integrals, before and after

Largest absolute difference from pyscf, whole matrix:

| basis | nao | overlap | kinetic, SIMD | kinetic, plain | nuclear, SIMD | nuclear, plain |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cc-pVTZ | 58 | 2.8e-11 | 1.1e-10 | 1.1e-10 | 4.7e-10 | 4.7e-10 |
| cc-pVQZ | 115 | 2.8e-11 | 1.7e-10 | 1.7e-10 | 5.3e-10 | 5.3e-10 |
| cc-pV5Z | 201 | 3.1e-11 | 3.0e-10 | **1.5e+01** | 4.8e-10 | **8.4e+00** |

The plain drivers are exact through g and return **zeros for h blocks**. The SCF took
its core Hamiltonian from them, so a basis with h functions was built on a wrong
Hamiltonian with no warning of any kind. Both are now taken from the SIMD drivers,
which are right at h and i.

### What that did not fix, which is the larger half

With the core Hamiltonian correct, the conventional SCF at cc-pV5Z is **still** wrong:

| basis | conventional | pyscf exact | difference |
| --- | ---: | ---: | ---: |
| def2-svp | -75.9609698336 | -75.9609698336 | +1.5e-12 |
| cc-pVTZ | -76.0570982357 | -76.0570982357 | +2.3e-12 |
| cc-pV5Z | -76.0724789681 | -76.0670116535 | **-5.5e-03** |

The four-center driver has no h kernels at all -- the highest in
`ElectronRepulsionFunc.hpp` is `...RecSSSG`. From an identical density, mapped from a
converged pyscf calculation, its matrices relative to pyscf's:

| basis | Coulomb | exchange |
| --- | ---: | ---: |
| cc-pVQZ | 1.4e-11 | 1.1e-11 |
| cc-pV5Z | **4.3e-01** | **1.4e-02** |

**Above g, RI-JK is the only correct route in this code.** Which also means the
comparisons this file made against the conventional build above g were measuring the
conventional build's gaps, not the driver's -- the trap of using a reference whose
coverage is narrower than the thing being checked.

### The driver, against a reference which does reach h

pyscf `RHF().density_fit(auxbasis=...)`, def2-universal-jkfit on both sides, so the
fitting error is common to the two and what remains is the implementation:

| basis | SIMD RI-JK | pyscf DF | difference | pyscf exact | RI error |
| --- | ---: | ---: | ---: | ---: | ---: |
| def2-svp | -75.9609138700 | -75.9609138700 | +1.8e-12 | -75.9609698336 | 5.6e-05 |
| cc-pVTZ | -76.0570953623 | -76.0570953623 | +2.6e-12 | -76.0570982357 | 2.9e-06 |
| cc-pVQZ | -76.0647542452 | -76.0647542452 | +2.4e-12 | -76.0647584041 | 4.2e-06 |
| cc-pV5Z | -76.0670057742 | -76.0670057742 | +2.5e-12 | -76.0670116535 | 5.9e-06 |

**Agreement is 2.5e-12 at every basis, h functions included.** The last column is the
resolution of the identity itself, which is the approximation being made and is four
to seven orders larger than the disagreement between the two implementations of it.

## The molecular gradient, where the diffuse functions decide the ratio

*Superseded for the timings of the resolution of the identity, which are about ten
times what the same calculation costs now. The shape of the argument, and every
four-center number, still holds: "Ninety-two per cent of the gradient was not the
integrals", at the end of this file.*

The RI-JK gradient is wired into `ScfGradientDriver` and this is its first
measurement. Caffeine, 24 atoms, def2-universal-jkfit, one rank of 14 threads on the
M4 Max, at `cd9cb941a`. Records in
`benchmarks/data/gradient/2026-09-16_m4max_caffeine.json`; the suite is
`benchmarks/scripts/grad_laptop.py`.

Gradient wall time alone -- the SCF before it is timed separately, because an RI-JK
calculation has already won on the energy before the gradient starts. Best of two,
both ways in one process per case, so the comparison is never made across runs.

| functional | basis | nao | method | gradient | speedup | SCF | vs four-centre |
| --- | --- | ---: | --- | ---: | ---: | ---: | ---: |
| HF | def2-svp | 246 | four-centre | 8.69 | 1.00 | 12.45 | |
| | | | RI-JK simd, in memory | 2.19 | 3.97 | 1.12 | 7.7e-05 |
| HF | def2-svpd | 366 | four-centre | 42.13 | 1.00 | 50.71 | |
| | | | RI-JK simd, in memory | 2.94 | **14.31** | 2.43 | 7.7e-05 |
| B3LYP | def2-svp | 246 | four-centre | 9.54 | 1.00 | 15.04 | |
| | | | RI-JK simd, in memory | 2.83 | 3.37 | 4.04 | 1.6e-05 |
| B3LYP | def2-svpd | 366 | four-centre | 43.79 | 1.00 | 56.72 | |
| | | | RI-JK simd, in memory | 4.24 | **10.33** | 8.80 | 1.7e-05 |

The repeats were tight throughout -- 42.13 against 42.18, 2.99 against 2.94 -- so none
of this is a single bad sample.

**The ratio is not a property of the driver, it is a property of the basis.** Adding
the diffuse shell costs the four-center gradient 8.69 to 42.13 seconds, a factor of
4.8 for a factor of 1.49 in the basis, which is the fourth power it is built on. It
costs the resolution of the identity 2.19 to 2.94, a factor of 1.3. Quoting a speedup
without the basis beside it says nothing: the same driver is 4.0 and 14.3 in the same
table.

### The functional dilutes the ratio, and does not slow anything down

B3LYP looks worse than Hartree-Fock -- 3.4 and 10.3 against 4.0 and 14.3 -- and the
reason is in the two columns rather than in either of them. Subtracting the
Hartree-Fock row from the B3LYP row of each method gives what the quadrature costs:

| basis | four-centre | RI-JK simd |
| --- | ---: | ---: |
| def2-svp | 0.85 | 0.64 |
| def2-svpd | 1.66 | 1.30 |

It is **very nearly the same work in both rows**, and it does not shrink when the
two-electron part does. At def2-svp it is added to a numerator of 8.69 and a
denominator of 2.19, and a near-constant added to both sides of a ratio pulls it
toward one. The RI-JK gradient is not slower at B3LYP than at Hartree-Fock for
anything it is responsible for -- it is carrying a fixed passenger that the
four-centre path barely notices and it cannot hide.

### What the gradient agrees with, and what it does not

Against the four-center gradient the difference is 7.7e-05 at Hartree-Fock and 1.6e-05
at B3LYP -- the fitting error, consistent with the 5.7e-04 the energies differ by, and
smaller at B3LYP because only a fifth of the exchange is fitted at all.

The sharper test is finite differences of the RI-JK energy itself, which has no
fitting error in it because both sides make the same approximation. def2-svp, central
differences at 1e-4 bohr, largest component:

| molecule | Hartree-Fock | B3LYP |
| --- | ---: | ---: |
| water | 1.2e-09 | 1.4e-06 |
| CH3NH.OH, 6 atoms, no symmetry | 6.8e-09 | 2.2e-06 |

**The B3LYP row is not the gradient.** Run the same finite differences against the
four-center path and it gives 1.389e-06 on water -- the same number to four figures,
where four-center Hartree-Fock gives 1.2e-09 against the resolution of the identity's
1.2e-09. The residual is the quadrature grid and the step, and it belongs to both
paths equally.

The same thing shows in translational invariance, which is a check that costs nothing
and travels with every record. At Hartree-Fock the gradient sums over the atoms to
1e-12. At B3LYP it sums to 1.7e-05 -- and the four-center path sums to 1.73939e-05
where the resolution of the identity sums to 1.73893e-05. An atom-centered grid is not
translationally invariant, and a check that looks like a failure of the driver is a
property of the quadrature that both drivers inherit.

### One rank, and why that is not a temporary omission

These numbers are OpenMP on one rank, and the gradient refuses to run on more. The
Fock build tolerates B vectors spread over the ranks because the factor of the metric
is folded into them, so the Coulomb matrix is a sum over the auxiliary basis which
factorizes and the ranks simply add their shares. The gradient contracts the
derivatives of the integrals themselves, so it needs the fitting coefficients in the
basis of those integrals, and the transposed factor which carries them there reaches
across the whole auxiliary basis. A rank holding a share of it cannot form them. The
right-hand side has to be complete before the solve -- the same coupling the direct
Fock build already handles with an explicit reduction before `solve_fitting`, and the
same thing a distributed gradient will have to do.

## Geometry optimization, where the step count is half the ratio

The gradient is wired into `ScfGradientDriver`, so a geometry optimization can be
carried with the resolution of the identity end to end. Caffeine from the geometry
in `benchmarks/geometries`, def2-universal-jkfit, one rank of 14 threads on the M4
Max, run to convergence with no cap on the iterations, at `ac2b22a5f`. The working
tree was dirty for the run: what was uncommitted were the benchmark scripts
themselves, which are now in `benchmarks/scripts/opt_laptop.py`, with the records in
`benchmarks/data/optimization/2026-09-16_m4max_caffeine.json`. Two hours and five
minutes of machine for the eight optimizations.

| functional | basis | nao | method | total | speedup | steps | s/step | energy |
| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| HF | def2-svp | 246 | four-centre | 619.8 | 1.00 | 33 | 18.78 | -675.83732476 |
| | | | RI-JK simd, in memory | 90.2 | 6.87 | 27 | 3.34 | -675.83673457 |
| HF | def2-svpd | 366 | four-centre | 2961.8 | 1.00 | 35 | 84.62 | -675.86719050 |
| | | | RI-JK simd, in memory | 174.2 | **17.01** | 35 | 4.98 | -675.86660602 |
| B3LYP | def2-svp | 246 | four-centre | 628.5 | 1.00 | 29 | 21.67 | -679.88509552 |
| | | | RI-JK simd, in memory | 185.6 | 3.39 | 29 | 6.40 | -679.88514659 |
| B3LYP | def2-svpd | 366 | four-centre | 2471.5 | 1.00 | 27 | 91.53 | -679.92921001 |
| | | | RI-JK simd, in memory | 322.4 | **7.67** | 27 | 11.94 | -679.92926228 |

### Why the record keeps the steps

An optimization is not one measurement, it is a step count times a cost per step,
and **only one of those two is the driver's doing**. Three of the four pairs
converged in the same number of steps as each other, and for those the total
speedup and the per-step speedup are the same number to three figures: 17.01 and
16.99, 3.39 and 3.39, 7.67 and 7.67. Those are clean.

The fourth is not. Hartree-Fock in def2-svp took the four-centre path 33 steps and
the resolution of the identity 27, so its 6.87 is a per-step speedup of **5.62**
multiplied by the optimizer happening to take a shorter route over a slightly
different surface. That is luck and not merit, and it could as easily have gone the
other way. A table which carried only the total would have reported the largest
Hartree-Fock def2-svp speedup in this file and been wrong about where it came from.

The shape of the rest is what the single-point gradients already said: the diffuse
shell decides the ratio, 6.9 to 17.0 at Hartree-Fock, and B3LYP is lower only
because the quadrature is a fixed cost in both columns.

### An outlier which was a methyl group

That same Hartree-Fock def2-svp row put an atom 4.57e-02 bohr away from where the
four-centre optimization put it -- sixty-six times the displacement of any other
row, and suspicious in exactly the row whose step count already disagreed.

It is not a structural disagreement. The four atoms which moved are H17, H15, H16
and H21, which is one methyl group, and the heavy atom framework agrees to 7.46e-05
bohr, indistinguishable from every other row:

| | four largest movers | heavy-atom bonds agree to |
| --- | --- | ---: |
| HF, def2-svp | H17, H15, H16, H21 | 7.46e-05 |
| HF, def2-svpd | H17, H16, H22, H14 | 7.97e-05 |
| B3LYP, def2-svp | H18, H20, H19, O11 | 3.86e-05 |
| B3LYP, def2-svpd | H15, H17, C9, H16 | 3.40e-05 |

Caffeine has three methyl groups and a methyl rotation costs almost nothing, so the
two surfaces put the rotor at slightly different angles for no energy worth
measuring -- and the optimizer spent its six extra steps chasing that flat
direction. **A displacement is not a disagreement until it is a heavy atom.** The
check is one line and it turns the alarming number in the table into the
uninteresting one it actually is.

### What the two surfaces differ by

| | four-centre | RI-JK simd | difference |
| --- | ---: | ---: | ---: |
| HF, def2-svp | -675.83732476 | -675.83673457 | +5.90e-04 |
| HF, def2-svpd | -675.86719050 | -675.86660602 | +5.84e-04 |
| B3LYP, def2-svp | -679.88509552 | -679.88514659 | -5.11e-05 |
| B3LYP, def2-svpd | -679.92921001 | -679.92926228 | -5.23e-05 |

The fitting error, at each method's own minimum rather than at a common geometry.
It is an order of magnitude smaller at B3LYP, where only a fifth of the exchange is
fitted at all, and it changes sign there: the fitted B3LYP minima lie **below** the
four-centre ones, which the Hartree-Fock rows do not.

## Four bases of the gradient, and the exponent each way is really running at

*Superseded. The exponents fitted here for the resolution of the identity are
distorted by a constant which was later removed, and its timings are nine-tenths
overhead. The four-center numbers stand: "Ninety-two per cent of the gradient was
not the integrals", at the end of this file.*

The gradient section above measured def2-svp and def2-svpd and concluded that the
basis decides the ratio. This extends the same table to the triple zeta pair, which
is enough points to fit an exponent instead of reasoning from two. Caffeine,
def2-universal-jkfit, one rank of 14 threads on the M4 Max, best of two, both ways
in one process per case, at `15b33b62e`. The records are in
`benchmarks/data/gradient/2026-09-16_m4max_caffeine.json` and the table and its
scaling page are rendered from them by `benchmarks/scripts/render_runs.py`.

The def2-svp and def2-svpd rows were re-measured rather than carried over, and they
reproduce the earlier run to within one or two per cent, so the two halves of this
table are comparable.

| functional | basis | nao | four-centre | RI-JK simd | speedup | vs four-centre |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| HF | def2-svp | 246 | 8.61 | 2.22 | 3.87 | 7.7e-05 |
| HF | def2-svpd | 366 | 41.68 | 2.91 | 14.33 | 7.7e-05 |
| HF | def2-tzvp | 494 | 136.34 | 4.13 | 33.00 | 7.1e-05 |
| HF | def2-tzvpd | 614 | 352.64 | 5.46 | **64.61** | 7.3e-05 |
| B3LYP | def2-svp | 246 | 9.26 | 2.78 | 3.33 | 1.6e-05 |
| B3LYP | def2-svpd | 366 | 42.80 | 4.09 | 10.47 | 1.7e-05 |
| B3LYP | def2-tzvp | 494 | 138.23 | 5.75 | 24.06 | 6.2e-05 |
| B3LYP | def2-tzvpd | 614 | 356.52 | 7.92 | **45.04** | 6.3e-05 |

Cost fitted as proportional to nao to the p, over all four bases:

| functional | method | p |
| --- | --- | ---: |
| HF | four-centre | 4.04 |
| HF | RI-JK simd | **0.98** |
| B3LYP | four-centre | 3.97 |
| B3LYP | RI-JK simd | **1.13** |

**The four-center gradient is at its textbook exponent and the resolution of the
identity is linear.** Four point oh four and three point nine seven are what a sum
over four indices costs, with no screening benefit visible across this range. Nought
point nine eight is the other side of it: over a factor of two and a half in the
basis the gradient of the resolution of the identity went from 2.22 seconds to 5.46.
The B3LYP exponent is higher at 1.13 because the quadrature is in that column and
does grow, which is the same effect that lowers its speedups.

### The exponent is in one dimension only, and the table says so

One fitting set serves every row: naux is 1242 in all eight of them while nao runs
246 to 614. **So the exponent measured for the resolution of the identity is in the
orbital dimension alone, and it is not the scaling of the method with the problem.**
The ratio between the columns widens from 3.9 to 64.6 partly because the denominator
is being held still, and a reader who takes 64.6 as a trend and extrapolates it will
be wrong. The naux column is in the rendered table for exactly this reason: the
constant is visible beside the thing that is growing.

What the numbers do support is narrower and still worth having. For a fixed fitting
set, which is how these calculations are actually run, enlarging the orbital basis
costs the four-center gradient a fourth power and costs this driver a first power.

### And they do not survive the molecule growing either

*Superseded. The collapse to 1.23 recorded here was one unthreaded loop and not a
property of the method; it is 11.51 now and the two grow alike: "Ninety-two per
cent of the gradient was not the integrals", at the end of this file.*

The exponents above are scaling **with the basis at a fixed geometry**. They are not
scaling with the size of the molecule, and the run which tested that broke both of
them at once. Tagrisso, 70 atoms against caffeine's 24, in the same def2-svp and with
the same fitting set, so nao goes 246 to 683 and naux 1242 to 3387:

| | caffeine | tagrisso | observed | the exponent predicts |
| --- | ---: | ---: | ---: | ---: |
| four-centre | 8.61 | 114.17 | 13.3x | 61.9x |
| RI-JK simd | 2.22 | 92.92 | **41.9x** | 2.7x |

**Both predictions are wrong, in opposite directions, for different reasons.** The
four-center gradient came in nearly five times cheaper than a fourth power says,
because seventy atoms spread out is a geometry where screening finally has something
to discard and a compact molecule with a bigger basis is not. The resolution of the
identity came in fifteen times dearer than a first power says, because that exponent
was fitted with naux held still and here it nearly tripled.

So the speedup goes with it. The gradient is 3.87 on caffeine in def2-svp and **1.23
on tagrisso in the same basis**, 3.33 and 1.22 at B3LYP. The advantage measured in
the table above belongs to a compact molecule with a large basis, which is the shape
this file has been measuring all along, and it does not carry to a large molecule
with a small one.

Two things do carry. The gradient is still right -- it agrees with the four-center
one to 6.0e-05 at Hartree-Fock and 1.7e-05 at B3LYP, the same fitting error as
caffeine, at 3387 auxiliary functions and 3.83 GB of B vectors. And the energy still
wins outright, 151 seconds against 30, so a single point and its gradient together
are 2.16 times quicker at Hartree-Fock and 1.88 at B3LYP. **It is the gradient
specifically, and at this shape of problem specifically, that has lost its lead** --
which is a statement about where to look next, not a retraction of the table above.

### One column that does not behave, and is not yet explained

The last column is the largest disagreement with the four-center gradient. At
Hartree-Fock it is flat across the whole range -- 7.7, 7.7, 7.1, 7.3, all e-05 --
which is what a fitting error should do. At B3LYP it is not:

| | def2-svp | def2-svpd | def2-tzvp | def2-tzvpd |
| --- | ---: | ---: | ---: | ---: |
| B3LYP, observed | 1.6e-05 | 1.7e-05 | 6.2e-05 | 6.3e-05 |
| a fifth of the Hartree-Fock row | 1.5e-05 | 1.5e-05 | 1.4e-05 | 1.5e-05 |

At double zeta B3LYP sits where a functional which fits a fifth of its exchange
should sit. At triple zeta it is four times that, and it steps rather than drifts:
flat, then a jump at the double to triple zeta boundary, then flat again. A fitting
error which is a fifth of another fitting error should not do that.

The likely candidate is the quadrature rather than the fit -- the two paths converge
to slightly different densities, so their exchange-correlation gradients are not
quite the same number, and that difference is not scaled by the fraction of exact
exchange -- but **this has not been checked and is written here as a question, not a
finding.** It is 6e-05 on a gradient whose largest component is order 0.1, so it
changes nothing about the numbers above; it is recorded because a column which steps
where nothing else does is worth returning to.

## Ninety-two per cent of the gradient was not the integrals

Everything above about the gradient measured a driver in which the derivative
integrals were under two per cent of the time. This section is the profile that
found that out, what was changed, and the numbers the two sections above have to
be read against now.

### The profile

Tagrisso, def2-svp, 70 atoms, 683 orbital and 3387 auxiliary functions. The phases
timed through the bindings the driver already exposes, replicating its own loop
rather than instrumenting it:

| phase | time | share |
| --- | ---: | ---: |
| forming the fitted densities | 84.4 s | **92%** |
| contraction, by remainder | 5.1 s | 5.6% |
| the three-center derivative integrals | 1.75 s | 1.9% |
| the per-atom sparsity patterns | 0.02 s | -- |
| the two-center (P\|Q) term | 0.01 s | -- |

**The SIMD derivative integrals, which is what the kernels were written for, were
one part in fifty of the gradient.** Everything this file has measured about them
was measuring a thing that was not the cost.

Inside that phase were two steps, each of them the square of the auxiliary basis
times the square of the orbitals, and neither of them threaded. The processor trace
says it plainly on fourteen cores: two hundred and seventy-seven per cent falling
to a hundred and ninety, and then a tail at ninety-nine. One core, for the last
third of it.

### What was changed

Four commits, no kernel touched:

| | fitted densities |
| --- | ---: |
| before | 84.4 s |
| expanding the metric once per phase and not once per element | 69.3 s |
| applying the transposed factor in one multiply over all elements | 27.6 s |
| taking the Gram product of the fitted densities as one multiply | **1.2 s** |

The first was a dense expansion of the metric, ninety-two megabytes, formed inside
a loop that ran once per element of a matrix. The second and third were a matrix
times a matrix written as a sum: one call to the library in place of eight thousand
calls of ours, and a Gram product in place of a quadruple loop. The elements are
taken in panels so the arrays are bounded by the budget and not by the problem,
which costs nothing -- one panel and five hundred and fifty-seven panels differ by
under two per cent and agree to 1e-12.

Every gradient is unchanged. They reproduce the measurements above to 1.2e-12 and
1.7e-12, the agreement with the four-center gradient does not move in any row, and
water against finite differences is back to 1.157e-09, the figure it had before any
of this.

### The corrected tables

Caffeine, gradient wall time, best of two, every four-center row re-measured in the
same session as a control and every one of them reproducing to between 0.05 and 2.7
per cent:

| functional | basis | four-centre | RI-JK simd | speedup | was |
| --- | --- | ---: | ---: | ---: | ---: |
| HF | def2-svp | 8.84 | 0.74 | 11.95 | 3.87 |
| HF | def2-svpd | 42.11 | 1.36 | 30.96 | 14.33 |
| HF | def2-tzvp | 136.48 | 2.52 | 54.16 | 33.00 |
| HF | def2-tzvpd | 352.46 | 3.84 | **91.79** | 64.61 |
| B3LYP | def2-svp | 9.17 | 1.25 | 7.34 | 3.33 |
| B3LYP | def2-svpd | 42.79 | 2.53 | 16.91 | 10.47 |
| B3LYP | def2-tzvp | 138.15 | 4.01 | 34.45 | 24.06 |
| B3LYP | def2-tzvpd | 353.66 | 6.64 | **53.26** | 45.04 |

And tagrisso, which is the row that mattered:

| functional | four-centre | RI-JK simd | speedup | was |
| --- | ---: | ---: | ---: | ---: |
| HF | 114.50 | 9.95 | **11.51** | 1.23 |
| B3LYP | 117.10 | 13.01 | **9.00** | 1.22 |

### Two conclusions above are now wrong, and this is how

**The exponent.** The section above records the gradient of the resolution of the
identity scaling as nao to the 0.98, and calls it linear. It is not:

| | recorded | now |
| --- | ---: | ---: |
| HF, four-centre | 4.04 | 4.01 |
| HF, RI-JK simd | **0.98** | **1.82** |
| B3LYP, four-centre | 3.97 | 3.97 |
| B3LYP, RI-JK simd | **1.13** | **1.78** |

The four-center exponents do not move, as they cannot. The other two nearly doubled,
and the reason is instructive rather than embarrassing: across the four bases of
caffeine, naux is 1242 and the occupied orbitals are 51 in every one of them, and
only nao grows. The term which cost the square of each was therefore **a constant
of about one and a half seconds added to every row**, and a constant added to a
power law flattens it. Take 1.5 off the recorded series and it reads 0.72, 1.41,
2.63, 3.96, which is the new series. **A fitted exponent is only the exponent of
the thing that varies; a large constant in the same column reads as a smaller
power.**

**The molecule.** The section above concludes that the advantage "does not carry to
a large molecule with a small one", on the evidence that caffeine's 3.87 became
tagrisso's 1.23. That conclusion was measuring the same constant, which is not
constant between molecules: naux and the orbitals both grow with the molecule, so
the term grew fifty-fold where everything else grew thirteen-fold. With it gone the
two methods grow at the same rate from caffeine to tagrisso in def2-svp:

| | caffeine | tagrisso | factor |
| --- | ---: | ---: | ---: |
| four-centre | 8.84 | 114.50 | 12.95x |
| RI-JK simd | 0.74 | 9.95 | 13.45x |

Thirteen and thirteen, where it was thirteen and forty-two. The divergence was the
Gram product and nothing about the method. End to end a single point and its
gradient on tagrisso is now 6.69 times quicker, where that section recorded 2.16.

**What both mistakes have in common** is that the measurement was sound and the
attribution was not. Every number in those sections is reproducible and none of
them has been withdrawn. What was wrong was reading a curve without asking which
of its terms was moving -- which a profile answers in twenty minutes and four
tables of timings do not answer at all.

## The excited states, where the ratio reaches a hundred

The resolution of the identity is wired into the Tamm-Dancoff approximation, by way
of a driver which is not the one the self consistent field uses. The two are asked
for different things: the field has one density per build, symmetric and idempotent,
and forms its exchange from the occupied orbitals alone, where a response
calculation has a batch of densities per build, none of them symmetric and each of
them living between the occupied orbitals and the virtual ones.

What makes it cheap is that such a density arrives already factorised. A trial
vector gives C(occupied) Z C(virtual) transposed, and with B(q) symmetric the
exchange of a density left times right transposed is

    K = sum over q of (B(q) left) (B(q) right) transposed

so both halves are the transformation the field already had, with different
coefficients. The virtual space never appears -- C(virtual) Z transposed is the
basis by the occupied orbitals -- and the left factor is the ground state occupied
orbitals for every trial vector of the batch, transformed once for all of them.

Caffeine, five states, def2-universal-jkfit, one rank of 14 threads on the M4 Max,
at `dfc9f00c6`. Records in `benchmarks/data/tda/2026-09-16_m4max_caffeine.json`.
Two hours and twenty minutes for the sixteen rows.

| functional | basis | nao | four-centre | RI-JK simd | speedup | iter | s/iter |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| HF | def2-svp | 246 | 58.50 | 3.37 | 17.34 | 15 | 3.90 to 0.22 |
| HF | def2-svpd | 366 | 289.93 | 7.81 | 37.13 | 17 | 17.05 to 0.46 |
| HF | def2-tzvp | 494 | 914.97 | 13.16 | 69.52 | 16 | 57.19 to 0.82 |
| HF | def2-tzvpd | 614 | 2554.92 | 23.22 | **110.03** | 18 | 141.94 to 1.29 |
| B3LYP | def2-svp | 246 | 51.64 | 12.24 | 4.22 | 10 | 5.16 to 1.22 |
| B3LYP | def2-svpd | 366 | 211.23 | 29.41 | 7.18 | 10 | 21.12 to 2.94 |
| B3LYP | def2-tzvp | 494 | 646.65 | 40.26 | 16.06 | 10 | 64.66 to 4.03 |
| B3LYP | def2-tzvpd | 614 | 1509.78 | 64.36 | **23.46** | 9 | 167.75 to 7.15 |

**Every pair converged in the same number of iterations as its partner**, which is
what makes each ratio a comparison of speed and not of luck. The iteration column is
in the table for that reason and is worth reading before the speedup column: it is
the thing which, when it differs, quietly turns a ratio into something else.

### The exponents, and a check across two calculations

Cost per iteration fitted against the orbital basis:

| | four-centre | RI-JK simd |
| --- | ---: | ---: |
| HF | 3.93 | **1.91** |
| B3LYP | 3.78 | **1.84** |

The gradient, measured the same day and through entirely separate code, gives 4.01
and 1.82 at Hartree-Fock. **Two calculations which share nothing but the driver
underneath agree on what that driver costs**, which is worth more than either number
alone. The iteration counts vary from fifteen to eighteen across the Hartree-Fock
series, so the per-iteration unit carries some noise from the batch of trial vectors
changing size; fitting the totals instead gives 4.08 and 2.06, the same picture.

### The fitting error falls as the basis grows

| | def2-svp | def2-svpd | def2-tzvp | def2-tzvpd |
| --- | ---: | ---: | ---: | ---: |
| HF | 2.2e-05 | 1.7e-05 | 1.3e-05 | 7.4e-06 |
| B3LYP | 6.0e-06 | 5.9e-06 | 5.6e-06 | 4.1e-06 |

The largest disagreement in an excitation energy, in hartree, which is four
hundredths of a millielectronvolt at its worst. It **improves** with the orbital
basis, and for the reason the gradient section gives in reverse: one fitting set
serves every row, and it is a better fit to a larger orbital basis than to a small
one.

### Why B3LYP is lower, which these numbers do not establish

B3LYP reaches 23 where Hartree-Fock reaches 110, and the natural reading is the
exchange-correlation quadrature: it is the same work in both columns, it does not
shrink when the two-electron part does, and a constant added to both sides of a
ratio pulls it toward one -- which is what the gradient section found for the same
functional.

**The reading is probably right and these numbers do not show it.** The obvious way
to extract the quadrature is to subtract the Hartree-Fock cost per iteration from
the B3LYP one at the same basis, and that gives two answers which disagree:

| | from the RI-JK column | from the four-centre column |
| --- | ---: | ---: |
| def2-svp | 1.00 s | 1.26 s |
| def2-svpd | 2.48 s | 4.07 s |
| def2-tzvp | 3.20 s | 7.48 s |
| def2-tzvpd | **5.86 s** | **25.81 s** |

A factor of four and a half apart at the largest basis. The subtraction is not
valid: Hartree-Fock converges in eighteen iterations there and B3LYP in nine, so the
two runs do not carry the same number of trial vectors per iteration and a second
per iteration is not the same unit in the two of them. **A quantity which can be
estimated two ways and gives two answers has been estimated zero ways.** What the
quadrature actually costs wants a profile, which is how the gradient's ninety-two
per cent was found and not something a table of totals can answer.

### The same calculation on a molecule three times the size

Tagrisso, 70 atoms against caffeine's 24, in the same def2-svp and with the same
fitting set, so nao goes 246 to 683, the auxiliary basis 1242 to 3387 and the
occupied orbitals 51 to 133. Five states, at `633e3a8ef`, records in
`benchmarks/data/tda/2026-09-16_m4max_tagrisso.json`.

| functional | four-centre | RI-JK simd | speedup | iter | s/iter | max dE |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| HF | 866.14 | 138.38 | **6.26** | 19 | 45.59 to 7.28 | 7.4e-06 |
| B3LYP | 693.42 | 167.62 | **4.14** | 13 | 53.34 to 12.89 | 5.0e-06 |

**It is right at this size**, which was the open question and not the speedup: the
response driver had never been asked for seventy atoms or three thousand auxiliary
functions before this run. The excitation energies agree with the four-center ones
to 7.4e-06 and 5.0e-06 hartree, which is **tighter** than the same molecule's
caffeine rows, for the reason the caffeine section gives -- one fitting set is a
better fit to a larger orbital basis. Both pairs converged in the same number of
iterations as their partners.

### Why the ratio falls, and why that is not the gradient's story

Hartree-Fock goes from 17.36 on caffeine to 6.26 here. The gradient section above
records a collapse of exactly this shape, from 3.87 to 1.23, and it turned out to be
an unthreaded loop which should not have been there. **This one is not that.** The
exchange of a factorised density is the auxiliary basis times the square of the
orbital basis times the occupied orbitals, and all three of those grow with the
molecule:

| | caffeine | tagrisso | factor |
| --- | ---: | ---: | ---: |
| naux times nao squared times nocc | 3.83e+09 | 2.10e+11 | **54.8x** |
| observed, per iteration | 0.22 s | 7.28 s | 32.4x |
| four-centre, per iteration | 3.90 s | 45.59 s | 11.7x |

The resolution of the identity grew by **less** than its own cost model says, so the
batching of the trial vectors is earning something. The four-center path grew by far
less than a fourth power, because seventy atoms spread out is where screening
finally has distant pairs to discard. Two methods with honest costs, moving apart.
The lesson of the gradient was to ask which term is growing before believing a
curve; asking it here gives an answer that exonerates the implementation instead of
indicting it.

End to end, a ground state and its five excited states is 6.05 times quicker at
Hartree-Fock and 3.94 at B3LYP.

### The two functionals converge, which is worth noticing and not yet explaining

| | caffeine | tagrisso |
| --- | ---: | ---: |
| HF | 17.36 | 6.26 |
| B3LYP | 4.22 | 4.14 |

Hartree-Fock falls by nearly three and B3LYP barely moves, so what was a fourfold
gap between the functionals on caffeine is under fifty per cent here. A fixed
quadrature cost mattering less as the two-electron work grows would produce exactly
this, and that is the natural reading.

**It is a reading and not a measurement**, for the same reason the caffeine section
sets out: nineteen iterations against thirteen means the two runs do not carry the
same trial vectors per iteration, so a second per iteration is not the same unit in
the two of them and the functionals cannot be subtracted. So it was profiled, and
the next section is the answer.

### What the iteration is actually made of

The solver's own profiler, on the two runs above. No subtraction of anything from
anything:

| | Hartree-Fock, 19 iter | B3LYP, 13 iter |
| --- | ---: | ---: |
| the two-electron build | 128.40 s, **93.7%** | 91.64 s, **54.6%** |
| the quadrature | -- | 66.84 s, **39.8%** |
| forming the B vectors, once | 7.77 s, 5.7% | 7.78 s, 4.6% |
| everything else | under 1 s | under 2 s |

**The quadrature is two fifths of a hybrid's iteration.** The guess this file was
about to record, from the invalid subtraction, was four fifths. It was wrong by
half, and wrong in the direction that would have sent the next piece of work at the
wrong target.

What makes the split believable is a number which appears twice: the two-electron
build costs 6.76 seconds an iteration at Hartree-Fock and 7.05 at B3LYP. **Those
should be equal** -- the exchange is formed whole and then scaled, so a fifth of it
costs exactly what all of it costs -- and two runs which share no timing apparatus
agree on it to four per cent. A quantity measured twice by accident is worth more
than one measured once on purpose.

### Inside the two-electron build

Timing the driver's two entry points against a batch of trial vectors of the size
the solver actually forms:

| batch | the whole build | the exchange | the Coulomb and the assembly | per density |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 2.14 s | 1.99 s | 0.15 s | 2.14 s |
| 5 | 6.19 s | 6.01 s | 0.18 s | 1.24 s |
| 10 | 11.26 s | 10.85 s | 0.41 s | 1.13 s |
| 20 | 21.98 s | 21.15 s | 0.84 s | 1.10 s |

**The exchange is ninety-three to ninety-six per cent of the build**, and the build
is ninety-four per cent of a Hartree-Fock calculation, so nine parts in ten of the
whole thing is one routine. The Coulomb, which the resolution of the identity is
usually introduced for, is under a tenth of it.

The last column is the batching working. The left factor is the ground state
occupied orbitals and is transformed once for the whole batch, so its cost is fixed
per batch and divided among the densities: half the time at one density, a few per
cent at twenty, and flat from ten onward. That is the one thing the response driver
does which the field's driver has no reason to.

### Inside the exchange, as far as it can be seen from outside

The driver's loop replicated through the bindings it calls, so that the two halves
are timed where they stand and nothing in the C++ is instrumented:

| batch | the exchange | the transforms | the rest, by remainder |
| ---: | ---: | ---: | ---: |
| 5 | 6.00 s | 3.07 s | 2.93 s |
| 20 | 21.32 s | 10.20 s | 11.12 s |

The totals reproduce the table above to within one per cent, so the replication is
the same work. **The transform side is nevertheless too large**, and by a knowable
amount: compute_w_vectors returns its matrices by value and pybind copies each of
them into a python object, which at a batch of twenty is some four gigabytes across
the boundary for every batch of the auxiliary basis and about fifty over the run --
seconds of memcpy which the driver calling itself never pays.

So what can be said is that **the transform is between a third and a half of the
exchange and the accumulation is the rest, and neither of them dominates**. That is
the useful part. There is no single routine here holding nine tenths of the time,
which is what the gradient turned out to have and what makes a profile worth
running: the two halves are both already products of matrices, and both are doing
work the cost model asks for.

The consequence is a negative result worth writing down. **There is no cheap win
left in this code.** The one redundancy -- the left factor is the ground state
occupied orbitals and is transformed again at every iteration of the Davidson,
though it never changes -- is worth one part in one plus the batch, a few per cent
where the batch is twenty, and would cost two and a half gigabytes to hold. What
remains is algorithmic: a smaller fitting set, or something which screens the
transformed vectors, neither of which is a change to this routine.

Pinning the split closer wants timing inside the driver. The technique which
answered the gradient in twenty minutes runs into the binding here, which is worth
knowing about the technique.

## Linear response, and what the second term of a density costs

The Tamm-Dancoff approximation drops the de-excitation block, so its trial vector
gives a density of one term. A full linear response vector gives two, and they are
not each other's transpose:

    D = C(occupied) (-Z) C(virtual) transposed + C(virtual) Y transposed C(occupied) transposed

The second has the virtual orbitals on the left, which is the expensive side to
transform. It need not be: the exchange of a transposed density is the transposed
exchange of it, so the second term is taken as the transpose of one which carries
the occupied orbitals on the left, and the virtual orbitals are never transformed.
Both terms then share one transformation of the occupied orbitals. **And the Coulomb
is taken once for the two of them**, because it sees only the symmetric part of a
density and the two terms together have the same symmetric part as the single term
whose right factor is the sum of theirs.

Caffeine, five states, def2-universal-jkfit, one rank of 14 threads, at `29419eceb`.
Records in `benchmarks/data/tda/2026-09-16_m4max_caffeine_rpa.json`.

| functional | basis | nao | four-centre | RI-JK simd | speedup | iter |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| HF | def2-svp | 246 | 58.19 | 5.32 | 10.94 | 17 |
| HF | def2-svpd | 366 | 281.97 | 11.75 | 24.00 | 19 |
| HF | def2-tzvp | 494 | 935.20 | 20.89 | 44.76 | 19 |
| HF | def2-tzvpd | 614 | 2666.00 | 39.17 | **68.06** | 22 |
| B3LYP | def2-svp | 246 | 49.25 | 12.94 | 3.81 | 11 |
| B3LYP | def2-svpd | 366 | 205.82 | 31.90 | 6.45 | 11 |
| B3LYP | def2-tzvp | 494 | 656.47 | 46.16 | 14.22 | 12 |
| B3LYP | def2-tzvpd | 614 | 1639.05 | 80.01 | **20.49** | 11 |

Every pair converged in the same number of iterations as its partner. The excitation
energies agree with the four-center ones to between 4.5e-06 and 2.3e-05 hartree, and
improve with the basis as the Tamm-Dancoff ones do.

### Measuring the second term without comparing two runs

The obvious way to price the second term is to divide the linear response time by
the Tamm-Dancoff one. **That comparison is not available**: the two runs converge in
different numbers of iterations -- seventeen against fifteen in the first row alone
-- so they do not carry the same trial vectors and their times are not of the same
thing. This file has already recorded one quantity ruined that way.

What is available is the **ratio of two ratios**. A speedup is measured inside one
run between two methods which took the same iterations, so it is clean; and to the
four-center path a density is a density, its cost not depending on whether it was
made of one term or two. So the speedup falls by exactly what the second term costs
the resolution of the identity, and nothing else:

| functional | basis | Tamm-Dancoff | linear response | what the term cost |
| --- | --- | ---: | ---: | ---: |
| HF | def2-svp | 17.34 | 10.94 | 1.59 |
| HF | def2-svpd | 37.13 | 24.00 | 1.55 |
| HF | def2-tzvp | 69.52 | 44.76 | 1.55 |
| HF | def2-tzvpd | 110.03 | 68.06 | 1.62 |
| B3LYP | def2-svp | 4.22 | 3.81 | 1.11 |
| B3LYP | def2-svpd | 7.18 | 6.45 | 1.11 |
| B3LYP | def2-tzvp | 16.06 | 14.22 | 1.13 |
| B3LYP | def2-tzvpd | 23.46 | 20.49 | 1.15 |

**A second term costs about three fifths of a first one, and not a whole one.** Two
exchanges are formed where one was, so the naive price is two; the measured price is
1.55 to 1.62 across a factor of two and a half in the basis. The difference is the
sharing: one Coulomb for both terms and one transformation of the occupied orbitals
for the whole batch.

At B3LYP it costs 1.11 to 1.15, which is the quadrature again. It does not care how
many terms a density has, so where it is a large part of the iteration the second
term is nearly free. **The same structure that caps the speedup at B3LYP also makes
B3LYP the place where linear response is cheapest over Tamm-Dancoff.** One fact,
cutting both ways.

### The scaling does not move

| | four-centre | RI-JK simd |
| --- | ---: | ---: |
| HF | 3.90 | 1.89 |
| B3LYP | 3.78 | 1.88 |

Against 3.93 and 1.91, and 3.78 and 1.84, for the Tamm-Dancoff approximation. Adding
a term to the density changes the constant in front and not the power, which is what
it should do and is worth having measured rather than assumed.

## Three more solvers, which needed no code

The linear response wiring was put into `_e2n_half_size_single_comm`, and that
method is not the eigensolver's alone: the polarizability solver, the damped one
and the C6 driver all reach the Fock build through it. Wiring it for excitation
energies wired it for all of them. This section is the check that this is true and
not merely plausible.

Caffeine and water, def2-svp, one rank of 14 threads. Frequencies zero and 0.1, and
for the damped solver a damping of 0.004556.

| solver | molecule | functional | four-centre | RI-JK simd | speedup | largest difference |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| polarizability | water | HF | 0.06 | 0.06 | 1.08 | 5.2e-04 |
| polarizability | water | B3LYP | 0.11 | 0.10 | 1.07 | 2.2e-04 |
| polarizability | caffeine | HF | 43.98 | 3.74 | **11.76** | 3.1e-03 |
| polarizability | caffeine | B3LYP | 45.74 | 11.59 | 3.95 | 6.0e-03 |
| damped | water | HF | 0.07 | 0.07 | 1.03 | 5.2e-04 |
| damped | water | B3LYP | 0.11 | 0.11 | 1.03 | 2.2e-04 |
| damped | caffeine | HF | 59.19 | 4.85 | **12.21** | 3.1e-03 |
| damped | caffeine | B3LYP | 57.36 | 14.29 | 4.01 | 6.0e-03 |
| C6 | water | HF | 0.14 | 0.11 | 1.19 | 5.1e-04 |
| C6 | water | B3LYP | 0.21 | 0.19 | 1.08 | -- |
| C6 | caffeine | HF | 94.00 | 8.52 | **11.03** | -- |
| C6 | caffeine | B3LYP | 100.04 | 24.58 | 4.07 | -- |

All eighteen components of each tensor, three directions by two frequencies, and for
the damped solver the imaginary parts with them: -0.006699 against -0.006697 and
-0.014854 against -0.014854 on water. The speedups sit where the excitation energies
of the same molecule and basis put them, 10.94 and 3.81, which is what a calculation
solving the same equations at fixed frequencies should give.

**The difference is the evidence, not a worry.** Caffeine's components are about 130
atomic units, so three thousandths is two parts in a hundred thousand -- the fitting
error, and smaller in relative terms than water's. Had the guard quietly sent these
solvers down the dense path the two columns would have agreed to 1e-14, and the
agreement would have proved nothing.

Water's speedup of about one is not a failure either. Twenty-four basis functions is
far below where resolving the identity pays for itself, and the whole calculation is
six hundredths of a second.

### The C6 coefficients, which are the strictest of the four

A C6 coefficient is a Gauss-Legendre quadrature over polarizabilities at imaginary
frequencies: five points by three directions is fifteen independent response solves,
each of them feeding a numerical integration. **An error with a sign to it would
accumulate through that rather than cancel**, which is what makes the coefficient a
better test than the components it is made of.

| | four-centre | RI-JK simd | difference | relative |
| --- | ---: | ---: | ---: | ---: |
| water, HF | 16.74651323 | 16.74405803 | -2.46e-03 | 1.47e-04 |
| water, B3LYP | 17.33060726 | 17.33133502 | +7.28e-04 | 4.20e-05 |
| caffeine, HF | 5255.33816139 | 5255.25395011 | -8.42e-02 | **1.60e-05** |
| caffeine, B3LYP | 5680.62292412 | 5680.85962728 | +2.37e-01 | 4.17e-05 |

Caffeine at sixteen parts in a million is the closest agreement of any property
measured today, and closer than the response components the quadrature is built
from. The differences also change sign between the rows -- minus, plus, minus, plus
-- which is what an error without a bias looks like and is the reason the quadrature
does not make things worse.

### A guard which does nothing, deliberately

The factors are formed only for real trial vectors. **That test passes every time it
is reached**: every solver which gets here works in real arithmetic, the complex
vectors of the damped response being carried as real blocks, and the complex path is
left for future development. So the test is dead code today.

It is kept, and the comment beside it now says why rather than implying complex
vectors are a live case. A driver which takes doubles should not be handed complex
data by a caller which has quietly changed underneath it, and the cost of the test
is one call per batch of trial vectors.

## Higher order response, where the density changes shape with the order

A response density is not one thing. Which blocks of it are nonzero alternates with
the order of the perturbation, and the pattern decides what the resolution of the
identity can do with it. Written out in the molecular orbitals, and checked by
building each one and looking rather than by trusting the algebra:

| density | nonzero blocks | rank |
| --- | --- | ---: |
| first order, a trial vector | ov, vo | twice the occupied |
| second order, a commutator of two first-order things | **oo, vv** | occupied, twice the occupied |
| third order | ov, vo | twice the occupied |

The first and the third have the occupied orbitals on one side of each term, which is
the shape the linear response driver already took. **The second does not**: it is
block diagonal, and its virtual block is carried by the virtual orbitals, of which
there are four times as many.

### The measurement which made the plan smaller

The virtual block factorises two ways, and which is cheaper is not obvious. Its
rank is twice the occupied orbitals, so a factorisation of that rank with both
factors different for every density has the fewest operations; taking the virtual
orbitals themselves as a shared left factor has rank nvir, four times more. Counting
operations says the first should win about two to one.

Measured on caffeine in def2-svp, with densities of the right shape and rank:

| densities | four-centre build | the occupied block | the virtual block, per density | the virtual block, shared |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 2.73 s | 0.13 s | 0.43 s | **0.34 s** |
| 10 | 6.84 s | 0.25 s | 1.08 s | **0.88 s** |
| 20 | 13.66 s | 0.42 s | 2.17 s | **1.70 s** |

**The count of operations had it backwards.** Sharing wins at every size, by about a
quarter, because a shared factor is one wide product of matrices where per-density
factors are many narrow ones, and the profile of the exchange had already said the
transformation is half of it. The same lesson as the gradient, in a different place:
the arithmetic is not the cost.

That made the plan smaller than it was going to be. With the virtual orbitals as a
shared factor the block diagonal density is **two calls of the interface which
already existed, added** -- the Fock matrix being linear in the density -- and no new
kernel, no new interface, and no per-density factors were needed at all.

### What it bought

Caffeine, def2-svp, one rank of 14 threads. Every number here is a whole calculation
against the same calculation built the four-center way.

| | four-centre | RI-JK simd | speedup | largest relative difference |
| --- | ---: | ---: | ---: | ---: |
| quadratic response, first hyperpolarizability | 27.81 | 4.06 | **6.84** | 2.6e-04 |
| second harmonic generation, reduced | 99.26 | 10.07 | **9.86** | 9.5e-04 |
| second harmonic generation, full | 103.33 | 10.83 | **9.54** | 9.5e-04 |
| two-photon absorption, reduced | 164.61 | 16.32 | **10.09** | 1.5e-03 |

Three drivers, one form. The second harmonic driver sends six real columns one way
and twelve real and imaginary ones the other, and the reduced two-photon driver sends
six real columns in its first pass and six real and imaginary in its second, so the
factors are built from the same branch which chooses the columns and never from a
rule of their own.

The disagreements are larger than the linear response ones, which are two parts in a
hundred thousand. That is the fitting error compounding: a hyperpolarizability is
built from products of first-order vectors, and a two-photon amplitude from products
of those, so each order carries the error of the one below it and adds its own.

### A setting which is not passed down is a setting which does nothing

The first quadratic response measurement was **1.03**, against six on the two-electron
build measured directly. The Fock build was not the problem. A response driver of
this kind drives linear solvers of its own, and the settings it hands them are a list
written out by name -- which had `ri_coulomb` in it and not `ri_jk`. So the linear
solves, which are most of the calculation, ran the four-center way in both columns
and the ratio measured almost nothing.

Adding three names to the list turned 1.03 into 6.84.

**The same list appears in nine other drivers**, and in each of them a user setting
ri_jk today would get a calculation which honours it in some places and not others,
with no warning of any kind and a ratio near one to show for it. Two of them are
fixed because their turn came. The rest are a trap for whoever measures them next,
which is why it is written here and not only in the commit.

The two-photon transition driver is the case which shows the difference cleanly. With
its linear solves accelerated and its own Fock builds still on the dense path it
reaches **5.63** on caffeine, against the reduced driver's **10.09** with both. Wiring
the outer loop is worth nearly a factor of two, and the inner solves alone are worth
more than half.

## A frequency sweep of a third order property

The reduced two-photon driver is wired, so the whole calculation goes through the
resolution of the identity and not only the linear solves inside it. Caffeine, five
frequencies from 0.050 to 0.150 in steps of 0.025, one rank of 14 threads. The gamma
of every frequency is in the records, in
`benchmarks/data/redtpa/2026-09-16_m4max_caffeine.json`, and drawn on the second
page of the table beside it.

| functional | basis | four-centre | RI-JK simd | speedup | gamma at 0.100 | largest relative difference |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| HF | def2-svp | 521.42 | 54.72 | 9.53 | 13978.0 + 2441.8i | 8.9e-04 |
| HF | def2-svpd | 2465.55 | 118.25 | **20.85** | 19061.8 + 7474.3i | 6.5e-04 |
| B3LYP | def2-svp | 515.02 | 136.57 | 3.77 | 32374.7 + 20396.0i | 5.3e-04 |
| B3LYP | def2-svpd | 2281.42 | 358.93 | **6.36** | 56171.7 + 26439.5i | 6.4e-04 |

The agreement holds across the whole sweep and not only at the frequency the table
can carry. **It holds where the real part passes through zero**, which it does
between 0.125 and 0.150 in three of the four cells: a relative error is at its worst
where the quantity is smallest, and a disagreement of structure rather than of
fitting would show there first. It does not.

### The sweep is cheaper than five calculations

Five frequencies cost **3.17** times one, not five, on the four-center side, and 3.35
on the other. The two scale alike, so the speedup barely moves from the single
frequency measurement -- 9.53 against 10.09 -- and the run cost two hours where three
were budgeted. Worth knowing before costing the next sweep: the frequencies of one
calculation share more than they look like they do.

### The quadrature, measured without subtracting anything

This file has twice recorded a quantity spoiled by subtracting one run from another.
Here the same question answers itself, because **the four-center columns are nearly
the same for the two functionals** -- 521 against 515 seconds, and 2466 against 2281
-- so the whole of the difference in the ratio sits on the other side:

| | Hartree-Fock | B3LYP | B3LYP over Hartree-Fock |
| --- | ---: | ---: | ---: |
| RI-JK, def2-svp | 54.72 | 136.57 | **2.50** |
| RI-JK, def2-svpd | 118.25 | 358.93 | **3.04** |

The two-electron work is identical in the two rows, the exchange being formed whole
and then scaled. So two and a half to three times is what the exchange-correlation
quadrature costs across five frequencies of Fock builds, read off two columns of one
table with nothing inferred. That is the number the two-photon and excitation
sections could not get, and it took a case where the reference happened to cost the
same either way.

## Cubic response, and a form which carries any density at all

The three-time perturbed calculations were left until last because their batch holds
two orders of density at once. Working out how to hand that to the driver turned out
to be the whole of the problem, and the answer made the problem disappear.

### One form for every density

A density in the molecular orbitals, transformed to the atomic ones, is

    D = C(occ) M(oo) C(occ)^T + C(occ) M(ov) C(vir)^T
      + C(vir) M(vo) C(occ)^T + C(vir) M(vv) C(vir)^T

Gather the two terms which carry the occupied orbitals on the left, and the two which
carry the virtual ones, and it is

    D = C(occ) r_a^T + C(vir) r_b^T
    r_a = C(occ) M(oo)^T + C(vir) M(ov)^T
    r_b = C(occ) M(vo)^T + C(vir) M(vv)^T

which is the four-factor shape the second-order densities already used, with the same
two shared halves -- **and it holds whatever blocks of M are nonzero**. Checked to
1e-14 against a density which is block diagonal and against one which is not.

So a batch which mixes the orders needs no telling apart of its densities. It costs
the basis times the orbitals where a density living only between the occupied
orbitals and the virtual ones would cost the basis times twice the occupied, so it is
the general form and not the cheapest one; for this driver two columns in eight pay
that, which is the right trade for one code path.

### The two batches are laid out differently

| | densities in the batch | cut with |
| --- | --- | --- |
| no functional | the two-time and the three-time ones **in one array** | one stride |
| a functional | two arrays | two strides of their own |

The first of those was misread when this was planned -- the Hartree-Fock path was
taken to send the three-time densities alone, and it concatenates both orders and
sends them as one. The correction is what made the general form necessary, since a
batch boundary can fall anywhere in a concatenated array, including between the
orders.

### What it bought

Caffeine, def2-svp, one frequency triple, one rank of 14 threads.

| | four-centre | RI-JK simd | speedup | gamma |
| --- | ---: | ---: | ---: | --- |
| Hartree-Fock | 67.39 | 10.39 | **6.49** | -4090.71 to -4093.39 |
| B3LYP | 111.88 | 42.38 | **2.64** | 47628.6 to 47672.4 |

Gamma agrees to 9.2e-04, and water to 4.5e-04 at Hartree-Fock and 1.1e-03 at B3LYP.
The second pass of the driver, which is two-time perturbed and was wired after the
first, is worth 5.92 to 6.49 and 2.59 to 2.64 -- nine per cent and two. It carries
one density per frequency triple where the first pass carries four, and the ratio
follows the count.

### Two small terms, and why they read as four per cent

The largest relative disagreement of any quantity the driver returns is 4.3e-02, on
the X3 term. That is not the error of anything worth having:

| term | four-centre | RI-JK simd | difference | relative |
| --- | ---: | ---: | ---: | ---: |
| X2 | 60678.4 | 60724.8 | 46.4 | 7.6e-04 |
| gamma | 47628.6 | 47672.4 | 43.8 | 9.2e-04 |
| E3 | -11157.4 | -11148.6 | 8.8 | 7.9e-04 |
| A2 | -2501.83 | -2502.53 | 0.7 | 2.8e-04 |
| A3 | 404.497 | 404.561 | 0.06 | 1.6e-04 |
| **X3** | **397.019** | **380.077** | **16.9** | **4.3e-02** |
| **T4** | **-192.161** | **-185.882** | **6.3** | **3.3e-02** |

X3 and T4 are the two smallest terms of a sum whose largest is sixty thousand, and
their differences are of the same size as everyone else's. **A relative error is a
statement about the denominator as much as the numerator**, and a term which is a
hundred and fiftieth of the sum it belongs to will always read badly by it.

### One thing which does not fit the pattern

Everywhere else in this file, B3LYP costs about what Hartree-Fock costs on the
four-center side, and the ratio falls only because the quadrature sits in the other
column and does not shrink. Here **the four-center calculation itself is sixty-six
per cent dearer at B3LYP** -- 112 seconds against 67 -- which none of the other
properties showed.

That made 2.64 hard to read, so it was profiled. It is two things and not one.

| | Hartree-Fock | B3LYP |
| --- | ---: | ---: |
| the whole calculation | 66.80 s | 111.95 s |
| the four-center build | 65.80 s over 96 calls | 82.37 s over 118 calls |
| **each of those calls** | **0.685 s** | **0.698 s** |
| the quadrature | -- | **24.05 s** |
| sigma builds of the inner solves | **71** | **94** |

**A Fock matrix costs the same at either functional** -- 0.685 against 0.698 seconds
-- so none of the difference is the build being dearer. Of the forty-five seconds
between them, twenty-four are the quadrature, which Hartree-Fock does not pay at all,
and fifteen are twenty-two more Fock matrices: the inner linear solves needed
ninety-four sigma builds at B3LYP against seventy-one, a third more, because the
calculation converged more slowly. The remainder is the quadrature of the nonlinear
part and the setting up.

So half of it is the functional's integration and half is the functional's
convergence, and **neither is visible in a table of totals**, which is why this was
left unexplained rather than guessed at.

### What is left after the two-electron part goes away

The same two calculations with the resolution of the identity:

| | Hartree-Fock, 10.51 s | B3LYP, 42.39 s |
| --- | ---: | ---: |
| the quadrature | -- | **24.67 s, 58%** |
| the two-electron build | 8.79 s, **84%** | 11.42 s, 27% |
| forming the B vectors | 0.73 s | 0.74 s |

**A cubic response calculation at B3LYP is a quadrature calculation now.** Fifty-eight
per cent of it is integrating the functional and a quarter is the thing the
resolution of the identity was brought in for, which is the whole of why the speedup
is 2.64 and not 6.49. The excitation sections guessed at this and could not measure
it; here it reads off one profile.

The B vectors were formed **three times** in both of them, 0.74 seconds here. That
was the driver and the solvers it drives each building their own, which is nothing at
this size and was 7.8 seconds a time on tagrisso.

**That is fixed.** The forming is now skipped where the vectors are already held, a
driver hands its own to the solvers it drives, and the drivers which made a solver
before setting up their integrals now form them first so there is something to hand
over. Counting the builds rather than reading a profile, because the two-photon
drivers nest profilers and the second one refuses to start:

| | times the forming was entered | times it built anything |
| --- | ---: | ---: |
| quadratic response | 3 | **1** |
| second harmonic generation | 3 | **1** |
| two-photon absorption, reduced | 3 | **1** |
| cubic response | 4 | **1** |

The cubic calculation above went from 10.51 to 9.94 seconds by it, and the forming
from 0.735 to 0.238. Small here and not small on a molecule where the transformation
is eight seconds.

## The remaining six drivers, and a check which is not a benchmark

Cubic response was the hard one. After it, six drivers were left whose Fock builds
still went through the four-centre integrals, and every one of them turned out to be
the same two shapes already written: the `second`/`third` dictionary for a three-time
perturbed pass, the four-factor tuple for a two-time one. **No new code was needed in
the solver or in the shared module** -- `_comp_two_el_int` already took both shapes
and both batch layouts, and `rijkresponse` was not touched. What each driver needed
was the factors collected where its densities are made and handed to its Fock call.

| driver | mode | two-time | three-time |
| --- | --- | ---: | ---: |
| two-photon absorption, full | `tpa` / `tpa_ii` | 24 / 6 | 6 |
| third harmonic generation | `thg` / `thg_ii` | 12 / 6 | 6 |
| third harmonic generation, reduced | `thgred` / `thgred_ii` | 6 / 3 | 3 |
| three-photon absorption | `3pa` / `3pa_ii` | 9 / 6 | 6 |
| two-photon transitions | `tpa_quad` | 4 | -- |
| excited state moments | `qrf` | 2 | -- |

Densities per frequency which become Fock matrices; the first-order ones beside them
are there for the quadrature alone.

### Two shapes which differ from the eight before

**Four of them are real.** The reduced third harmonic, three-photon absorption and
the two-photon transition driver build with `fock_flag='real'` and store only real
columns, so each density is one factor and not the two a complex one is carried by.
Read off the column counts rather than assumed.

**One of them transformed in place.** The six two-time densities of the three-photon
quadratic pass were formed straight into the atomic orbitals inside the `multi_dot`
call, with no name in the molecular orbitals to take factors from. They are now named
first and both the factors and the transformation read the one name, rather than the
expression being written twice where two copies could drift apart.

### What was checked, and how it was made to prove something

Water, def2-svp, 24 basis functions, one frequency, one excited state where the
driver needs one. Every value the driver returns, compared against the same
calculation through the four-centre integrals.

An agreement figure alone proves nothing here: a path which silently fell back to the
four-centre integrals would agree perfectly. So each run also reports what actually
happened inside, in two ways -- the densities in every batch the resolution of the
identity built, and a probe on the Fock build itself which says the mode and whether
the factors arrived.

| | densities per batch, predicted | observed | HF | B3LYP |
| --- | --- | --- | ---: | ---: |
| two-photon, full | 30 joined, or 24 then 6 | as predicted, plus 6 | 2.803e-04 | 2.906e-04 |
| third harmonic | 18 joined, or 12 then 6 | as predicted, plus 6 | 2.778e-04 | 2.867e-04 |
| third harmonic, reduced | 9 joined, or 6 then 3 | as predicted, plus 3 | 2.778e-04 | 2.867e-04 |
| three-photon | 15 joined, or 9 then 6 | as predicted, plus 6 | 2.324e-03 | 3.645e-03 |
| two-photon transitions | 4 | 4 | 3.741e-04 | 1.711e-04 |
| excited state moments | 2 | 2 | 3.741e-04 | 1.711e-04 |

"Joined" is the Hartree-Fock layout, "then" the one with a functional; the trailing
number is the two-time pass. For the last two the probe is the plainer evidence: the
four-centre run reports `('tpa_quad', False)` and `('qrf', False)`, the other one
`True`.

Two of those numbers need qualifying rather than reading straight.

**Three-photon absorption looks ten times worse and is not.** Its worst component is
the `xxx` transition moment, 1.0066 at Hartree-Fock and 0.479 at B3LYP where the
largest component of the same tensor is 46 and 61. The absolute difference is 2.3e-03
and 1.7e-03; the components which are actually determined agree to 1.5e-04. A
relative error on a near-zero number is a statement about the number, not the path.

**Every three-photon value came back sign-flipped.** The phase of an excited state
vector is arbitrary and the two runs pick it differently, so anything odd in that
vector changes sign and means the same thing. The check compares magnitudes. The same
appears in the two-photon transition driver, which reports three flips out of nine.

### What this is not

**There is no timing here on purpose.** Twenty-four basis functions is a correctness
size, not a benchmark size; the wall clocks came out between 0.97 and 1.65 and that
range is scheduling noise, not a result. Nothing in this section should be read as a
speedup, and the tables above deliberately have no such column. What these drivers
cost on a real molecule is not yet measured.

The two-photon transition benchmark which was started and abandoned earlier is now
worth running: **its outer loop was the gap**, the `tpa_quad` build, and that is what
made the earlier attempt meaningless rather than merely slow.

### Where the response path now stands

Every `_comp_nlr_fock` call in `src/pymodule` is handed the factors of its densities
-- checked by walking the call sites, not by memory. Ten nonlinear drivers and every
linear solver.

What is still outside it, unchanged and not a temporary omission in any of these
cases:

| | |
| --- | --- |
| unrestricted references | not covered anywhere in response |
| range-separated functionals | **now covered; see the section on them below.** This row read "asserted against" until the attenuated B vectors were added |
| more than one rank | refused, so the multi-rank half-size path never runs with it on |
| core excitations, restricted subspaces | fall back to the dense route, excluded in the solver's guard |
| complex trial vectors | fall back; the guard is a no-op today, since damped response carries them as real blocks |

## Two-photon transitions, the calculation which was not worth measuring before

This is the benchmark which was started and abandoned, on the grounds that the
resolution of the identity was not reaching the outer loop. It was not: the
`tpa_quad` build, the two-time perturbed Fock matrices the driver forms once its
solvers have finished, went through the four-centre integrals in both columns. The
solvers were fast and the thing they fed was not, so the ratio measured a fraction of
the calculation and called it the calculation. It is wired now, and this is the run.

Caffeine, five excited states, one rank of 14 threads, `def2-universal-jkfit` for all
three orbital bases.

| | basis | nao | four-centre | RI-JK simd | speedup | SCF |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Hartree-Fock | def2-svp | 246 | 154.08 s | 15.68 s | **9.83** | 12.41 -> 1.12 |
| | def2-svpd | 366 | 713.21 s | 34.74 s | **20.53** | 50.38 -> 2.42 |
| | def2-tzvp | 494 | 2362.74 s | 65.06 s | **36.32** | 173.57 -> 4.45 |
| B3LYP | def2-svp | 246 | 150.87 s | 41.38 s | **3.65** | 14.89 -> 4.04 |
| | def2-svpd | 366 | 626.89 s | 101.27 s | **6.19** | 56.22 -> 8.66 |
| | def2-tzvp | 494 | 1946.67 s | 150.07 s | **12.97** | 180.24 -> 13.33 |

**Forty minutes becomes a minute** at Hartree-Fock in def2-tzvp, and half an hour
becomes two and a half minutes at B3LYP. The whole grid took 1 h 50 min where the
four-centre half alone was 1 h 40 of it.

### The exponent, and what it is an exponent in

Fitted over the three bases of each row, cost proportional to nao to the p:

| | four-centre | RI-JK simd |
| --- | ---: | ---: |
| Hartree-Fock | 3.91 | **2.04** |
| B3LYP | 3.66 | **1.87** |

The ground state underneath them fits 3.77 and 1.97 at Hartree-Fock, 3.56 and 1.72 at
B3LYP -- the excited state part scales the way the ground state does, which is what
should happen when both are made of the same Fock builds.

**The second column is not the scaling of the method.** One fitting set serves all
three orbital bases here, 1242 functions throughout, so the auxiliary dimension is
held fixed while the orbital one grows and the exponent is in the orbital dimension
alone. The same caveat was written into the gradient tables and it has not changed.

### B3LYP gives up two thirds of the ratio, again

| basis | HF | B3LYP | HF / B3LYP |
| --- | ---: | ---: | ---: |
| def2-svp | 9.83 | 3.65 | 2.69 |
| def2-svpd | 20.53 | 6.19 | 3.32 |
| def2-tzvp | 36.32 | 12.97 | 2.80 |

Nothing new in it: what is left after the two-electron part goes away is the
quadrature, and the cubic response profile measured that directly -- 58 per cent of a
B3LYP calculation was integrating the functional. The ratio here is consistent with
that and does not independently establish it.

### The ground state is now a rounding error

Of the whole RI-JK run, ground state and excited state together, the ground state is
**6.4 to 6.7 per cent** at Hartree-Fock and **7.9 to 8.9** at B3LYP, and the fraction
does not grow with the basis. Four-centre it is a fifth to a quarter. There is nothing
left to win there.

### What the two columns agree on

The comparison is over everything the driver reports for all five states -- circular
and linear two-photon strengths, oscillator strengths, photon energies -- and the
column in the table is the worst of them, so a quantity which agrees cannot hide one
which does not.

| | def2-svp | def2-svpd | def2-tzvp |
| --- | ---: | ---: | ---: |
| Hartree-Fock | 1.9e-03 | 1.3e-03 | 8.9e-04 |
| B3LYP | 6.7e-04 | 1.2e-03 | 2.1e-03 |

It is a strength which sets the worst figure in five rows of six and an oscillator
strength in the other two; **the photon energies agree an order of magnitude better
throughout**, 8.99e-05 down to 2.84e-05. That ordering is the expected one -- an
energy is an eigenvalue of the linear problem and a strength is a product of response
vectors and a quadratic Fock matrix, so the fitting error enters it more times.

The two rows move in opposite directions with basis size, Hartree-Fock improving and
B3LYP worsening. Six points across two functionals is not enough to call that a
trend, and no explanation is offered here.

### What was fixed to get this table out

Two defects in the benchmark suite, neither of them in the physics, both found by
this run and not by reading:

- **The record asked for a cross section.** `TpaTransitionDriver` reports strengths;
  a cross section belongs to the full two-photon driver. Every row of the first run
  stored `None` and the table printed an empty column. The strengths were being
  saved all along, so the table above was rendered from that same data with nothing
  recomputed.
- **The output path was rebuilt every iteration**, so it followed the calendar. This
  run crossed midnight, the last two rows went to a second file, and the first was
  left as a stale prefix of itself. The path is formed once now.

## The open shell, where the driver had been closed shell only

The simd RI-JK driver built one Fock matrix from one density and one set of
orbitals. An unrestricted calculation asking for it did not fall back and did not
refuse: `_prepare_for_ri_fock_build` swapped in the simd driver, the open shell
branch then called a method only the conventional driver has, and the calculation
died part way with `AttributeError: 'SimdRIJKFockDriver' object has no attribute
'compute_screened_j_fock'`. A restricted open shell run did the same, sharing the
path. The **conventional** RI-JK driver had served both all along.

The closed shell assumption was two lines, not a design:

```cpp
auto fock = _drv.compute_fock_matrix(_bq_vectors, _basis, _aux_basis, density);
fock.scale(2.0);      // the density is one spin's, so the Coulomb enters twice
```

and one set of coefficients giving one exchange. What an open shell wants is the
Coulomb of the **total** density undoubled, and an exchange from each spin's own
occupied orbitals.

### What was added

One routine, taking the total density and a set of coefficients for each spin and
returning both matrices:

    compute(density, coefficients_alpha, coefficients_beta, exchange_scaling_factor)
        -> (F_alpha, F_beta)

Both spins run inside **one** pass over the ranges of the auxiliary basis, so a
range's B vectors are read once and serve both rather than being swept twice. Two
things follow from the spins occupying different numbers of orbitals: the W matrices
need two storages rather than one reused, which would otherwise be formed again at
every range of every build; and a function of a range costs the two spins together,
so the same memory buys about half the range. That halving is what holding two
spins' W matrices costs and is not a penalty of doing them in one pass.

### What it agrees with

def2-svp with the universal jkfit throughout, converged to 1e-8. The energy of the
four centre path, then what each of the two resolutions of the identity differs from
it by, and then -- the column that matters -- what the two of them differ from
**each other** by:

| | nao | naux | alpha/beta | four centre | conv vs 4c | simd vs 4c | **simd vs conv** |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| CH3 doublet, UHF | 29 | 129 | 5/4 | -39.5329552732 | 4.63e-06 | 4.63e-06 | **9.24e-14** |
| CH3 doublet, UB3LYP | 29 | 129 | 5/4 | -39.8091355099 | 1.90e-05 | 1.90e-05 | **8.53e-14** |
| O2 triplet, UHF | 28 | 154 | 9/7 | -149.4904009308 | 1.97e-04 | 1.97e-04 | **5.40e-13** |
| H doublet, UHF | 5 | 18 | 1/0 | -0.4992784057 | 0.00e+00 | 1.67e-16 | **1.67e-16** |
| CH3 doublet, ROHF | 29 | 129 | 5/4 | -39.5288782290 | 4.57e-06 | 4.57e-06 | **5.68e-14** |

The middle two columns being equal row by row is the signature to look for: the two
resolutions of the identity are the same approximation and should miss the four
centre answer by the same fitting error, which is 4.6e-06 here and 2.0e-04 on O2 --
a property of the fitting set and the molecule, not of either implementation. The
last column is the implementation, and at **5.4e-13 and below** it is convergence
noise.

Three of the rows are there for branches rather than for chemistry. **H has one
alpha orbital and no beta**, which is the spin with nothing to transform; **O2 has
spins differing by two**, where a single W storage reused between them would be
formed again at every range; and **ROHF** shares the open shell path and would have
been missed by testing the unrestricted driver alone.

### What is not here

**No timing.** Twenty-nine basis functions is a correctness size and nothing above
is a benchmark -- there is no speedup column on purpose, and what an open shell
build costs against the four centre way on a real molecule is not yet measured. The
structural expectation, which is a statement about the code and not a measurement,
is that an open shell build is two exchanges where a closed shell one is a single
exchange doubled, over B vectors formed once for both.

**The direct way is not served.** It accumulates the right hand side of its fitting
from the integrals during the same sweep which builds the exchange, and splits a
build into three calls so a rank can gather that fitting between the first and the
last. Two spins there means two exchanges and one fitting summed over both inside
that sweep, which those three calls have no shape for. It refuses with a sentence
saying so, which is also what removed the `AttributeError` above.

**Range separated functionals** remain refused for RI-JK, open shell included, by
the assertion which was already there.

## The gradient of an open shell, and the largest ratio in this file

The gradient driver was closed shell in the same way the Fock driver had been: one
density, one set of occupied orbitals, one exchange. Unlike the Fock case it
refused cleanly rather than dying in a binding, so an unrestricted calculation
stopped with a sentence -- it simply had no gradient to give.

### One routine, and where the factors come from

A separate routine and not an overload, because the two differ in **which density
the Coulomb half is of** and a caller passing the wrong one would get a gradient
which is merely wrong rather than an error:

    compute_open_shell(molecule, basis, aux_basis, bq_vectors, metric,
                       density, coefficients_alpha, coefficients_beta, a_x, ...)

The three-centre term underneath is not a second copy of its two hundred line
kernel loop. It takes a list of spins, one for a closed shell and two for an open
one, each carrying the orbitals it occupies and the fitted densities formed from
them. The factors follow from writing the closed shell expressions in spin summed
form, the closed shell ones being those expressions specialised to two identical
spins:

| | closed shell | open shell |
| --- | --- | --- |
| three-centre Coulomb | 4 c_P D(one spin) | c_P D(total) |
| three-centre exchange | -2 a_x, one spin | -a_x, each spin |
| two-centre Coulomb | 2 c_P c_Q | (1/2) c_P c_Q |
| two-centre exchange | -a_x, one Gram | -(a_x/2), each spin's Gram |

### Two defects the checks caught, both of which give a wrong gradient silently

**The multiply overwrites.** `_multiply` passes `beta = 0.0` to the library, so the
second spin's exchange was replacing the first rather than adding to it. The first
check below failed at 9.8e-02 relative, and the failure scaled with the fraction of
exact exchange, which placed it in the exchange half within one run. An accumulating
form was added beside the one the Gram product already used.

**A spin can occupy nothing.** The hydrogen atom is one alpha orbital and no beta.
The transformation asked the library for a product of no columns and it refused.
Guarded at the source rather than at each of the three places which consume it.

### What it agrees with

Three checks, in increasing independence. def2-svp throughout.

**Identical spins must return the closed shell answer.** Feed the two spins the same
orbitals and half the density each; the open shell routine has to reproduce the
restricted gradient exactly. This pins every factor in the table above with no
second implementation involved, and is the check which caught the multiply.

| | agreement |
| --- | ---: |
| Hartree-Fock, a_x = 1 | **4.4e-16** |
| B3LYP, a_x = 0.2 | **1.3e-15** |

**Against the four centre gradient, and against finite differences.**

| | four centre vs RI-JK | RI-JK vs numerical |
| --- | ---: | ---: |
| CH3 doublet, UHF | 9.047e-06 | **2.208e-07** |
| CH3 doublet, UB3LYP | 1.577e-05 | 1.710e-05 |
| O2 triplet, UHF | 1.214e-04 | **5.716e-07** |
| O2 triplet, UB3LYP | 3.750e-05 | 1.906e-05 |
| H doublet, UHF | 1.759e-34 | 5.551e-14 |
| H doublet, UB3LYP | 3.556e-17 | 1.110e-13 |

The right column is the one which tests this code: the analytic gradient against
finite differences of **its own** energy, to 2e-07 at Hartree-Fock. The left column
is the fitting error and is a property of the fitting set. The B3LYP rows of the
right column are larger because a numerical gradient of a functional carries the
quadrature grid's sensitivity to displacement; their Hartree-Fock counterparts at
the same geometry are two orders better, which is what says the difference is the
grid and not the term.

Three of the rows are branches rather than chemistry: **H has a spin which occupies
nothing**, **O2 has spins differing by two**, and both are in because the arithmetic
of a second spin is where this could go wrong quietly.

### What it bought

The caffeine cation, doublet, against the neutral's table row for row. One rank of
14 threads, `def2-universal-jkfit` throughout, the gradient timed twice and the
better kept.

| | four centre | RI-JK simd | speedup | the neutral's |
| --- | ---: | ---: | ---: | ---: |
| HF def2-svp | 21.17 | 0.84 | 25.3 | 11.9 |
| HF def2-svpd | 94.75 | 1.59 | 59.6 | 31.0 |
| HF def2-tzvp | 333.08 | 2.89 | 115.3 | 54.2 |
| HF def2-tzvpd | 834.48 | **4.41** | **189.2** | 91.8 |
| B3LYP def2-svp | 22.19 | 1.88 | 11.8 | 6.4 |
| B3LYP def2-svpd | 96.79 | 3.88 | 25.0 | 14.6 |
| B3LYP def2-tzvp | 335.82 | 5.84 | 57.5 | 34.5 |
| B3LYP def2-tzvpd | 842.87 | **9.37** | **89.9** | 47.8 |

**Fourteen minutes becomes four and a half seconds**, and 189 is the largest ratio
anywhere in this file.

### Why the open shell is *better* served than the closed one

Every row is about twice its neutral counterpart, which is the opposite of what the
self consistent field did, where the cation's extra iterations diluted the ratio.
The gradient has no iterations to dilute it, and the two ways pay differently for
the second spin:

| what the second spin costs | def2-svp | def2-svpd | def2-tzvp | def2-tzvpd |
| --- | ---: | ---: | ---: | ---: |
| four centre, HF | 2.39 | 2.25 | 2.44 | 2.37 |
| RI-JK simd, HF | **1.14** | **1.17** | **1.15** | **1.15** |
| four centre, B3LYP | 2.42 | 2.26 | 2.43 | 2.38 |
| RI-JK simd, B3LYP | **1.50** | **1.53** | **1.46** | **1.41** |

The four centre way builds a second exchange from scratch and pays about 2.4 for it.
The resolution of the identity forms the B vectors once, for both spins, and the
second spin costs only the contraction against them -- 15 per cent at Hartree-Fock.
**The dearest thing the method forms does not depend on spin**, which is the whole
of why the ratio doubles.

B3LYP pays more for the second spin than Hartree-Fock does, 1.46 against 1.15, and
the reason is the quadrature: it is spin resolved, so the functional's integration
genuinely doubles where the exchange contraction is the only part which grows in the
Hartree-Fock rows.

### The exponents do not move

Fitted over the four bases, cost proportional to nao to the p:

| | cation | neutral |
| --- | ---: | ---: |
| four centre, HF | 4.02 | 4.01 |
| RI-JK simd, HF | 1.82 | 1.82 |
| four centre, B3LYP | 3.98 | 3.97 |
| RI-JK simd, B3LYP | 1.71 | 1.78 |

An open shell is a constant times the work and not a different scaling, which is
what the arithmetic says it should be and is worth having measured rather than
assumed. The same caveat as every other table here: one fitting set serves all four
orbital bases, so the second column of each pair is an exponent in the orbital
dimension alone and not the scaling of the method.

### Unchanged

The l = 4 ceiling on the orbital centres is in the 175 generated derivative
kernels and not in this driver. The single rank restriction is the same coupling of
the fitting across the whole auxiliary basis. Range separated functionals are
refused, open shell included. The conventional resolution of the identity has no
open shell gradient at all, so `ri_jk` without `ri_jk_simd` is refused rather than
measured.

## Optimizing a radical, where the two ways agree on the path

The open shell gradient made an unrestricted geometry optimization work without a
line being written for it: the optimizer builds its gradient driver from whatever
self consistent field driver it was handed, so the gradient was the only thing
missing. What was checked rather than assumed is that the arrangements made for the
closed shell case still hold -- the banner is said once over an optimization and not
once a step, and a mode which cannot differentiate is refused before the first step
rather than after it.

The molecule is a nitronyl nitroxide core radical, C7H13N2O2, 24 atoms, a doublet
with 43 alpha and 42 beta electrons. One rank of 14 threads,
`def2-universal-jkfit` throughout, run to convergence with no cap on the steps.

| | nao | four centre | RI-JK simd | speedup | steps |
| --- | ---: | ---: | ---: | ---: | --- |
| HF def2-svp | 219 | 347.0 s | **15.2 s** | **22.8** | 8 and 8 |
| B3LYP def2-svp | 219 | 411.2 s | **60.2 s** | **6.8** | 9 and 9 |
| HF def2-svpd | 330 | 1198.6 s | **26.9 s** | **44.5** | 7 and 7 |
| B3LYP def2-svpd | 330 | 1521.9 s | **135.7 s** | **11.2** | 9 and 9 |

The optimized energies agree to 3.9e-04 at Hartree-Fock and 8e-05 at B3LYP, which is
the fitting error and not a difference in where the two paths stopped.

### The cleanest ratio in this file, and why

**Both ways take the same number of steps in all four rows.** That makes the whole
run ratio and the per step ratio the same number, to two decimals. Nothing here is
a total divided by a total over different amounts of work.

It is worth saying why that matters, because the closed shell table does not have
it. There, caffeine at def2-svp took 33 steps by the four centre way and 27 by the
resolution of the identity, so its 6.87 is a ratio of two different journeys and the
per step figure is 5.62. A speedup which moves when the optimizer takes a different
path is a weaker measurement than one which does not, and this suite happens to give
the stronger kind.

### Against the closed shell table

| | nitroxide, a radical | caffeine, closed shell |
| --- | ---: | ---: |
| HF def2-svp | **22.8** | 6.87 |
| B3LYP def2-svp | **6.8** | 3.39 |
| HF def2-svpd | **44.5** | 17.01 |
| B3LYP def2-svpd | **11.2** | 7.67 |

Two to three times better for the open shell, the same way round as the gradient
table and for the same reason: the four centre way builds a second exchange from
nothing while the resolution of the identity forms its B vectors once for both
spins. The molecules differ -- 219 and 330 functions against 246 and 366 -- so this
is not a controlled comparison of the two spin cases, and the gradient section,
where the same molecule is measured both ways, is the one which establishes the
effect. This table is consistent with it.

### The estimate, and two defects in the suite

**The first estimate was 1 to 1.5 hours and the run took 15 minutes.** The cause is
worth recording because the information was in the file: the geometry's own comment
line says "Optimized", so it converges in eight steps where the caffeine table's
thirty was the number used to predict it. A step count taken from a different
molecule's table is not an estimate.

**The runner overwrote its own record.** The output path carried the date, the
machine, the molecule and the spin state but not the basis, so running def2-svpd
after def2-svp on the same day replaced the first file rather than writing beside
it. The def2-svp numbers above survive because they had been committed. The path
now carries the bases, and a run which would overwrite an existing record refuses
to start instead: one file per run is what makes every ratio inside a file
comparable, and a record replaced in silence is worse than one appended to.

The optimizer also drops its checkpoints beside whatever ran it, which nearly went
into a commit. They are ignored now.

## The size of a block, when the buffer doubles

The range separated three-center driver inherited its three sizing constants from
the unattenuated one without measurement: a budget of 256 MB for the buffer a
thread holds, at most 256 atom pairs to a block and at least 8. That inheritance is
not obviously safe, because the buffer of a combination here is the larger of the
two -- it carries two chains of Boys values, one for each operator -- so the same
number of atom pairs costs twice the working set. The question is whether 256 still
sits inside the flat part of the curve at that size.

The table of rows is **1.86 to 1.99 times** the unattenuated one across its whole
range, which is the doubling and the shared prefactors. What that does to the block
size the budget computes:

| | l bra / aux | rows plain | rows rs | block plain | block rs | MB a thread, rs |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| caffeine, def2-TZVP + jfit | 3 / 4 | 15103 | 29570 | 256 | 256 | 60.6 |
| tagrisso, def2-SVP + jkfit | 2 / 4 | 4997 | 9718 | 256 | 256 | 19.9 |
| c60, cc-pVDZ + RIFIT | 2 / 3 | 3012 | 5808 | 256 | 256 | 11.9 |
| `(gg\|i)` | 4 / 6 | 80964 | 160167 | 256 | **209** | 267.8 |
| `(ii\|l)` | 6 / 8 | 580552 | 1154910 | 57 | **29** | 267.9 |

**The budget is not binding for any basis these drivers are used with.** It first
bites at `(gg\|i)`, and everything to `(ff\|g)` reaches the ceiling of 256 in both
drivers long before the 256 MB is spent. So the only constant the doubling can
reach is the ceiling, and the ceiling is what was swept.

### The doubled buffer does not move the optimum

Milliseconds, fourteen threads, best of three, threshold 1e-12, omega 0.3. Each
case runs in its own process. The two columns of a case were taken in one session,
so they are an A/B and not two logs read against each other.

| atom pairs | caffeine plain | caffeine rs | tagrisso plain | tagrisso rs | c60 plain | c60 rs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 130.2 | 260.2 | 662.6 | 1370.8 | 1633.9 | 3453.5 |
| 32 | 86.4 | 175.8 | 362.5 | 747.7 | 941.1 | 2001.4 |
| 128 | 84.1 | 169.1 | 270.3 | 559.3 | 700.3 | 1434.5 |
| 256 | 84.0 | **168.3** | 250.1 | **531.2** | 687.0 | **1403.7** |
| 512 | 84.8 | 167.7 | 248.9 | 528.4 | 682.3 | 1419.5 |

The range separated curve has the same shape as the plain one on all three cases --
steep below 128, flat from 256 -- so the inherited ceiling sits in the flat range
for both and needs no change.

The ratio is the sharper reading, because it is free of everything the two drivers
share:

| block size | 8 | 32 | 128 | 256 | 512 |
| --- | ---: | ---: | ---: | ---: | ---: |
| caffeine, def2-TZVP + jfit | 2.00 | 2.04 | 2.01 | 2.00 | 1.98 |
| tagrisso, def2-SVP + jkfit | 2.07 | 2.06 | 2.07 | 2.12 | 2.12 |
| c60, cc-pVDZ + RIFIT | 2.11 | 2.13 | 2.05 | 2.04 | 2.08 |

**It does not drift with the block size.** The second operator costs what its second
Boys chain costs, at every size, and nothing further through the working set. Were
60 MB a thread hurting the caches where 30 did not, this ratio would grow towards
the large sizes; it is flat to within the scatter of a best of three.

### Caffeine cannot answer the question, and c60 can

Caffeine's block count is 40 at 128, at 256 and at 512 alike. It has four distinct
atom bases, so ten pairs of them, and the atom pairs of its groups run out before
the ceiling does -- its flat tail is a property of the molecule and not evidence
about the ceiling. Tagrisso goes 100, 64, 48 blocks over those three sizes and c60
goes 14, 7, 4, so those two are the cases which actually probe it, and both are flat
there too. A sweep on caffeine alone would have concluded nothing while appearing
to.

### How it was measured, and why the plain column was re-measured

The ceiling is a `static constexpr`, so a sweep means a rebuild for each value. Both
drivers were patched to read it from the environment instead, which makes it one
rebuild and five free points; the patch was reverted afterwards and the library
rebuilt, and the revert was checked by confirming that the variable no longer moves
the block count.

The plain column above is **not** the one in "The size of a block" earlier in this
file, which reads 92.2, 310.8 and 837.6 at 256 against the 84.0, 250.1 and 687.0
here. The driver has changed since that table was taken and the two are not
comparable, which is the whole reason for measuring plain again beside the range
separated driver rather than reading the new numbers against the old ones.

## The range separated hybrids, where the split costs one way twice and the other a quarter

A hybrid range separated functional splits its exchange between the plain operator
and the attenuated one, and the two ways of building pay for that split very
differently. The four-centre way makes a second full sweep of its own kernels on
every iteration, `kx_rs` on top of `2jkx`. The simd resolution of the identity holds
a second set of B vectors, formed once, and adds a second exchange inside the same
pass over the auxiliary basis it was already making.

The conventional RI-JK driver is not a column here. It has no attenuated B vectors
and refuses a range separated functional, so it would be a column of refusals.

`OMP_NUM_THREADS=14`, one rank, `def2-universal-jkfit` throughout, convergence 1e-8.
The records are `benchmarks/data/scf/2026-09-18_m4max_caffeine_rs_closed.json` and
`..._nitroxide_rs_m2.json`, and the runner is `benchmarks/scripts/scf_rs_laptop.py`.
The provenance says the tree was dirty: the only thing uncommitted was that runner,
which was written for this measurement and is committed with it.

### Caffeine, closed shell

**plain** is the speedup of B3LYP on the same molecule and basis, from the table
earlier in this file, and is what the range separated column is to be read against.

| functional | basis | nao | four-centre | RI-JK simd | speedup | build only | plain |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CAM-B3LYP | def2-svp | 246 | 27.54 | 5.20 | 5.30 | 18.0 | 3.75 |
| | def2-svpd | 366 | 107.72 | 11.63 | 9.26 | 35.1 | 6.48 |
| | def2-tzvp | 494 | 361.56 | 16.95 | 21.34 | 68.4 | 13.47 |
| | def2-tzvpd | 614 | 868.74 | 28.31 | **30.69** | 99.0 | 19.20 |
| wB97X-D4 | def2-svp | 246 | 27.18 | 5.12 | 5.31 | 17.7 | 3.75 |
| | def2-svpd | 366 | 106.92 | 11.46 | 9.33 | 34.9 | 6.48 |
| | def2-tzvp | 494 | 360.31 | 17.54 | 20.54 | 64.9 | 13.47 |
| | def2-tzvpd | 614 | 868.04 | 29.23 | **29.70** | 95.0 | 19.20 |

### Nitroxide, a doublet radical, unrestricted

C7H13N2O2, 24 atoms, 43 alpha and 42 beta electrons. There is no plain column: the
nitroxide tables elsewhere in this file are optimizations and not self consistent
field runs, and a speedup taken from a different kind of calculation is not a
comparison.

| functional | basis | nao | four-centre | RI-JK simd | speedup | build only |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| CAM-B3LYP | def2-svp | 219 | 55.49 | 7.72 | 7.18 | 26.7 |
| | def2-svpd | 330 | 207.22 | 16.45 | 12.59 | 51.2 |
| | def2-tzvp | 419 | 647.43 | 23.37 | 27.71 | 99.8 |
| | def2-tzvpd | 530 | 1464.14 | 39.65 | **36.93** | 135.5 |
| wB97X-D4 | def2-svp | 219 | 55.48 | 7.67 | 7.23 | 26.9 |
| | def2-svpd | 330 | 197.45 | 17.22 | 11.47 | 46.6 |
| | def2-tzvp | 419 | 621.66 | 24.38 | 25.50 | 91.7 |
| | def2-tzvpd | 530 | 1403.55 | 41.06 | **34.18** | 124.8 |

### The split is served better than the plain hybrid, by half again

Caffeine's range separated rows run **1.4 to 1.6 times** the speedup of B3LYP on the
same molecule and the same basis: 5.30 against 3.75, 9.26 against 6.48, 21.34 against
13.47, 30.69 against 19.20. That is the whole point of the arrangement and it is
worth saying why it is not a paradox that a harder Fock matrix is built relatively
faster.

The obvious explanation is wrong and worth writing down, because it is the one
anybody would reach for: it is **not** that the attenuated operator is cheap for the
resolution of the identity and dear for the four-centre way. The two-electron time
roughly doubles on both sides. At def2-tzvpd on caffeine, going from B3LYP to
CAM-B3LYP takes the four-centre build from 410.2 seconds to 847.8, a factor of 2.07,
and the simd build from 4.8 to 8.6, a factor of 1.78. Those are close enough to each
other that they cannot be where a half again comes from.

Where it comes from is **what fraction of the run the build is**. The two-electron
part is 95 to 98 per cent of a four-centre run at these bases and 17 to 30 per cent
of a simd one:

| basis | 2e share, four-centre | 2e share, simd |
| --- | ---: | ---: |
| def2-svp | 74.9% -> 85.1% | 17.9% -> 25.1% |
| def2-svpd | 86.0% -> 92.1% | 17.0% -> 24.3% |
| def2-tzvp | 93.2% -> 96.4% | 21.9% -> 30.1% |
| def2-tzvpd | 95.4% -> 97.6% | 21.4% -> 30.3% |

Doubling something which is 95 per cent of the wall doubles the wall; doubling
something which is a fifth of it adds a quarter. The four-centre run goes 430.1 to
868.7 seconds, 2.02 times, and the simd run 22.4 to 28.3, 1.26 times. **That is the
whole of the effect**, and it says the gain belongs to the resolution of the identity
having made the build small in the first place rather than to anything the attenuated
path does especially well.

There is a second, smaller effect underneath it which is the attenuated path's own:
1.78 against 2.07 on the builds, so the ratio of the builds alone improves from 85.4
for B3LYP to 99.0 for CAM-B3LYP at def2-tzvpd. Forming both operators over one set of
primitive pairs and transforming B vectors which are already resident is worth about
fifteen per cent of the build. It is real and it is not what the table above is
mostly showing.

### The two functionals are one measurement

Every pair of rows agrees within a few per cent although the functionals do not:
CAM-B3LYP is 0.190 of exact exchange at short range and 0.650 at long, wB97X-D4 is
0.167 and 1.000. The cost follows the number of passes and not the coefficients,
which is what both implementations say it should, and a pair which disagreed would
have meant one of them was doing arithmetic that depended on the numbers.

### Every ratio here is like for like

All thirty two runs converged, in 20 to 23 iterations. No speedup in either table is
a short journey divided by a long one, which the closed shell optimization table
earlier in this file cannot say of all of its rows.

The radical was the better behaved of the two molecules, which was not the
expectation: a doublet with long-range corrected functionals looked like the
convergence risk of the pair, and the tagrisso cation and triplet had failed to
converge in a hundred iterations by every path. Nitroxide took 21 or 22 everywhere,
against caffeine's 20 to 22 and the caffeine cation's 30 to 45.

### What now bounds the simd path is the quadrature

The exchange-correlation grid is **54 to 62 per cent** of every simd run on caffeine
and **65 to 70 per cent** on nitroxide. At def2-tzvpd on nitroxide it is 26.0 seconds
of a 39.7 second run against 10.5 seconds of two-electron work. The Fock build is no
longer the thing to work on for these calculations; the grid is.

### The estimate, and the third time the analogy lost to the first measurement

Quoted at 3 to 3.5 hours, re-estimated at 2 hours 36 after the first measured pair,
and it took **2 hours 8**. The first figure was about 60 per cent high because the
cost of the second four-centre pass was guessed at 2 to 2.5 times a plain hybrid and
is 1.84 at def2-svp, rising to 2.02 at def2-tzvpd as the fixed costs shrink beside
the integrals. Re-estimating from one measured pair was within 10 per cent per row,
which is the third time in these notes that the analogy was wrong and the first
measurement was right.

## The range separated hybrids in TDA and TD-DFT

The excited state solvers reach the resolution of the identity through the factorised
Fock build, which takes a batch of trial vectors as the factors they were made from
and transforms the shared left factor once for all of them. A hybrid range separated
functional adds a second set of B vectors there as it does in the ground state, and
its exchange is accumulated into the same matrices inside the same pass over the
auxiliary basis.

Five states, caffeine, `def2-universal-jkfit`, `OMP_NUM_THREADS=14`, one rank. The
records are `benchmarks/data/tda/2026-09-18_m4max_caffeine_rs.json` and
`..._rs_rpa.json`, and the runner is `benchmarks/scripts/tda_laptop.py`, which now
takes the functionals as an argument and writes a range separated grid to a file of
its own. **B3LYP is measured in the same file and the same run as the rows it is the
control for**, rather than being read out of the plain tables earlier in this file.

**This is def2-svp alone.** The four basis grid the plain tables carry was estimated
at eight hours once this row had been measured, and was not run. What is here
establishes the ratio; how it grows with the basis is not measured for these
functionals.

Seconds for the excited state part, with the ground state beneath it.

| solver | functional | four-centre | RI-JK simd | speedup | vs control | SCF |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| TDA | B3LYP | 51.12 (10 it) | 12.18 (10 it) | 4.20 | control | 15.3 / 4.1 |
| | CAM-B3LYP | 102.32 (12 it) | 17.07 (12 it) | **5.99** | 1.43 | 27.6 / 5.3 |
| | wB97X-D4 | 108.94 (13 it) | 18.12 (13 it) | **6.01** | 1.43 | 27.3 / 5.2 |
| TD-DFT | B3LYP | 49.65 (11 it) | 13.63 (11 it) | 3.64 | control | 15.2 / 4.1 |
| | CAM-B3LYP | 99.30 (14 it) | 20.55 (14 it) | **4.83** | 1.33 | 27.3 / 5.2 |
| | wB97X-D4 | 105.51 (14 it) | 21.66 (14 it) | **4.87** | 1.34 | 27.0 / 5.2 |

**Every row took the same number of iterations both ways**, so each speedup is a
ratio of two equal amounts of work and not of two different journeys.

### What the split costs each way

Per iteration, against the B3LYP control of the same solver:

| | four-centre | RI-JK simd |
| --- | ---: | ---: |
| TDA, CAM-B3LYP | 1.67x | 1.17x |
| TDA, wB97X-D4 | 1.64x | 1.14x |
| TD-DFT, CAM-B3LYP | 1.57x | 1.18x |
| TD-DFT, wB97X-D4 | 1.67x | 1.25x |

The same shape as the ground state and for the same reason set out there: the two
ways pay similar multiples on the build itself, and what differs is how much of the
run the build is. It is a smaller effect here -- 1.43 and 1.33 against the 1.4 to 1.6
of the ground state -- because an excited state calculation carries more that is
neither Coulomb nor exchange.

### The iterations are not the same, and that is the functional's doing

The range separated runs take **12 and 13 iterations where B3LYP takes 10** in the
Tamm-Dancoff approximation, and 14 against 11 in linear response. Neither path causes
it: the four-centre and the simd rows of a given functional agree on the count
exactly. It is worth recording because it is what made the estimate of the four basis
grid grow: a twenty to twenty seven per cent longer journey multiplies whatever the
per iteration cost is, and an estimate built from per iteration multipliers alone is
that much low.

### These rows were measured twice, and the first set was thrown away

The first measurement of this table ran while a code generator held one core at a
hundred per cent for its whole duration, so it had thirteen of the fourteen threads
it asked for. The rows above are the second measurement, taken after that finished,
with the machine otherwise idle.

It is worth recording what that cost and what it did not, because the numbers happen
to say it cleanly. The contended four-centre rows were **11 to 13 per cent** slow and
the contended simd rows **4 to 5**, which is the thread bound path losing more from
losing a thread. So the **speedups were inflated**: 6.41 where it is 5.99 for
CAM-B3LYP in the Tamm-Dancoff approximation. The **ratio against the control was
not**: 1.42 contended against 1.43 clean, because the control was measured under the
same load, minutes apart.

The first reading of this was wrong in a way worth naming. The contended B3LYP row
came out at 57.88 seconds against the 51.64 of the plain table taken two days before,
and that was written down here as the machine drifting between runs. It was not
drift: measured again on a quiet machine the same row is **51.12**, within one per
cent of the two day old number. There was nothing to drift. **A benchmark which
disagrees with an older one by ten per cent is a reason to look at what else is
running, not a reason to write a sentence about drift.**

## What the range separated path now covers, and how each part was checked

The attenuated operator reaches every Fock build the resolution of the identity
serves, ground state and response. This section is the record of what that is and of
what established it, so a reader does not have to infer coverage from which sections
happen to carry timings.

Everything below is the restricted reference. The unrestricted reference is covered
in the ground state and nowhere in response, for the reason the table above gives:
the unrestricted response path has no factorised Fock build at all, so there is
nothing there to add an operator to.

| what | wired in | checked against | agreement |
| --- | --- | --- | ---: |
| SCF, closed shell | `scfdriver` | four-centre Fock matrices | 0.86 to 1.14 x control |
| | | converged energies | 0.99 to 1.31 x control |
| SCF, open shell | `scfdriver` | the same, two spins | as above, within six per cent |
| TDA | `linearsolver` | excitation energies, five states | 8.1e-06 to 9.2e-06 a.u. |
| TD-DFT | `linearsolver` | as above | 8.1e-06 to 9.2e-06 a.u. |
| linear response, CPP | `linearsolver` | polarizabilities at two frequencies | 1.00 x control |
| ten nonlinear drivers | `nonlinearsolver` | the results each driver returns | 0.98 to 1.12 x control |

**Every comparison is against the plain hybrid of the same molecule and basis**, run
in the same session, and not against a fixed tolerance. A fitting set carries an
error of its own and a fixed number would measure that rather than the code. The
question which has an answer is whether fitting the attenuated operator is as good as
fitting the plain one, and across every row of every one of those checks it is.

### The nonlinear drivers, one builder and sixteen modes

`nonlinearsolver._comp_two_el_int` is the only place the factorised build is chosen
for a nonlinear calculation, and the sixteen mode strings which reach it are covered
exactly by the two sets it dispatches on -- five cubic, eleven quadratic, none left
over. So wiring is a property of that one function and not of the drivers.

Checked driver by driver even so, water in def2-SVP, CAM-B3LYP against a B3LYP
control, four-centre against the simd path:

| driver | control | CAM-B3LYP |
| --- | ---: | ---: |
| quadratic response | 4.57e-04 | 4.65e-04 (1.02) |
| cubic response | 6.92e-05 | 7.05e-05 (1.02) |
| two-photon absorption, full | 1.09e-04 | 1.13e-04 (1.04) |
| two-photon absorption, reduced | 9.08e-05 | 9.32e-05 (1.03) |
| two-photon transitions | 3.03e-05 | 2.98e-05 (0.98) |
| three-photon transitions | 7.12e-05 | 7.20e-05 (1.01) |
| third harmonic generation | 1.24e-04 | 1.24e-04 (1.00) |
| third harmonic, reduced | 1.24e-04 | 1.24e-04 (1.00) |
| second harmonic generation | 1.49e-04 | 1.50e-04 (1.01) |
| second harmonic, reduced | 1.49e-04 | 1.50e-04 (1.01) |

### Three ways this check passed while proving nothing

Worth writing down, because each of them produced green rows.

**The settings were dropped in silence.** The resolution of the identity keys of a
nonlinear driver live in its **method** settings, where the linear solver keeps them
in its **response** settings. Passed in the wrong dictionary they are discarded
without a word, `ri_jk` comes back false, the auxiliary basis reverts to the default,
and both columns of the comparison run four centres. Six rows then agreed to
**2.4e-15** and reported success. The tell was that the agreement was far too good:
the fit carries an error near 1e-4 here, so machine precision is impossible if one
side used it. The check now asserts that the driver holds B vectors at the omega its
functional asks for, and that the two columns differ by more than 1e-10.

**A phase read as a catastrophe.** The transition drivers reported a relative
difference of exactly **2.00**, which is what a sign flip gives at the largest
element. A transition moment is defined up to the phase of an eigenvector and the two
ways of building can land on either. Nothing was wrong: the three-photon strengths
were 1688.54 and 1688.38. Those drivers are compared as magnitudes now; the others
keep signed comparison, their response functions having a sign which means something.

**A flattener which never reached the answer.** It walked the top level of a result
and stopped, so the strengths of a transition driver, nested two dictionaries deep,
were never compared at all -- while the excitation energies and dipoles beside them
were, and passed.

None of the three was a defect in the code under test. All three were defects in the
thing measuring it, and two of them presented as success.

### What is still outside

The **gradient** refuses a range separated functional and should keep refusing: there
are no attenuated derivative kernels, so the long-range term has no derivative on this
path and a run which proceeded would differentiate the wrong energy expression in
silence. That assert matters more now than it did, because the self consistent field
it would be differentiating converges.

The **unrestricted response path** and the **exciton driver** are untouched, the first
because it has no factorised build to extend and the second because it refuses range
separation for reasons of its own.

## The unrestricted excited states, where the open shell helps rather than hurts

The response path had no resolution of the identity for an unrestricted reference at
all -- not a range separated gap but an absent path, and one which did not refuse
either: a calculation which asked for RI-JK formed the B vectors, paid for them, and
built every Fock matrix from four-centre integrals. The factorised build now serves
both spins.

What it forms, per density, is the Coulomb of the two spins' densities added and
**undoubled**, less each spin's own exchange. The two spins share that Coulomb and
share nothing else: their left factors are the occupied orbitals of each of them and
differ both in what they are and in how many there are, so each exchange is formed
from its own transformation.

Nitroxide, the doublet radical of the optimization tables, def2-SVP,
`def2-universal-jkfit`, five states, `OMP_NUM_THREADS=14`, one rank. Twelve runs in
27 minutes. The records are
`benchmarks/data/tda/2026-09-18_m4max_nitroxide_rs_m2.json` and `..._rpa.json`, and
B3LYP is measured in the same file and the same run as the rows it is the control
for.

**restricted** is the same functional and solver on caffeine from the section above,
and is there for the comparison the next subsection makes. It is a different
molecule -- 219 functions against 246 -- so it is not a controlled comparison.

| solver | functional | four-centre | RI-JK simd | speedup | vs control | restricted |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| TDA | B3LYP | 151.40 (13 it) | 30.43 (13 it) | 4.98 | control | 4.20 |
| | CAM-B3LYP | 213.84 (11 it) | 31.04 (11 it) | **6.89** | 1.38 | 5.99 |
| | wB97X-D4 | 199.21 (10 it) | 28.82 (10 it) | **6.91** | 1.39 | 6.01 |
| TD-DFT | B3LYP | 127.08 (12 it) | 28.39 (12 it) | 4.48 | control | 3.64 |
| | CAM-B3LYP | 206.47 (12 it) | 35.36 (12 it) | **5.84** | 1.30 | 4.83 |
| | wB97X-D4 | 205.92 (12 it) | 35.27 (12 it) | **5.84** | 1.30 | 4.87 |

Every row took the same number of iterations both ways. The excitation energies of
the two paths agree to **7.1e-06 to 7.9e-06 hartree** over the five states.

### A prediction which was wrong, and why

Before this was measured it was written down that the unrestricted speedup would come
out **clearly lower** than the restricted one, because the two spins cannot share the
transformation of the left factor and an unrestricted batch therefore costs two where
a restricted one costs one.

**All six rows came out higher**, by fifteen to twenty per cent.

The reasoning was half an argument. The part which is true is that the resolution of
the identity pays two transformations instead of one. The part which was left out is
that the four-centre way pays worse: it forms the Coulomb and two exchanges where a
restricted build forms a single `2jkx` matrix, while the fitted path still shares one
Coulomb between the spins. The open shell is harder on the dense way than on the
fitted one, so the ratio rises. What was costed was one side of a ratio and what was
stated was a conclusion about the ratio.

The same effect is in the ground state tables, where the open shell speedups also run
ahead of the closed shell ones, and the gradient section states the mechanism
correctly -- the four-centre way builds a second exchange from nothing while the
resolution of the identity forms its B vectors once for both spins. It was not
carried across to the excited states, where the arithmetic is the same.

**How far it goes is not established here.** Nitroxide and caffeine are different
molecules, so the six pairs above agree on the direction and say nothing reliable
about the size. The measurement which would settle it is one molecule in both spin
states, as the gradient section did for the ground state, and it has not been made.

### What the split costs, unchanged from everywhere else

Per iteration, against the B3LYP control of the same solver: **1.63 to 1.67 times**
for the four-centre way and about **1.2** for the simd one, which are the same two
numbers the ground state and the restricted excited states gave. The two range
separated functionals agree with each other to within three per cent in every row,
here as everywhere: what the split costs follows the number of passes and not the
coefficients.

## The range separated gradient, where the split costs less than it does in the build

The gradient of a hybrid range separated functional carries a second exchange, and
both ways pay for it. The four-centre way makes a whole further pass of its own
derivative kernels, `kx_rs` on top of `2jkx`, once per atom. The resolution of the
identity contracts a second set of B vectors against a second derivative tensor which
the same kernel wrote, on one sparsity pattern and in one sweep of the auxiliary
basis.

Caffeine, `def2-universal-jkfit`, CAM-B3LYP, one rank of 14 threads on the M4 Max,
best of two, one SCF and then the gradients in one process per case, at `255140804`.
The records are `benchmarks/data/gradient/2026-09-18_m4max_caffeine_rs.json` and
`..._b3lyp_control.json`, rendered beside them as `.md` and `.pdf` by
`benchmarks/scripts/render_runs.py`, and the runner is
`benchmarks/scripts/grad_rs_laptop.py`. The provenance of both says the tree was
dirty: the only thing uncommitted was that runner, which was written for this
measurement and is committed with it.

| functional | basis | nao | four-centre | RI-JK simd | speedup | vs four-centre |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| CAM-B3LYP | def2-svp | 246 | 16.72 | 1.97 | 8.49 | 1.5e-05 |
| | def2-svpd | 366 | 77.01 | 4.07 | 18.92 | 1.7e-05 |
| | def2-tzvp | 494 | 257.77 | 6.90 | 37.36 | 6.1e-05 |
| | def2-tzvpd | 614 | 655.69 | 11.01 | **59.55** | 6.2e-05 |

### The control, and why one was run at all

**What the split costs is an absolute wall time divided by an absolute wall time**,
and that division only means something when both halves were measured in one sitting
on one tree. The B3LYP rows already in this file are from a different day and a
different tree, and dividing by them was how a 1.21 once became a 1.42. So B3LYP was
re-measured here, immediately after the range separated run exited, from the same
file:

| functional | basis | nao | four-centre | RI-JK simd | speedup | vs four-centre |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| B3LYP | def2-svp | 246 | 9.29 | 1.29 | 7.20 | 1.6e-05 |
| | def2-svpd | 366 | 43.67 | 2.68 | 16.29 | 1.7e-05 |
| | def2-tzvp | 494 | 139.76 | 4.29 | 32.58 | 6.2e-05 |
| | def2-tzvpd | 614 | 359.18 | 6.99 | 51.38 | 6.3e-05 |

It also answers a question nobody asked. Against the corrected B3LYP table earlier in
this file, the **four-centre column reproduces to between 1.2 and 2.1 per cent**,
which is the usual spread. The **RI-JK simd column is 3.2 to 7.0 per cent slower**
than it was. That is outside the spread of the four-centre column measured beside it,
so it is more likely real than noise: the range separated work restructured
`_compute_three_center` to loop over a vector of derivative tensors where it had one,
and the plain path now goes through the same general body. Seven per cent of four and
a quarter seconds is not worth undoing, but it should not be discovered later and
mistaken for a regression of something else.

### What the split costs each way

| basis | four-centre | RI-JK simd | speedup, CAM over B3LYP |
| --- | ---: | ---: | ---: |
| def2-svp | 1.80 | 1.53 | 1.18 |
| def2-svpd | 1.76 | 1.52 | 1.16 |
| def2-tzvp | 1.84 | 1.61 | 1.15 |
| def2-tzvpd | 1.83 | 1.58 | 1.16 |

**The gradient gets off more lightly than the build does.** In the SCF the four-centre
way pays 2.07 for the split and the fitted way about 1.2; here it is 1.80 to 1.84 and
1.52 to 1.61. The two columns move toward each other because the gradient is not only
integrals: quadrature, the one-electron terms and the fitted densities are in both
totals and none of them doubles, so the pass which does double is a smaller share of a
larger whole on each side.

What survives is the ratio of the ratios. **The range separated gradient is served 15
to 18 per cent better by the resolution of the identity than a plain hybrid is**,
which is the same fifteen to twenty per cent the ground state and the excited states
reported, arrived at by a different route.

### Splitting the exchange adds nothing to the fitting error

The last column of the two tables is the same column twice: 1.5e-05, 1.7e-05, 6.1e-05,
6.2e-05 for CAM-B3LYP against 1.6e-05, 1.7e-05, 6.2e-05, 6.3e-05 for B3LYP. The
attenuated half is fitted in a metric of its own, and **that metric is singular by
construction** -- the Fourier transform of the attenuated operator carries
`exp(-k^2/4 omega^2)`, so its smallest eigenvalues run to 1e-15. None of that reaches
the gradient, for the same reason it never reached the Fock matrix: the attenuated
integrals are zero in the directions the metric is blind in. This is the first time
that has been visible in a derivative rather than argued from the operator.

### The exponent does not move when the exchange splits

Fitted over the four bases by the renderer:

| functional | method | p |
| --- | --- | ---: |
| CAM-B3LYP | four-centre | 4.00 |
| CAM-B3LYP | RI-JK simd | 1.86 |
| B3LYP | four-centre | 3.97 |
| B3LYP | RI-JK simd | 1.81 |

Four oh oh against three nine seven, and one eight six against one eight one. **The
second exchange is a constant on both sides and not a power**, which is what the pass
structure says it should be: the four-centre way repeats a pass it already makes and
the fitted way contracts a second tensor over the same pattern. Neither adds an index.
The fitted exponent near 1.8 rather than 1 is the naux caution below, not the split.

### The caution the earlier table states applies here unchanged

One fitting set serves every row: naux is 1242 while nao runs 246 to 614. So the
speedup widening from 8.5 to 59.6 is partly the denominator being held still, and a
reader who takes 59.6 as a trend and extrapolates it will be wrong.

### What the gradient was checked against

The agreement column above is bounded by the fitting error and cannot settle a
coefficient: the two codes converge different densities, and a factor wrong by a few
per cent would sit under 1e-05. What settles it is differencing the SCF energy of the
**same** method, where nothing of the fitting cancels.
`benchmarks/scripts/rs_scf_gradient_check.py` does that for water and a hydroxyl
radical in def2-svp, B3LYP as a control:

| case | shell | functional | analytic vs four-centre | analytic vs finite difference | four-centre vs its own |
| --- | --- | --- | ---: | ---: | ---: |
| water | closed | B3LYP | 1.84e-05 | 2.09e-06 | |
| water | closed | CAM-B3LYP | 1.86e-05 | **2.35e-06** | 2.36e-06 |
| water | closed | wB97X-D4 | 1.87e-05 | 2.70e-06 | |
| hydroxyl | open | B3LYP | 1.47e-05 | 4.36e-06 | |
| hydroxyl | open | CAM-B3LYP | 1.55e-05 | **4.16e-06** | 4.11e-06 |
| hydroxyl | open | wB97X-D4 | 1.48e-05 | 4.31e-06 | |

The range separated rows land on the four-centre noise floor to two digits, closed
shell and open, and the control sits at the same level. Two things that check needs:
`conv_thresh` of 1e-9 and not tighter, because plain B3LYP water does not converge at
1e-10 in two hundred iterations; and an off-equilibrium geometry, or the whole
gradient is 1e-02 and its terms cancel.

## Water clusters, where the fitted path is the slower one

Every molecule measured above this point is compact: caffeine, tagrisso, taxol, c60.
This section is the first extended, sparse system in the file, and it reverses the
result. **On water clusters the resolution of the identity is overtaken by the
four-centre build**, at 47 waters in def2-svp and at 76 in def2-tzvp, and the margin
widens from there.

B3LYP, `def2-universal-jkfit`, convergence 1e-6, one node of the EPYC, 8 ranks of 32
threads, `in_memory` throughout. The runner is `node_b3lyp.py`, which is **not in this
repository** -- it lives beside the geometries on the node -- and the rows below are
from its stdout rather than from a records file, so unlike every other table here they
cannot be re-rendered. That is a gap, not a choice.

### def2-svp, where naux is 4.71 times nao

| waters | nao | naux | RI-JK Fock | 4c Fock | RI-JK whole | 4c whole | who wins |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 5 | 120 | 565 | 1.08 | 1.66 | 1.87 | 2.57 | RI 1.37x |
| 10 | 240 | 1130 | 1.14 | 1.45 | 2.58 | 2.80 | RI 1.09x |
| 20 | 480 | 2260 | 2.01 | 2.38 | 4.93 | 5.13 | RI 1.04x |
| 32 | 768 | 3616 | 4.35 | 5.94 | 10.35 | 11.16 | RI 1.08x |
| 47 | 1128 | 5311 | 13.99 | 13.83 | 26.65 | 23.08 | **4c 1.15x** |
| 76 | 1824 | 8588 | 60.10 | 42.74 | 104.64 | 63.62 | **4c 1.64x** |
| 100 | 2400 | 11300 | 159.43 | 83.81 | 262.15 | 120.50 | **4c 2.18x** |
| 139 | 3336 | 15707 | 521.70 | 191.41 | 832.28 | 255.35 | **4c 3.26x** |
| 190 | 4560 | 21470 | 1581.51 | -- | 2559.56 | -- | |

### def2-tzvp, where it is 2.63 times nao

| waters | nao | naux | RI-JK Fock | 4c Fock | RI-JK whole | 4c whole | who wins |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 5 | 215 | 565 | 1.27 | 1.68 | 2.44 | 2.68 | RI 1.10x |
| 10 | 430 | 1130 | 1.47 | 2.47 | 3.13 | 4.30 | RI 1.37x |
| 20 | 860 | 2260 | 3.79 | 11.99 | 9.61 | 19.21 | **RI 2.00x** |
| 32 | 1376 | 3616 | 13.98 | 33.37 | 27.69 | 45.13 | RI 1.63x |
| 47 | 2021 | 5311 | 49.30 | 84.00 | 82.70 | 106.10 | RI 1.28x |
| 76 | 3268 | 8588 | 207.72 | 260.66 | 338.90 | 310.22 | **4c 1.09x** |

### The advantage peaks and then decays, in both bases

Read the last column of either table downward. It is not a ratio which is simply small;
it is a ratio which **rises, turns and goes through one**. The larger basis moves the
turn later -- the peak from 5 waters to 20, the crossing from 47 to 76 -- and changes
nothing about its shape.

The reason is in the exponents, fitted between the last two clusters of each series:

| | four-centre | RI-JK | gap |
| --- | ---: | ---: | ---: |
| def2-svp, Fock | 2.51 | 3.60 | 1.09 |
| def2-svp, whole | 2.28 | 3.51 | 1.23 |
| def2-tzvp, Fock | 2.36 | 2.99 | 0.63 |
| def2-tzvp, whole | **2.23** | **2.93** | **0.70** |

**The four-centre exponent is 2.2 in both bases.** It is a property of the screening on
a hydrogen-bonded cluster and not of the basis, and it is nothing like the 3.97 to 4.04
the compact molecules gave: on a sparse system the dense build is most of two powers
better than its textbook form.

**The fitted exponent fell from 3.5 to 2.9 when naux/nao fell from 4.71 to 2.63.** So a
larger orbital basis improves the power and not merely the prefactor, which is the one
encouraging number here. It is still 0.70 above the dense path, and a positive gap
postpones a crossing without ever removing it.

### What is losing it is not the integrals

At 76 waters in def2-tzvp the **Fock build is still 1.25 times faster fitted** and the
whole calculation is 1.09 times slower. Two things in between:

| | four-centre | RI-JK |
| --- | ---: | ---: |
| iterations, every def2-tzvp row | 12 | **14** |
| `outside` as a share of the wall, 76 waters | 12.4% | **34.9%** |

`outside` is the whole calculation less the Fock builds, and for the fitted path it
holds the metric and the B vectors along with the ordinary one-electron work. Its
exponent at def2-svp climbs the whole way -- 2.86, 3.24, 3.47, **3.74** at the largest
pair -- and by 190 waters it is 959 seconds of a 2560 second calculation. **At the top
of the def2-svp series the setup is growing faster than the Fock build it exists to
serve.** Its internal split is taken apart in "Ninety-four per cent of the B vectors was not
the integrals either", at the end of this file: the metric is one to two per cent of
it, the three-center integrals five to eight, and everything else is one matrix
product applying the metric -- of which two thirds is multiplying structural zeros
at this size.

### What the runs cost in accuracy, and what they say about the node

The fitting error is flat per water and does not drift with cluster size:

| basis | hartree per water |
| --- | ---: |
| def2-svp, 20 waters and up | 6.5 to 6.7e-06 |
| def2-tzvp, 32 waters and up | **1.91e-05**, to three digits in three clusters |

Three times larger in the bigger basis, which is the same auxiliary set fitting a 1.79
times larger orbital basis.

One cluster was run twice by accident and is worth more than that suggests: 32 waters,
def2-tzvp, fitted path, on n260 and on n353. **Identical energy to every digit, 26.52
seconds against 27.69 -- 4.4 per cent apart.** That is the only direct measure of
node-to-node variation in this file, and it bounds the cross-node pairs below.

### Which rows are cross-node

The series was measured over four nodes as allocations came free. Pairs at 32 waters and
above in def2-tzvp are same-node; in def2-svp everything above 47 waters is not.

| rows | node |
| --- | --- |
| def2-svp: 4c 5-47 waters, RI 10-47 | n322 |
| def2-svp: RI 5 waters, 76-190 | n248 |
| def2-svp: 4c 76-139. def2-tzvp: both ways 5-20 | n260 |
| def2-tzvp: both ways 32-76 | n353 |

The curve shows no discontinuity at any join, and the 4.4 per cent above says what the
join is worth. It is still a consistency argument and not a control.

### What this does and does not say

It says that on an extended, sparse, hydrogen-bonded system the fitted path has a
bounded advantage which is spent by the time the calculation is large enough to care
about, and that the bound moves with the basis but the shape does not.

It does not say the resolution of the identity is slower in general, and the rest of
this file is the counter-example: caffeine gives 3.75 at def2-svp and 19.20 at
def2-tzvpd, gradients up to 59.6, because a compact molecule gives the dense build
nothing to screen. **Water clusters are the adversarial case for this method and were
chosen as one.** The def2-tzvpd and def2-qzvp columns, where naux/nao falls to 1.95 and
0.97, have not been measured and are where the gap should be narrowest.

### A prediction which was wrong three times, in the same direction

Every four-centre extrapolation made while this series was running came in high: about
10x for a def2-svp to def2-tzvp step which measured 5.04x, and 56 seconds for the 32
water four-centre run which measured 45.13. The error each time was reaching for the
p near 4.0 that the compact molecules gave, on a system whose screening puts it at 2.2
to 2.4. The fitted path's extrapolations were accurate to 2 per cent over the same
range. See "Benchmark runtime estimates": the failure is analogy across molecule types,
and the fix is to fit the exponent on the series being measured.

### The dense control, and the factor of eight that separates the two

c60 was run at the same sizes on the same node as the clusters above, n353, both ways
and both bases. It is the opposite system: sixty heavy atoms in a shell, nothing at a
distance from anything, the least screenable molecule in this file.

| | nao | naux | iterations | Fock | whole | `outside` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| c60 def2-svp, RI-JK | 840 | 4500 | 17 | 9.13 | 19.96 | 5.88 |
| c60 def2-svp, four-centre | 840 | -- | 15 | 35.08 | 44.03 | 4.69 |
| c60 def2-tzvp, RI-JK | 1860 | 4500 | 18 | 37.92 | 79.42 | 22.49 |
| c60 def2-tzvp, four-centre | 1860 | -- | 15 | **907.37** | **935.83** | 14.53 |

| at comparable size | nao | naux | Fock | whole |
| --- | ---: | ---: | ---: | ---: |
| **c60, def2-tzvp** | 1860 | 4500 | **RI 23.93x** | **RI 11.78x** |
| water 32-mer, def2-tzvp | 1376 | 3616 | RI 2.39x | RI 1.63x |
| water 47-mer, def2-tzvp | 2021 | 5311 | RI 1.70x | RI 1.28x |
| **c60, def2-svp** | 840 | 4500 | **RI 3.84x** | **RI 2.21x** |
| water 32-mer, def2-svp | 768 | 3616 | RI 1.37x | RI 1.08x |
| water 47-mer, def2-svp | 1128 | 5311 | 0.99x | 4c 1.15x |

**At the same nao and the same naux, the dense molecule gives the fitted path eight to
nine times more than the sparse one does.** Nothing else differs. The c60 def2-tzvp
four-centre run is the most integral bound calculation in this file: 907 seconds of
Fock build in a 936 second calculation, ninety-seven per cent.

### The mechanism, in one number

The step from def2-svp to def2-tzvp holds naux fixed at 4500 while nao goes 840 to
1860, which is the one experiment the cluster series cannot do -- there nao and naux
always move together.

| p in nao, at fixed naux | four-centre | RI-JK |
| --- | ---: | ---: |
| c60, dense | **4.09** | 1.79 |
| water clusters, sparse | 2.77 | 1.09 |

**4.09 is the textbook fourth power.** On c60 the dense build gets nothing from
screening and runs at its formal cost, which is also what every compact molecule in
this file gave: 3.97 to 4.04. On water it runs at 2.77. That difference of one and a
third powers is the screening, and it is the entire water cluster result. The fitted
path moves by 0.7 of a power between the two systems, which is small by comparison:
**this is a story about what the denominator can do, not about what the numerator
cannot.**

### The preparation, separated at last

`outside` for the fitted path less `outside` for the dense one leaves the RI
preparation, everything else in there being common work:

| case | naux | prepare | share of the fitted wall |
| --- | ---: | ---: | ---: |
| water 32-mer, def2-svp | 3616 | 0.75 s | 7.2% |
| c60, def2-svp | 4500 | 1.19 s | 6.0% |
| water 47-mer, def2-svp | 5311 | 3.29 s | 12.3% |
| c60, def2-tzvp | 4500 | 7.96 s | 10.0% |
| water 47-mer, def2-tzvp | 5311 | 10.53 s | 12.7% |
| water 76-mer, def2-tzvp | 8588 | 79.57 s | 23.5% |
| water 139-mer, def2-svp | 15707 | **246.73 s** | **29.6%** |

**It tracks naux and not the system.** Six to thirteen per cent while naux is under six
thousand, on either kind of molecule, and a third of the calculation once naux reaches
five figures -- which is what a large sparse system in a small basis produces. So the
`outside` share reported for the clusters above is not a property of water; it is a
property of having twenty thousand auxiliary functions.

Two cautions. The subtraction is an upper bound: the fitted path runs two to three more
iterations and `outside` holds per-iteration work as well, which for the small rows is
most of the difference. And at fixed naux the preparation still grows with nao -- c60
pays 6.7 times more for 2.2 times the orbital basis -- so it is not metric-bound
either, the metric depending on naux alone. What fits is the B vectors, which the profile at the end of this file confirms
directly.

### The iteration count is general, and nobody has explained it

| | four-centre | RI-JK |
| --- | ---: | ---: |
| c60, both bases | 15 | 17 to 18 |
| water clusters, every def2-tzvp row | 12 | 14 |

**Two to three more iterations everywhere**, on the densest molecule here and on the
sparsest, at both bases. Convergence is 1e-6 in all of them and the guess is the same.
That is a flat penalty of roughly fifteen per cent which is not integrals, not
screening and not scaling, and it is the one cost in this section that is paid by every
calculation in the file rather than by the adversarial ones. It has not been
investigated.

### What the two systems together say

The fitted path has two separate problems and they belong to different calculations.

**The exchange exponent is a sparse-system problem.** On c60 the build is 23.93 times
ahead and there is nothing to fix. It only bites where the dense build screens down to
2.2, which is where the ratio of exponents, not the ratio of times, decides the race.

**The preparation is a large-naux problem**, reached by big systems in small bases on
either kind of molecule.

**The iterations are everybody's problem** and the cheapest of the three to look at.

The fitting error, for the record: 7.0e-06 hartree per atom for c60 in def2-svp and
5.3e-06 in def2-tzvp.

## Ninety-four per cent of the B vectors was not the integrals either

The water clusters above put a third of the fitted calculation in `outside`, the
part which is not a Fock build, and named the setup as the likely cause. This is the
profile which took that apart. It reaches the same shape of answer as "Ninety-two
per cent of the gradient was not the integrals", one phase over, and this time the
phase it lands on cannot be fixed.

Everything here is switched on by `VLX_RIJK_PROFILE` in the environment, which was
already in the driver and had never been run on a system this size.

### The setup is the B vectors, and the metric is nothing

| | ranks | naux | metric | B vectors | setup as a share of the wall |
| --- | ---: | ---: | ---: | ---: | ---: |
| 20 waters, def2-svp | 1 | 2260 | **0.02 s** | 1.71 s | 13.0% |
| 32 waters, def2-svp | 1 | 3616 | **0.08 s** | 10.90 s | 21.9% |
| 76 waters, def2-svp | 8 | 8588 | **0.32 s** | 24.47 s | 24.3% |

The metric is formed and inverted **on the master rank alone** and broadcast
(`scfdriver.py:1945`), with every other rank idle through it, and the two-center
integrals behind it are unscreened by design -- the Coulomb operator does not fall
off, so no pair of auxiliary atoms is negligible. It was the obvious suspect and it
was named as such twice. It is **one to two per cent of the setup** at every size
measured, on one rank and on eight. Nothing about it is worth changing.

### Inside the B vectors

| | integrals | contract | rest |
| --- | ---: | ---: | ---: |
| 20 waters, 1 rank | 0.160 s, 8.8% | 1.615 s, **89.3%** | 1.9% |
| 32 waters, 1 rank | 0.554 s, 5.1% | 10.203 s, **93.9%** | 1.0% |
| 76 waters, 8 ranks | 2.02 s, 8.2% | 21.83 s, **90.4%** | 1.4% |

**The SIMD three-center integrals are five to eight per cent of forming the B
vectors**, and the whole of the rest is one matrix product applying the metric to
them. The kernels this file spends hundreds of lines on are not the cost of the
setup any more than they were the cost of the gradient.

### The contraction is not badly written

At 32 waters it moves 8.05 GB of gathered values against 3616 metric rows, which is
7.27e12 flops, in 10.20 seconds: **713 Gflop/s where Accelerate's `dgemm` peaks at
846 on those eight threads, 84 per cent.** It is threaded, it is on the right
routine, and it is close to what the machine can do.

So there is no arithmetic to win back by running it better. Only by not doing it.

### Two thirds of it is multiplying zeros

A chunk is gathered to the width of the widest auxiliary group and multiplied whole,
while each group fills only the atom pairs it keeps. The rest is zeroed and
multiplied:

| | naux | padding | of which absent groups | of which tails |
| --- | ---: | ---: | ---: | ---: |
| 10 waters | 1130 | 18.7% | 0.0% | 18.7% |
| 20 waters | 2260 | 31.8% | 0.0% | 31.8% |
| 32 waters | 3616 | 44.7% | 0.0% | 44.7% |
| 76 waters | 8588 | **66.9%** | -- | -- |
| 139 waters | 15707 | **78.7%** | -- | -- |

**All of it is tails and none of it is absent groups.** No auxiliary group is ever
missing outright; they run out at different places, so the filled region is a
staircase and the product multiplies its bounding rectangle. The share grows with
the system, which is why the setup's exponent outruns the Fock build's at the top of
the def2-svp series.

### What cutting the chunk would save, and why it cannot be had

Cutting the chunk across the pairs into bands, each multiplied against only the rows
which reach it, would recover most of that. The driver counts what each choice would
save, in the order the metric is already permuted into and in a per block order which
would need the metric gathered block by block:

| bands | 32 waters | 76 waters | 139 waters | per block, 139 waters |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 22.7% | 37.7% | 46.8% | 49.3% |
| 4 | 31.2% | 52.4% | 64.1% | 66.9% |
| 8 | 36.2% | 58.5% | 71.1% | 74.1% |
| 16 | 39.0% | 61.5% | 74.3% | 77.3% |

The global order captures 93 to 96 per cent of the per block one, so the expensive
half of the idea is unnecessary.

**It was implemented and it is slower.** At 32 waters the contraction goes 10.299 s
at one band to 10.077 at two, **13.618 at four and 25.083 at eight** -- against a
predicted 1.28x at four. The change was reverted.

**The reason is in the counters and was there all along.** The profile prints
`chunk 579`, which is `nchunk`, the largest a chunk may be. What the products
actually get is `ngathered / (8 x ncols x nproducts)`:

| | ngathered | ncols | products | average count |
| --- | ---: | ---: | ---: | ---: |
| 32 waters, laptop | 8.05 GB | 3616 | 9191 | **30.3** |
| 76 waters, node | 107 GB | 8588 | 51885 | **30.1** |

**The average product is 3616 x 30 x 3616.** Thirty columns do not divide into four
bands, let alone eight. Every BLAS measurement taken to choose a band count was made
at 579 columns, a shape which essentially never occurs in this loop, and every one of
them was therefore answering a question about a different computation.

### The BLAS measurements, which are sound and do not apply

They are kept because they say something about the library which will matter again.
`blas_trmm_check.py`, one rank, the shapes named rather than the shapes that occur:

| naux | library | `dtrmm` | 2 bands | 4 bands | 8 bands | 16 bands |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 3616 | Accelerate, 8 threads | 98% | 94% | 88% | 74% | 53% |
| 8588 | OpenBLAS, 32 threads | **12%** | 86% | 79% | 67% | 47% |
| 15707 | OpenBLAS, 32 threads | **11%** | 84% | 62% | 39% | 21% |

**`dtrmm` is barely threaded in the node's OpenBLAS** -- 338 Gflop/s against
`dgemm`'s 2776, flat across every size, which is about four cores' worth. Applying a
triangular factor with the routine written for it would be two and a half to four
times slower than the general one. On Accelerate it is 96 to 99 per cent, so the
laptop would have given a clean and entirely misleading green light.

### The exchange build has no sparsity to exploit either

The same profile answered a separate question. `CSimdRIFockDriver` carries a sparse
alternative to the dense half transformation, chosen by `_dense_threshold`, which
defaults to zero -- "expand always" -- and which nothing in `pymodule` had ever set.
It was exposed as `ScfDriver.ri_dense_threshold` and run:

| | B vector density |
| --- | ---: |
| 5 waters | 0.9142 |
| 10 waters | 0.8053 |
| 20 waters | 0.6750 |
| 32 waters | **0.5468** |
| c60 | **0.8679** |

**The B vectors are half full on the sparsest system and seven eighths full on a
compact one.** The walk saves at most one over the density in arithmetic, under
three times, and gives up the matrix unit to do it. Measured at 76 waters on the
node the two are a tie: 57.24 s against 56.35 for the Fock builds, with the energies
identical to every digit -- the first time that code path had executed at all.

So the dense default was right, and it now has a number behind it rather than an
inference from laptop sized problems.

### Three routes closed, and what is left

| route | why not |
| --- | --- |
| apply the triangular metric with `dtrmm` | 11 to 12% efficiency on the node's OpenBLAS |
| exploit the triangle by row panels of `dgemm` | needs the metric un-permuted, and only 1.18x even then |
| band the chunk across the pairs | the chunks are thirty columns wide |
| the sparse exchange half transformation | the B vectors are 55 to 87% dense |

What remains is not an optimisation. The contraction applies a **dense** metric, of a
dimension 4.7 times the orbital basis at def2-svp, to every surviving pair; that is
the naux squared per pair which sets the exponent, and no arrangement of the same
arithmetic removes it. Either the products are made wider first -- batching tasks
which share a pattern, which would improve the `n = 30` shape and create columns
worth banding -- or the metric is made sparse, which means local or robust fitting
and is a different method rather than a faster one.

### Two mistakes worth keeping

**A benchmark of the wrong shape passed every check.** `blas_trmm_check.py` verified
its own arithmetic, reported plausible speedups, and was measured on the node rather
than by analogy. It was still wrong, because the shape it was given came from a
number in the profile which does not mean what it looks like. The counter which would
have caught it -- the average `count` -- was in the same output the whole time.

**A stale build reported a plausible wrong answer.** After the revert the per block
column read equal to the global one, which the source cannot produce. `make` had not
rebuilt: the object was older than the restored source. It was caught only because
those two numbers had been measured before and a wrong one was recognisable. See
"Detecting make failures".

## One product of five hundred columns instead of seventeen of thirty

The profile above leaves the contraction which forms the B vectors as the largest
phase of a fitted calculation, 90 to 94 per cent of a setup which is a quarter of the
wall, running at a good fraction of what the library gives and doing two thirds of its
arithmetic against structural zeros. Four ways at those zeros were measured and all
four failed. **The change which did work is not about the zeros at all.**

### What was wrong with the shape

A task is one combination of basis functions of one pair block, and the pairs of a
block are swept in chunks of at most `nchunk`. The gathered buffer is allocated for
`nchunk` columns -- 579 at 32 waters, 434 on the node -- because the widest block
needs that many. **The average block has thirty atom pairs.** The buffer was
ninety-five per cent empty on every call, and the library was being handed a matrix
product of thirty columns, which is a matrix product in name only:

    ngathered / (8 x ncols x nproducts)

| | before | after |
| --- | ---: | ---: |
| 20 waters, laptop | 30.4 | **836** |
| 32 waters, laptop | 30.3 | **528** |
| 76 waters, node | 30.1 | **392** |

The chunks of several tasks are now gathered side by side into the one buffer and
handed to a single product. Every task of a batch shares the metric and the number of
rows, so their chunks are column slices of one matrix and always were; nothing else
had to change to put them there. It costs no memory which was not already reserved.

### What it buys

Both builds run back to back on one machine in one session, the unbatched one made by
capping a group at one item:

| | products | one per group | batched | |
| --- | ---: | ---: | ---: | ---: |
| 20 waters, laptop | 3551 -> 131 | 1.527 s | 1.294 s | **1.18x** |
| 32 waters, laptop | 9191 -> 527 | 10.267 s | 8.165 s | **1.26x** |
| 76 waters, node | 51885 -> 3978 | 22.03 s | 17.195 s | **1.28x** |

It grows with the system, which is what it should do: a larger molecule has more tasks
to gather together and more to gain from the wider shape.

**The arithmetic is identical and the counters prove it.** The gathered volume is
unchanged to the hundredth of a gigabyte -- 1.98, 8.05 and 107.04 GB before and after
-- the padding is unchanged at 31.8, 44.7 and 66.9 per cent, the energies of four
clusters agree to 1.3e-11 or better, and every iteration count is the same. The same
work, reshaped.

On the node the setup fell from 24.47 s to 20.37. The wall clock did not move, because
the Fock build came out ten per cent slower on a different node in the same run, which
batching does not touch; that is node to node variation larger than the 4.4 per cent
measured elsewhere in this file, and it is why the table above is the laptop's.

### Batching and banding are alternatives, not a sequence

The counter which says what cutting the chunk into bands would save was rewritten to
ask the question of a group rather than of a chunk. Against a batched group it answers
**zero, at every band count, in both orderings**.

The reason is the pooling itself. A band of a batched group spans several tasks, and
their staircases do not line up: for any band, some task in it still has an entry
reaching, so the band needs every row. What made banding look attractive was a
structure which batching destroys.

That also settles which of the two to have. The 1.66x once quoted for banding
multiplied savings measured at thirty columns by an efficiency measured at a hundred
and eight -- two numbers from different computations. Banding at the shape which
actually occurs was implemented and measured at **0.76x**, and batching at **1.28x**.

### Five routes to the padding, and where it stands

| route | result |
| --- | --- |
| `dtrmm` on the triangular metric | 11 to 12% efficiency on the node's OpenBLAS |
| row panels of `dgemm` | needs the metric un-permuted, 1.18x at best |
| bands across the pairs, unbatched | 0.76x: thirty columns do not divide |
| bands across the pairs, batched | no saving available at all |
| the sparse exchange half transformation | the B vectors are 55 to 87% dense |

**The padding is untouched and remains 45 to 79 per cent of the contraction.** What
was won is the shape, not the sparsity. Reaching the sparsity needs either the pooled
columns sorted by how deep they reach, which would restore a staircase across the
group and has not been measured, or a kernel of our own, which must beat 553 Gflop/s
-- a fifth of the node's `dgemm` peak -- now that the baseline is 1.28 times faster
than it was.

Neither is likely to matter as much as the exponent. `prepare` is a fifth of the wall;
the exchange scaling is what decides whether the fitted path is usable on an extended
system at all, and none of this touches it.
