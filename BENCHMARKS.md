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

