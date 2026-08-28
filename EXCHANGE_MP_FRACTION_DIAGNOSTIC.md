# Exchange Mixed-Precision Work-Fraction Diagnostic

## Purpose

The exchange work-fraction diagnostic reports how the mixed-precision cut
builder partitions exchange tile pairs into FP64, FP32, and screened work. It
is disabled by default and does not alter the production workflow.

Enable it for a diagnostic run with:

```bash
VLX_EXCHANGE_FRACTION_STATS=1 vlx input.inp
```

The default run, with the variable unset or set to any value other than `1`,
does not allocate counters, launch the reduction kernel, copy counters,
synchronize for the diagnostic, or print diagnostic output.

## Definitions

For every `m` tile in every exchange `(i,k)` block, the cut builder produces:

- `prec`: number of `n` tiles assigned to FP64;
- `screen`: number of `n` tiles not removed by ERI screening;
- `n_n`: total number of `n` tiles before screening.

The diagnostic accumulates the following rank-local tile-pair counts across
all exchange combinations, GPUs, and exact-exchange interactions in one Fock
build:

```text
FP64     = sum(prec)
FP32     = sum(screen - prec)
screened = sum(n_n - screen)
computed = FP64 + FP32
total    = computed + screened
```

The reported mixed-precision work fraction is:

```text
FP32 fraction of computed work = FP32 / (FP64 + FP32)
```

This is an unweighted fraction of exchange tile pairs that survive screening.
It is not a runtime-weighted fraction, an instruction fraction, or a count of
primitive arithmetic operations. It should be used to describe the precision
partition selected by a threshold; measured timings are still required to
calculate achieved speedup.

## Implementation

`build_exchange_cuts_kernel` is unchanged. When the diagnostic is enabled,
`build_exchange_cuts_device` launches a separate shared-memory reduction after
the cut builder. Each reduction block performs only three global atomic adds,
one for each counter. The host copies and aggregates one three-counter buffer
per GPU after exchange computation has completed.

The output is labelled `rank-local`. A run using one MPI rank and multiple GPUs
already covers all GPUs in that rank. A multi-rank global result requires an
additional MPI reduction or post-processing of the output from all ranks.

## Validation

The diagnostic was validated on `guanine-8` using one GH200 node with four
GPUs and an exchange threshold of `1e-6`:

- the default run printed no fraction diagnostics;
- the default and flagged runs both converged in 10 SCF iterations to
  `-4624.0280060107 a.u.`;
- all 17 reported Fock builds satisfied `computed = FP64 + FP32` and
  `total = computed + screened`;
- the final Fock build reported `630649942` FP64 tile pairs, `3063439940`
  FP32 tile pairs, and an FP32 fraction of computed work of `0.829281`.
