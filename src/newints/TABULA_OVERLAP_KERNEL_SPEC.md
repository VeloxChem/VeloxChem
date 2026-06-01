# Tabula spec — newints spherical overlap kernels (current, no screening)

Authoritative description of the 49 auto-generated two-center overlap kernels in
`src/newints/` (namespace `ovlab`, files
`OverlapABRec{S,P,D,F,G,H,I}{S,P,D,F,G,H,I}.{hpp,cpp}`), as built on branch
`new-os-ri` (HEAD 80b6cd67d). Regenerate from this; do not hand-edit the emitted
files. Supersedes the earlier per-step specs (codegen-bug fix, buffer interface,
primitive-screening) — in particular **primitive-pair screening has been REMOVED**
(see the note at the end).

## Signature (all 49)

```cpp
auto overlap_<a>_<b>(
    const CBasisFunction &bra,
    const CBasisFunction &ket,
    const TPoint<double> &bra_center,
    const TPoint<double> &ket_center,
    double                *buffer) -> void;
```

- `<a>`,`<b>` in {s,p,d,f,g,h,i}; the kernel writes the (2*l_a+1) x (2*l_b+1)
  block row-major into `buffer`. Returns void.
- No `cutoff` / `threshold` parameter (no screening — see end).
- Headers include `BasisFunction.hpp`, `Point.hpp`, and (for l_a+l_b > 0)
  `RealSolidHarmonicAB.hpp`; not `SparseMatrix.hpp` (no `newints::Block` use).

## Buffer contract

- `buffer` points to at least `(2*l_a+1)*(2*l_b+1)` doubles, row-major.
- The kernel writes **every** entry (verified: each kernel assigns exactly
  nrows*ncols positions). The caller need not pre-zero — and indeed the driver's
  per-thread arena is non-zeroing, so a kernel that leaves an entry unwritten
  would emit garbage. Do not rely on zero-initialization.

## Phase 1 — geometry

From `bra_center`/`ket_center`: AB_x/y/z, the needed AB_* products, and
`R2 = AB_x^2 + AB_y^2 + AB_z^2`.

## Phase 2 — V evaluation + primitive contraction

Plain double loop over primitive pairs (no screening, no early exit). Use the
non-copying const-ref accessors (NOT the by-value `get_*`):

```cpp
const auto &exps_a  = bra.exponents();
const auto &coefs_a = bra.normalization_factors();
const auto &exps_b  = ket.exponents();
const auto &coefs_b = ket.normalization_factors();
const auto pi = mathconst::pi_value();
std::array<double, K> V = { ... };          // K = min(l_a,l_b)+1 reduction terms
for i in bra primitives:
    alpha = exps_a[i]; ca = coefs_a[i];      // + alpha powers as needed
    for j in ket primitives:
        beta = exps_b[j]; cb = coefs_b[j];   // + beta powers
        p = alpha+beta; pinv = 1/p;          // + pinv powers
        ss     = (pi*pinv) * sqrt(pi*pinv) * exp(-alpha*beta*pinv*R2);  // (s|s)
        cab_ss = ca*cb*ss;
        V[k] += cab_ss * (alpha,beta,pinv powers per the V[k] definition)
```

SS base case has no V array: accumulate `sab += ca*cb*ss` and write
`buffer[0] = sab;`.

## Phase 3 — fused M.V -> spherical block

Evaluate the needed `harm::Y_ll_<l>_m_<m>(AB_x,AB_y,AB_z)` solid harmonics once,
then write every block entry as the fused contraction of the harmonics with `V`:

```cpp
auto *d = buffer;
d[..] = <M.V expression in the Y_ll_* and V[k]>;   // all nrows*ncols entries
```

`RealSolidHarmonicAB.hpp` (namespace `harm`, inline) must define all real solid
harmonics up to the maximum order any kernel references (= up to l_a+l_b; for the
S..I grid that is l=12). Every `sqrt(N)` constant used in Phase 3 must be declared
in the kernel preamble.

## REMOVED: primitive-pair screening

Earlier revisions added a per-primitive-pair screen inside Phase 2 — first a
fixed `cutoff` on the exponent argument (`if (arg > cutoff) continue;`), then a
reverse-order loop with an early `break` on the full `c_i c_j (pi/p)^{3/2} R^L
exp(-mu R^2)` magnitude, with a `cutoff`/`threshold` kernel parameter. Both were
benchmarked and **reverted** because they gave no net speedup (geomean ~1.00x)
and the reverse-break also loosened accuracy to ~1e-13. Reason: at production
thread counts the overlap build is throughput / memory-bandwidth / merge bound,
not exp-bound, so cutting per-primitive-pair work does not move wall-clock; the
shell-pair coarse screener (in the driver, not the kernels) already removes whole
distant blocks. Do NOT reintroduce a kernel-level screen or `cutoff`/`threshold`
parameter without a profile showing a compute-bound regime.

## Validation after regeneration

`make -C src release`; the newints tests + `tests/test_basis_function.py` must
pass; `newints.OverlapDriver(...).to_dense()` must match the legacy
`veloxchemlib.OverlapDriver` to machine precision (~1e-15) across STO-3G and
cc-pVDZ..cc-pV6Z. See [[newints-integrals-subsystem]] and
[[newints-overlap-benchmark]].
