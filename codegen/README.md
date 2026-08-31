# Generators of the SIMD integral kernels

The kernels under `src/simd_func` and `src/simd_t2c_kinetic_energy` are generated.
These are the generators which produced the ones added here; the generator of the
overlap kernels under `src/simd_t2c_overlap` is not in this repository.

Run them from anywhere: each finds the checkout from its own location and writes
into `src/`.

## make_solid_harmonics.py

Emits one order of the solid harmonics recursion, in the shape of the orders
which already existed. The coefficients follow

    fact     = sqrt((2l - 1) / (2l))                        for m = +-l
    fz_(l-1) = sqrt(2l - 1)                                 for m = +-(l - 1)
    fz_m     = sqrt((2l - 1)^2 / ((l + m)(l - m)))          for |m| <= l - 2
    fr_m     = sqrt((l - 1 + m)(l - 1 - m) / ((l + m)(l - m)))

**It reproduces the existing order twelve byte for byte**, source and header,
which is how it was validated before orders thirteen and fourteen were emitted
with it. Regenerate order twelve and diff against `SimdHarmonicsRecQ.cpp` after
any change to it.

## make_kinetic_formula_kernels.py

Emits the thirteen kinetic energy kernels which follow a closed formula:

    (s|T|s)     = [ 3 mu - 2 mu^2 R^2 ] (s|S|s)
    (s|T|l)_m   = [ (2l + 3) mu - 2 mu^2 R^2 ] ( alpha / p)^l S_{l,m}(AB) (s|S|s)
    (l|T|s)_m   = [ (2l + 3) mu - 2 mu^2 R^2 ] (-beta  / p)^l S_{l,m}(AB) (s|S|s)

The ratio belongs to the side carrying no angular momentum, which is easy to get
backwards: `(s|.|l)` takes `alpha / p` and `(l|.|s)` takes `-beta / p`, as the
overlap kernels of the same combinations do.

## make_kinetic_recursion_kernels.py

Emits the thirty six kernels of two non-zero angular momenta from the recursion
data. Pass the combinations as arguments and point `KINETIC_DATA` at the
directory holding `kinetic_harmonics_*.txt`:

    KINETIC_DATA=/path/to/data python3 make_kinetic_recursion_kernels.py pp pd dd ...

Each term of the data is `(num/den) sqrt(radicand) E G (s|S|s)`, where `E` is
`alpha^a beta^b mu^m p^q` and `G` is `R^{2k} S_{L,M}(AB)`. Neither `R^{2k}` nor
`S_{L,M}` depends on the pair of primitives, so one accumulator per `E` suffices
and the terms are combined against the harmonics afterwards, which is the shape
of the overlap `(p|p)` kernel.

`LOAD_BUDGET` is the one number worth understanding. The angular components are
formed in several loops, and a loop is closed once the values it loads reach the
budget. Grouping by a fixed count of components instead let the large kernels
load sixty values against the thirty two vector registers of the machine, and the
vectorizer then emitted **scalar code for the whole combine**: `(h|T|h)` and
`(i|T|h)` fell to eight vector operations in three thousand lines. Check the
emitted code after changing it, and see the note on how to check vectorization,
because the obvious way reports zero for code which is fully vectorized.

## What every generator must keep

The buffer of a kernel is sized from `*std::ranges::max_element(dimensions)`, not
from `dimensions.back()`. The primitives are sorted by descending exponent, but a
pair's reach carries its prefactor as well as its decay, so the last pair is not
always the furthest reaching, and taking the last one writes past the buffer. That
bug reached the repository once through the overlap kernels and was found only by
AddressSanitizer on a diffuse basis.
