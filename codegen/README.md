# Generators of the SIMD integral kernels

The kernels under `src/simd_func`, `src/simd_t2c_kinetic_energy` and
`src/simd_t2c_electron_repulsion` are generated,
and these are the generators which produced them. The generator which produced the
overlap kernels under `src/simd_t2c_overlap` is not in this repository, and the
recursion generator here reads the overlap data but does not reproduce those files
byte for byte; see below.

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

## make_coulomb_formula_kernels.py

Emits the thirteen two-center electron repulsion kernels which follow a closed
formula:

    (s|J|s)     = (s|J|s)^(0)
    (s|J|l)_m   = ( alpha / p)^l S_{l,m}(AB) (s|J|s)^(l)
    (l|J|s)_m   = (-beta  / p)^l S_{l,m}(AB) (s|J|s)^(l)

where the auxiliary integral of order n is

    (s|J|s)^(n) = 2 pi^(5/2) N_a N_b F_n(mu R^2) / (alpha beta sqrt(p))

with `p = alpha + beta` and `mu = alpha beta / p`. The same rule for the ratio
applies as for the kinetic energy: it belongs to the side carrying no angular
momentum.

These differ from the overlap and the kinetic energy in two ways. Nothing is
screened, so there are no column dimensions, no `nmax` and no threshold, and
every row spans the atom pairs of the block. And the exponential of a pair of
primitives is replaced by the Boys function, which `simdfunc::compute_boys_function`
evaluates for all pairs of primitives of a kernel in one call, on a
`CSimdVariableMatrix` with the arguments in block zero and the orders zero to `l`
in the blocks above it. Only order `l` is read; the lower orders are formed on
the way to it.

Verified against `CTwoCenterElectronRepulsionDriver` on water in cc-pVQZ, cc-pV5Z
and cc-pV6Z, which reaches every angular momentum up to six: the largest
deviation is 8.5e-13 on values up to 1.9e2, i.e. 5e-15 relative.

## make_recursion_kernels.py

Emits the kernels of two non-zero angular momenta from the recursion data in
`data/`, for the overlap, the kinetic energy and the two-center Coulomb
integrals: the three share the format of the data and the shape of the generated
code, and the file reports which it is on its `KERNEL` line.

    RECURSION_KIND=kinetic python3 make_recursion_kernels.py pp pd dd ...
    RECURSION_KIND=overlap python3 make_recursion_kernels.py pp pd dd ...
    RECURSION_KIND=coulomb python3 make_recursion_kernels.py pp pd dd ...

Each term of the data is `(num/den) sqrt(radicand) E G (s|.|s)`, where `E` is
`alpha^a beta^b mu^m p^q` and `G` is `R^{2k} S_{L,M}(AB)`. Neither `R^{2k}` nor
`S_{L,M}` depends on the pair of primitives, so one accumulator per `E` suffices
and the terms are combined against the harmonics afterwards.

The overlap and the kinetic energy differ only in the name, the namespace and the
bound they screen the pairs of primitives with. The Coulomb data is `FORMAT 2`
and differs in more: its source integral `(s|J|s)^(m)` carries an order, so every
term has the order of its Boys function as a sixth field and the file declares
the range of orders it uses on its `ORDERS` line. Every exponent factor of the
data is used with one order alone, so an accumulator per exponent factor still
suffices and the order follows from the factor. Nothing is screened, so the
generated kernel has no column dimensions, no `nmax` and no threshold, and the
exponential of a pair of primitives is replaced by the Boys function, which
`simdfunc::compute_boys_function` evaluates for all pairs of primitives in one
call.

Verified against `CTwoCenterElectronRepulsionDriver` on ethanol in cc-pV6Z, where
three heavy atoms carry every angular momentum up to six and all thirty six
kernels are exercised: no combination deviates by more than 1.4e-13 relative.

**It reproduces the thirty six kinetic energy kernels in the repository byte for
byte**, which is how it is checked: regenerate them and diff.

**It reproduces the thirty six kinetic energy kernels after the Coulomb kernels
were added to it**, which is the check to run after any change: regenerate them
and diff.

**It does not reproduce the overlap kernels in the repository**, which came from
another generator which is not here. What it emits for the overlap is equivalent
term for term, and differs in the wording of the comments, in naming the
accumulated factors `ff_0` rather than `f_0`, in writing `a / (p * p)` rather
than `a / p / p`, which can differ in the last place, and in how the angular
components are grouped into loops. The overlap kernels were left as they are
rather than regenerated, so for that kernel this generator is a reference and not
the source of truth.

`LOAD_BUDGET` is the one number worth understanding. The angular components are
formed in several loops, and a loop is closed once the values it loads reach the
budget. Grouping by a fixed count of components instead let the large kernels
load sixty values against the thirty two vector registers of the machine, and the
vectorizer then emitted **scalar code for the whole combine**: `(h|T|h)` and
`(i|T|h)` fell to eight vector operations in three thousand lines. Check the
emitted code after changing it, and note that the obvious way of checking
vectorization reports zero for code which is fully vectorized.

## data/

The recursion data of both kernels, thirty six files each, as
`<kernel>_harmonics_<la><lb>.txt`. The generator reads them from here unless
`RECURSION_DATA` says otherwise.

## What every generator must keep

The buffer of a kernel is sized from `*std::ranges::max_element(dimensions)`, not
from `dimensions.back()`. The primitives are sorted by descending exponent, but a
pair's reach carries its prefactor as well as its decay, so the last pair is not
always the furthest reaching, and taking the last one writes past the buffer. That
bug reached the repository once through the overlap kernels and was found only by
AddressSanitizer on a diffuse basis.
