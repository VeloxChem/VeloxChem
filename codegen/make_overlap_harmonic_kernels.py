"""Emits the two-center overlap kernels of one S type function and one function of
higher angular momentum, in both orders.

The Gaussian product center is displaced from the atom carrying the harmonic by
(a / p) times the vector between the atoms on the ket side and by -(b / p) times it
on the bra side, and the harmonic is homogeneous, so the integrals are the overlap
of the S type functions with that ratio raised to the power l folded into the
prefactor, times the solid harmonic of the vector between the atoms. The harmonic
does not depend on the pair of primitives, so it multiplies the accumulated
prefactor once, after the sum over the pairs of primitives is complete.

The harmonic is written out per component as a polynomial in the components of the
vector between the atoms and in its squared length, all four of which the
coordinates carry, so no recursion and no harmonics matrices are needed.
"""
import io
import os
from fractions import Fraction as F

import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src", "simd_t2c_overlap") + "/"
LICENSE = io.open(SRC + 'SimdOverlapRecSS.cpp').read().split('#include "SimdOverlapRecSS.hpp"')[0]
LIC_HPP = io.open(SRC + 'SimdOverlapRecSS.hpp').read().split('#ifndef')[0]

LETTER = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G', 5: 'H', 6: 'I'}
MOMENT = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six'}

X, Y, Z, R2 = sp.symbols("x y z r_2")


def isqrt_exact(n):
    r = int(n ** 0.5)
    for c in (r - 1, r, r + 1):
        if c >= 0 and c * c == n:
            return c
    return None


def terminating(fr):
    d = fr.denominator
    while d % 2 == 0:
        d //= 2
    return d == 1


def literal(fr):
    if terminating(fr):
        s = repr(float(fr))
        return s if '.' in s or 'e' in s else s + '.0'
    return f"{fr.numerator}.0 / {fr.denominator}.0"


def factor(fr):
    """sqrt(fr) as the files write it: exactly when it is rational."""
    p, q = isqrt_exact(fr.numerator), isqrt_exact(fr.denominator)
    if p is not None and q is not None:
        return literal(F(p, q))
    return f"std::sqrt({literal(fr)})"


def name(m):
    return "0" if m == 0 else (f"p{m}" if m > 0 else f"m{-m}")


def harmonics(lmax):
    """The real regular solid harmonics as polynomials in x, y, z and r_2, by the
    recursion the solid harmonics factory of simd_func runs."""
    s = {(0, 0): sp.Integer(1), (1, -1): Y, (1, 0): Z, (1, 1): X}
    for l in range(2, lmax + 1):
        for m in range(-l, l + 1):
            if m == l:
                c = sp.sqrt(sp.Rational(2 * l - 1, 2 * l))
                e = c * (X * s[(l - 1, l - 1)] - Y * s[(l - 1, -(l - 1))])
            elif m == -l:
                c = sp.sqrt(sp.Rational(2 * l - 1, 2 * l))
                e = c * (Y * s[(l - 1, l - 1)] + X * s[(l - 1, -(l - 1))])
            elif abs(m) == l - 1:
                e = sp.sqrt(sp.Integer(2 * l - 1)) * Z * s[(l - 1, m)]
            else:
                fz = sp.sqrt(sp.Rational((2 * l - 1) ** 2, (l + m) * (l - m)))
                fr = sp.sqrt(sp.Rational((l - 1 + m) * (l - 1 - m), (l + m) * (l - m)))
                e = fz * Z * s[(l - 1, m)] - fr * R2 * s[(l - 2, m)]
            s[(l, m)] = sp.expand(e)
    return s


def coefficient(c):
    """A coefficient of an expanded harmonic, which is a signed square root of a
    rational, as the kernels write it."""
    sign = "-" if c < 0 else ""
    sq = sp.Rational(sp.nsimplify(sp.simplify(c ** 2)))
    return sign, factor(F(sq.p, sq.q))


def monomial(powers):
    """The product of the variables of a monomial, by repeated multiplication."""
    out = []
    for sym, p in zip(("x", "y", "z", "r_2"), powers):
        out += [sym] * p
    return " * ".join(out)


def expression(e):
    """One component of the harmonic as a C++ expression."""
    poly = sp.Poly(e, X, Y, Z, R2)
    terms = []
    for powers, c in sorted(poly.terms(), key=lambda t: tuple(-p for p in t[0])):
        sign, mag = coefficient(c)
        mono = monomial(powers)
        if mag == "1.0":
            terms.append(f"{sign}{mono}" if mono else f"{sign}1.0")
        elif not mono:
            terms.append(f"{sign}{mag}")
        else:
            terms.append(f"{sign}{mag} * {mono}")
    out = " + ".join(terms).replace("+ -", "- ")
    return out


# NOTE: the harmonic of the higher orders does not fit one vectorized loop, as the
# 2l + 1 components together with the components of the vector between the atoms,
# its squared length and the accumulated prefactor exceed the vector registers and
# the loop starts spilling. The components are independent expressions, so they are
# split into passes of at most this many, each of which reloads only the variables
# its own components use.
CHUNK = 4


def chunks(l):
    """The angular components grouped into the passes of the harmonic loop."""
    ms = list(range(-l, l + 1))
    if len(ms) <= 9:
        return [ms]
    n = -(-len(ms) // CHUNK)
    size, rest = divmod(len(ms), n)
    out, i = [], 0
    for g in range(n):
        take = size + (1 if g < rest else 0)
        out.append(ms[i:i + take])
        i += take
    return out


def variables(exprs):
    """The variables the expressions of one pass read, in the order the loop
    declares them."""
    free = set()
    for e in exprs:
        free |= e.free_symbols
    return [(sym, nm) for sym, nm in ((X, "x"), (Y, "y"), (Z, "z"), (R2, "r_2")) if sym in free]


def emit(l, side):
    """side is 'ket' for (s|l) and 'bra' for (l|s)."""
    tag = f"S{LETTER[l]}" if side == "ket" else f"{LETTER[l]}S"
    fn = tag.lower()
    ncomps = 2 * l + 1
    lbra, lket = (0, l) if side == "ket" else (l, 0)
    momenta = f"{MOMENT[lbra]} and {MOMENT[lket]}"

    s = harmonics(l)

    # the prefactor of a pair of primitives, with the ratio raised to the power l
    if side == "ket":
        ratio = "const auto fr = pair.aexp / fexp;"
        why = ("the harmonic sits on the ket side, so the Gaussian product center is\n"
               "        // displaced from it by (a / p) times the vector between the atoms and the\n"
               f"        // prefactor carries that ratio raised to the power {MOMENT[l]}.")
    else:
        ratio = "const auto fr = -pair.bexp / fexp;"
        why = ("the harmonic sits on the bra side, so the displacement is -(b / p) times\n"
               "        // the vector between the atoms and the prefactor carries that ratio, sign\n"
               f"        // included, raised to the power {MOMENT[l]}.")

    hpp = LIC_HPP + f"#ifndef SimdOverlapRec{tag}_hpp\n#define SimdOverlapRec{tag}_hpp\n\n"
    hpp += "#include <cstddef>\n\n"
    hpp += '#include "BasisFunction.hpp"\n#include "SimdMatrix.hpp"\n\n'
    hpp += "namespace simdovl {  // simdovl namespace\n\n"
    hpp += ("/// @brief Computes the overlap integrals of a combination of basis functions of\n"
            f"/// angular momenta {momenta} on bra and ket sides.\n"
            "/// @param values The values of the combination of basis functions in the values\n"
            "/// block of the sparsity pattern.\n"
            "/// @param nvalues The number of values to compute, i.e. the number of atom pairs\n"
            "/// surviving the screening of the combination of basis functions.\n"
            "/// @param bra The basis function on bra side.\n"
            "/// @param ket The basis function on ket side.\n"
            "/// @param coordinates The coordinates of the atom pairs, as ten rows ordered by\n"
            "/// ascending interatomic distance, holding the vector between the atoms in rows\n"
            "/// six to eight and its squared length in row nine.\n"
            "/// @param threshold The screening threshold of the integrals.\n"
            "/// @note The harmonic does not depend on the pair of primitives, so it multiplies\n"
            "/// the accumulated prefactor once rather than every contribution.\n"
            f"/// @note The {ncomps} spherical components of the harmonic are written out as\n"
            "/// polynomials in the components of the vector between the atoms and in its\n"
            "/// squared length, which the coordinates carry, so no recursion is run here.\n")
    hpp += (f"auto compute_{fn}_overlap(double               *values,\n"
            f"                        const size_t          nvalues,\n"
            f"                        const CBasisFunction &bra,\n"
            f"                        const CBasisFunction &ket,\n"
            f"                        const CSimdMatrix    &coordinates,\n"
            f"                        const double          threshold) -> void;\n\n")
    hpp += "}  // namespace simdovl\n\n"
    hpp += f"#endif /* SimdOverlapRec{tag}_hpp */\n"

    rows = "\n".join(f"    auto *out_{name(m)} = buffer.data({i + 1});" for i, m in enumerate(range(-l, l + 1)))

    groups = chunks(l)

    loops = ""
    for ms in groups:
        used = variables([s[(l, m)] for m in ms])
        aligned = ", ".join([f"out_{name(m)}" for m in ms] + ["prim"] +
                            [("ab_2" if nm == "r_2" else "ab_" + nm) for _, nm in used])
        loops += f"#pragma omp simd aligned({aligned} : simd::cache_line_size())\n"
        loops += "    for (size_t k = 0; k < nmax; k++)\n    {\n"
        for _, nm in used:
            src = "ab_2" if nm == "r_2" else "ab_" + nm
            loops += f"        const auto {nm} = {src}[k];\n"
        loops += "\n"
        for m in ms:
            loops += f"        out_{name(m)}[k] = prim[k] * ({expression(s[(l, m)])});\n\n"
        loops = loops.rstrip("\n") + "\n    }\n\n"

    cpp = LICENSE + f'#include "SimdOverlapRec{tag}.hpp"\n\n'
    cpp += "#include <algorithm>\n#include <cmath>\n#include <ranges>\n#include <string>\n\n"
    cpp += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "ScreeningFunc.hpp"\n'
            '#include "SimdAlign.hpp"\n#include "SimdDimensions.hpp"\n#include "SimdPrimitives.hpp"\n'
            '#include "SimdStorage.hpp"\n\n')
    cpp += "namespace simdovl {  // simdovl namespace\n\n"
    cpp += (f"auto\ncompute_{fn}_overlap(double               *values,\n"
            f"                   const size_t          nvalues,\n"
            f"                   const CBasisFunction &bra,\n"
            f"                   const CBasisFunction &ket,\n"
            f"                   const CSimdMatrix    &coordinates,\n"
            f"                   const double          threshold) -> void\n{{\n")
    cpp += (f"    if ((bra.get_angular_momentum() != {lbra}) || (ket.get_angular_momentum() != {lket}))\n"
            "    {\n        errors::assertMsgCritical(\n"
            f'            false, std::string("SimdOverlapRec{tag}.compute_{fn}_overlap: Basis functions must be of angular momenta {momenta}"));\n'
            "    }\n\n")
    cpp += ("    if (nvalues > coordinates.number_of_columns())\n    {\n        errors::assertMsgCritical(\n"
            f'            false, std::string("SimdOverlapRec{tag}.compute_{fn}_overlap: Number of values exceeds number of atom pairs"));\n'
            "    }\n\n    if (nvalues == 0) return;\n\n")
    cpp += "    const auto nprims = bra.exponents().size() * ket.exponents().size();\n\n"
    cpp += ("    // NOTE: the pairs of primitives are screened with the threshold of the\n"
            "    // integrals divided by their number, as their contributions accumulate into\n"
            "    // a single value and the error of the sum is bounded by the number of terms.\n\n")
    cpp += ("    const auto dimensions = simdfunc::make_column_dimensions(\n"
            "        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, "
            "threshold / static_cast<double>(nprims));\n\n")
    cpp += ("    // NOTE: the buffer holds the prefactor shared by the angular components in\n"
            "    // its first row and the integrals of the components in the rows which follow,\n"
            "    // as the harmonic factors out of the sum over the pairs of primitives and\n"
            "    // multiplies the accumulated prefactor once.\n\n")
    cpp += f"    auto buffer = simdfunc::make_primitive_buffer(dimensions, {ncomps + 1});\n\n"
    cpp += ("    if (buffer.number_of_columns() == 0)\n    {\n"
            f"        simdfunc::store_components(values, nvalues, buffer, 1, {ncomps});\n\n        return;\n    }}\n\n")
    cpp += "    const auto nmax = buffer.number_of_columns();\n\n    auto *prim = buffer.data(0);\n\n"
    cpp += rows + "\n\n"
    cpp += ("    // NOTE: the components of the vector between the atoms and its squared length\n"
            "    // are carried by the coordinates, so the harmonic below is formed from rows\n"
            "    // which are already in place.\n\n")
    cpp += ("    const auto *ab_x = coordinates.data(6);\n    const auto *ab_y = coordinates.data(7);\n"
            "    const auto *ab_z = coordinates.data(8);\n\n    const auto *ab_2 = coordinates.data(9);\n\n")
    cpp += "    constexpr auto fpi = mathconst::pi_value();\n\n"
    cpp += "    // accumulate the prefactor of each pair of primitives\n\n"
    cpp += ("    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {\n"
            "        const auto ncols = pair.ncols;\n\n"
            "        const auto fexp = pair.aexp + pair.bexp;\n\n"
            "        const auto fmu = pair.aexp * pair.bexp / fexp;\n\n"
            "        const auto fovl = fpi / fexp;\n\n"
            f"        // NOTE: {why}\n\n"
            f"        {ratio}\n\n"
            "        const auto ffact = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl)"
            + " * fr" * l + ";\n\n"
            "#pragma omp simd aligned(prim, ab_2 : simd::cache_line_size())\n"
            "        for (size_t k = 0; k < ncols; k++)\n        {\n"
            "            prim[k] += ffact * std::exp(-fmu * ab_2[k]);\n        }\n    });\n\n")
    cpp += ("    // NOTE: the integrals of the angular components are the accumulated prefactor\n"
            "    // times the components of the harmonic, formed in one pass over the rows of\n"
            "    // the buffer and of the coordinates, all of which start at a cache line\n"
            "    // boundary.\n\n")
    if len(groups) > 1:
        cpp += (f"    // NOTE: the components are formed in {len(groups)} loops, as the vectorizer runs out\n"
                f"    // of registers with all {ncomps} of them in one. Only the accumulated prefactor and\n"
                "    // the vector between the atoms are loaded by more than one loop.\n\n")
    cpp += loops
    cpp += f"    simdfunc::store_components(values, nvalues, buffer, 1, {ncomps});\n}}\n\n"
    cpp += "}  // namespace simdovl\n"

    io.open(SRC + f"SimdOverlapRec{tag}.hpp", "w").write(hpp)
    io.open(SRC + f"SimdOverlapRec{tag}.cpp", "w").write(cpp)
    return tag


if __name__ == "__main__":
    import sys
    ls = [int(a) for a in sys.argv[1:]] or [2, 3]
    for l in ls:
        for side in ("ket", "bra"):
            print("wrote", emit(l, side))
