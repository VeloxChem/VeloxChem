"""Emits the two-center overlap kernels of two functions of non-zero angular
momentum.

With u = r - P, PA = alpha AB and PB = beta AB where alpha = -b / p and
beta = a / p, the integral is

    (la ma|lb mb) = (pi/p)^(3/2) exp(-mu AB^2) * I(alpha, beta, h; AB)

where I is the polynomial obtained by expanding the two harmonics about u and
replacing each monomial of u by its Gaussian moment, which introduces h = 1/(2p).
Collected by powers of h, the result is

    I = sum_k  alpha^(la-k) beta^(lb-k) h^k  P_k(AB)

with P_k a polynomial in the components of the vector between the atoms. Only the
factor in front of P_k depends on the pair of primitives, so the kernel carries
min(la, lb) + 1 accumulators, one per k, and the angular half combines them with
the P_k. That the collection has exactly this shape is asserted while generating.

When la equals lb the component block is symmetric, since the two harmonics enter
symmetrically, so only its triangle is computed and the rest is copied.
"""
import io
import os
import sys
from fractions import Fraction as F

import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src", "simd_t2c_overlap") + "/"
LICENSE = io.open(SRC + 'SimdOverlapRecSS.cpp').read().split('#include "SimdOverlapRecSS.hpp"')[0]
LIC_HPP = io.open(SRC + 'SimdOverlapRecSS.hpp').read().split('#ifndef')[0]

LETTER = {1: 'P', 2: 'D', 3: 'F', 4: 'G', 5: 'H', 6: 'I'}
MOMENT = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six'}

UX, UY, UZ = sp.symbols("u_x u_y u_z")
RX, RY, RZ = sp.symbols("x y z")
AL, BE, HH = sp.symbols("alpha beta h")

_solid = {}


def solid(l, m, X=RX, Y=RY, Z=RZ):
    """Real regular solid harmonic, the recursion of make_solid_harmonics.py."""
    if (l, m) not in _solid:
        if l == 0:
            e = sp.Integer(1)
        elif l == 1:
            e = {-1: RY, 0: RZ, 1: RX}[m]
        else:
            r2 = RX * RX + RY * RY + RZ * RZ
            if m == l:
                f = sp.sqrt(sp.Rational(2 * l - 1, 2 * l))
                e = f * (RX * solid(l - 1, l - 1) - RY * solid(l - 1, -(l - 1)))
            elif m == -l:
                f = sp.sqrt(sp.Rational(2 * l - 1, 2 * l))
                e = f * (RY * solid(l - 1, l - 1) + RX * solid(l - 1, -(l - 1)))
            elif abs(m) == l - 1:
                e = sp.sqrt(sp.Integer(2 * l - 1)) * RZ * solid(l - 1, m)
            else:
                fz = sp.sqrt(sp.Rational((2 * l - 1) ** 2, (l + m) * (l - m)))
                fr = sp.sqrt(sp.Rational((l - 1 + m) * (l - 1 - m), (l + m) * (l - m)))
                e = fz * RZ * solid(l - 1, m) - fr * r2 * solid(l - 2, m)
        _solid[(l, m)] = sp.expand(e)
    return _solid[(l, m)].subs({RX: X, RY: Y, RZ: Z}, simultaneous=True)


def dfact(n):
    """(n - 1)!! with (-1)!! = 1, the Gaussian moment of one axis."""
    out = 1
    k = n - 1
    while k > 1:
        out *= k
        k -= 2
    return sp.Integer(out)


def integrate(e):
    """Replaces every monomial of u by its Gaussian moment over exp(-p u^2),
    in units of (pi/p)^(3/2), which introduces h = 1 / (2 p)."""
    poly = sp.Poly(sp.expand(e), UX, UY, UZ)
    out = sp.Integer(0)
    for (i, j, k), c in poly.terms():
        if i % 2 or j % 2 or k % 2:
            continue
        out += c * dfact(i) * dfact(j) * dfact(k) * HH ** ((i + j + k) // 2)
    return sp.expand(out)


def terms(la, lb, ma, mb):
    """P_k for one component, as a dict k -> polynomial in the vector components."""
    bra = solid(la, ma, UX + AL * RX, UY + AL * RY, UZ + AL * RZ)
    ket = solid(lb, mb, UX + BE * RX, UY + BE * RY, UZ + BE * RZ)
    e = integrate(sp.expand(bra * ket))
    out = {}
    poly = sp.Poly(e, AL, BE, HH)
    for (i, j, k), c in poly.terms():
        assert i == la - k and j == lb - k, \
            "unexpected power of alpha or beta for l = (%d, %d), m = (%d, %d): %s" % (la, lb, ma, mb, (i, j, k))
        out[k] = out.get(k, sp.Integer(0)) + c
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


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
    p, q = isqrt_exact(fr.numerator), isqrt_exact(fr.denominator)
    if p is not None and q is not None:
        return literal(F(p, q))
    return f"std::sqrt({literal(fr)})"


def coefficient(c):
    """A coefficient as the kernels write it.

    A single coefficient of a solid harmonic is a square root of a rational, but a
    coefficient of an expanded product is a sum of such, which need not be one. The
    square root is written exactly when the square is rational, and the coefficient
    is written as a decimal of full double precision otherwise. No approximation is
    made in deciding which: sympy is asked for the exact square."""
    sign = "-" if c < 0 else ""
    a = sp.simplify(abs(c))
    sq = sp.simplify(sp.expand(a ** 2))
    if sq.is_Rational:
        return sign, factor(F(int(sq.p), int(sq.q)))
    lit = "%.17g" % float(a.evalf(50))
    if "." not in lit and "e" not in lit and "E" not in lit:
        lit += ".0"
    return sign, lit


def expression(e):
    """A polynomial in the vector components as a C++ expression."""
    poly = sp.Poly(e, RX, RY, RZ)
    out = []
    for powers, c in sorted(poly.terms(), key=lambda t: tuple(-p for p in t[0])):
        sign, mag = coefficient(c)
        mono = " * ".join(sum(([s] * p for s, p in zip(("x", "y", "z"), powers)), []))
        if not mono:
            out.append(f"{sign}{mag}")
        elif mag == "1.0":
            out.append(f"{sign}{mono}")
        else:
            out.append(f"{sign}{mag} * {mono}")
    return " + ".join(out).replace("+ -", "- ")


def name(m):
    return "0" if m == 0 else (f"p{m}" if m > 0 else f"m{-m}")


def sizes(la, lb):
    """The number of terms the component block costs, for planning."""
    n = 0
    for ma in range(-la, la + 1):
        for mb in range(-lb, lb + 1):
            if la == lb and mb < ma:
                continue
            for k, p in terms(la, lb, ma, mb).items():
                n += len(sp.Poly(p, RX, RY, RZ).terms())
    return n


# NOTE: the angular half writes one row of the values per component and reads the
# accumulators and the vector between the atoms, so it runs out of registers well
# before all components fit one loop. The components are independent, so they are
# split into passes of at most CHUNK of them, and of at most BUDGET terms in total,
# whichever binds first. Counting components alone is not enough: a component of the
# highest orders carries hundreds of terms on its own, and four of them together
# defeat the vectorizer.
CHUNK = 4
BUDGET = 48


def slots(la, lb):
    """The components the kernel computes, and the map from every component of the
    block to the one it is copied from."""
    ncomp = (2 * la + 1) * (2 * lb + 1)
    index = lambda ma, mb: (ma + la) * (2 * lb + 1) + (mb + lb)
    computed, source = [], [0] * ncomp
    for ma in range(-la, la + 1):
        for mb in range(-lb, lb + 1):
            if la == lb and mb < ma:
                source[index(ma, mb)] = index(mb, ma)
            else:
                source[index(ma, mb)] = index(ma, mb)
                computed.append((ma, mb, index(ma, mb)))
    return computed, source


def weight(la, lb, k):
    """The C++ expression of the weight of accumulator k."""
    out = ["fbase"] + ["fal"] * (la - k) + ["fbe"] * (lb - k) + ["fh"] * k
    return " * ".join(out)


def emit(la, lb):
    tag = LETTER[la] + LETTER[lb]
    fn = tag.lower()
    nacc = min(la, lb) + 1
    ncomp = (2 * la + 1) * (2 * lb + 1)
    momenta = f"{MOMENT[la]} and {MOMENT[lb]}"
    computed, source = slots(la, lb)

    per = {}
    for ma, mb, idx in computed:
        per[(ma, mb)] = terms(la, lb, ma, mb)

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
            f"/// @note The integrals carry {nacc} terms, one for each order of the harmonic which\n"
            "/// survives the integration over the Gaussian product center. Only their\n"
            "/// prefactors depend on the pair of primitives, so the buffer holds those alone and\n"
            "/// the integrals of the angular components are formed straight into the values.\n")
    if la == lb:
        hpp += ("/// @note The two harmonics enter symmetrically, so the block of the angular\n"
                "/// components is symmetric and only its triangle is computed, the rest copied.\n")
    hpp += (f"auto compute_{fn}_overlap(double               *values,\n"
            "                        const size_t          nvalues,\n"
            "                        const CBasisFunction &bra,\n"
            "                        const CBasisFunction &ket,\n"
            "                        const CSimdMatrix    &coordinates,\n"
            "                        const double          threshold) -> void;\n\n")
    hpp += "}  // namespace simdovl\n\n"
    hpp += f"#endif /* SimdOverlapRec{tag}_hpp */\n"

    accs = "\n".join(f"    auto *pe_{k} = buffer.data({k});" for k in range(nacc))
    acc_aligned = ", ".join([f"pe_{k}" for k in range(nacc)] + ["ab_2"])
    acc_f = "\n\n".join(f"        const auto f_{k} = {weight(la, lb, k)};" for k in range(nacc))
    acc_add = "\n".join(f"            pe_{k}[k] += f_{k} * fss;" for k in range(nacc))

    def count(ma, mb):
        return sum(len(sp.Poly(v, RX, RY, RZ).terms()) for v in per[(ma, mb)].values())

    groups, cur, load = [], [], 0
    for c in computed:
        n = count(c[0], c[1])
        if cur and (len(cur) >= CHUNK or load + n > BUDGET):
            groups.append(cur)
            cur, load = [], 0
        cur.append(c)
        load += n
    if cur:
        groups.append(cur)
    ptrs = "\n".join(f"    auto *pc_{i} = values + {idx} * nvalues;"
                     for i, (ma, mb, idx) in enumerate(computed))

    loops = ""
    for gi, grp in enumerate(groups):
        first = computed.index(grp[0])
        ks = sorted({k for ma, mb, _ in grp for k in per[(ma, mb)]})
        used = set()
        for ma, mb, _ in grp:
            for p in per[(ma, mb)].values():
                used |= p.free_symbols
        vs = [(s, n) for s, n in ((RX, "x"), (RY, "y"), (RZ, "z")) if s in used]
        aligned = ", ".join([f"pe_{k}" for k in ks] + ["ab_" + n for _, n in vs])
        loops += f"#pragma omp simd aligned({aligned} : simd::cache_line_size())\n"
        loops += "    for (size_t k = 0; k < nmax; k++)\n    {\n"
        for _, n in vs:
            loops += f"        const auto {n} = ab_{n}[k];\n"
        loops += "\n" if vs else ""
        for k in ks:
            loops += f"        const auto e_{k} = pe_{k}[k];\n"
        loops += "\n"
        for i, (ma, mb, _) in enumerate(grp):
            parts = [f"e_{k} * ({expression(p)})" for k, p in sorted(per[(ma, mb)].items())]
            loops += f"        pc_{first + i}[k] = " + " + ".join(parts) + ";\n\n"
        loops = loops.rstrip("\n") + "\n    }\n\n"

    cpp = LICENSE + f'#include "SimdOverlapRec{tag}.hpp"\n\n'
    cpp += "#include <algorithm>\n#include <cmath>\n#include <cstddef>\n#include <ranges>\n#include <string>\n\n"
    cpp += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "ScreeningFunc.hpp"\n'
            '#include "SimdAlign.hpp"\n#include "SimdDimensions.hpp"\n#include "SimdPrimitives.hpp"\n\n')
    cpp += "namespace simdovl {  // simdovl namespace\n\n"
    cpp += (f"auto\ncompute_{fn}_overlap(double               *values,\n"
            "                   const size_t          nvalues,\n"
            "                   const CBasisFunction &bra,\n"
            "                   const CBasisFunction &ket,\n"
            "                   const CSimdMatrix    &coordinates,\n"
            "                   const double          threshold) -> void\n{\n")
    cpp += (f"    if ((bra.get_angular_momentum() != {la}) || (ket.get_angular_momentum() != {lb}))\n"
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
    cpp += ("    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the\n"
            "    // integrals of the angular components are formed straight into the values and\n"
            "    // are not written a second time.\n\n")
    cpp += f"    auto buffer = simdfunc::make_primitive_buffer(dimensions, {nacc});\n\n"
    cpp += ("    if (buffer.number_of_columns() == 0)\n    {\n"
            f"        std::fill(values, values + {ncomp} * nvalues, 0.0);\n\n        return;\n    }}\n\n")
    cpp += "    const auto nmax = buffer.number_of_columns();\n\n"
    cpp += accs + "\n\n"
    cpp += ("    // NOTE: the components of the vector between the atoms and its squared length\n"
            "    // are carried by the coordinates, so the angular half below reads rows which\n"
            "    // are already in place.\n\n")
    cpp += ("    const auto *ab_x = coordinates.data(6);\n    const auto *ab_y = coordinates.data(7);\n"
            "    const auto *ab_z = coordinates.data(8);\n\n    const auto *ab_2 = coordinates.data(9);\n\n")
    cpp += "    constexpr auto fpi = mathconst::pi_value();\n\n"
    cpp += "    // accumulate the prefactor of each term over the pairs of primitives\n\n"
    cpp += ("    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {\n"
            "        const auto ncols = pair.ncols;\n\n"
            "        const auto fexp = pair.aexp + pair.bexp;\n\n"
            "        const auto fmu = pair.aexp * pair.bexp / fexp;\n\n"
            "        const auto fovl = fpi / fexp;\n\n"
            "        const auto fbase = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl);\n\n"
            "        // NOTE: the Gaussian product center is displaced from the atom on bra side\n"
            "        // by fal times the vector between the atoms and from the atom on ket side by\n"
            "        // fbe times it, and fh is the second moment the integration over that center\n"
            "        // leaves behind.\n\n"
            "        const auto fal = -pair.bexp / fexp;\n\n"
            "        const auto fbe = pair.aexp / fexp;\n\n"
            "        const auto fh = 0.5 / fexp;\n\n"
            + acc_f + "\n\n"
            "        // NOTE: the exponential depends on the pair of primitives alone, so it is\n"
            "        // evaluated once and shared by the prefactors of all terms.\n\n"
            f"#pragma omp simd aligned({acc_aligned} : simd::cache_line_size())\n"
            "        for (size_t k = 0; k < ncols; k++)\n        {\n"
            "            const auto fss = std::exp(-fmu * ab_2[k]);\n\n"
            + acc_add + "\n        }\n    });\n\n")
    cpp += ("    // NOTE: the rows of the values are not aligned, as they start at the offset of\n"
            "    // this combination of basis functions in the values block, so they are kept out\n"
            "    // of the aligned clauses below.\n\n")
    cpp += ptrs + "\n\n"
    if len(groups) > 1:
        cpp += (f"    // NOTE: the components are formed in {len(groups)} loops, as the vectorizer runs out\n"
                "    // of registers with all of them in one. Only the prefactors and the vector\n"
                "    // between the atoms are loaded by more than one loop.\n\n")
    cpp += loops
    if la == lb:
        cpp += ("    // NOTE: the values of a combination of angular components are stored as one\n"
                "    // row of nvalues columns, with the component on bra side running slowest. The\n"
                "    // rows which the symmetry relates to an already formed one are copied from it,\n"
                "    // and the atom pairs beyond the reach of every pair of primitives are set to\n"
                "    // zero.\n\n")
        cpp += f"    const size_t sources[{ncomp}] = {{{', '.join(str(s) for s in source)}}};\n\n"
        cpp += (f"    for (size_t m = 0; m < {ncomp}; m++)\n    {{\n"
                "        auto *pv = values + m * nvalues;\n\n"
                "        const auto *pc = values + sources[m] * nvalues;\n\n"
                "        if (pv != pc) std::copy(pc, pc + nmax, pv);\n\n"
                "        std::fill(pv + nmax, pv + nvalues, 0.0);\n    }\n")
    else:
        cpp += ("    // NOTE: the atom pairs beyond the reach of every pair of primitives have no\n"
                "    // contribution and are set to zero.\n\n")
        cpp += (f"    for (size_t m = 0; m < {ncomp}; m++)\n    {{\n"
                "        auto *pv = values + m * nvalues;\n\n"
                "        std::fill(pv + nmax, pv + nvalues, 0.0);\n    }\n")
    cpp += "}\n\n}  // namespace simdovl\n"

    io.open(SRC + f"SimdOverlapRec{tag}.hpp", "w").write(hpp)
    io.open(SRC + f"SimdOverlapRec{tag}.cpp", "w").write(cpp)
    return tag


if __name__ == "__main__":
    ls = [int(a) for a in sys.argv[1:]] or [1, 2, 3]
    for la in ls:
        for lb in ls:
            print("wrote", emit(la, lb), flush=True)


def emit_single(l):
    """Emits the kernel of one S type function and one function of angular momentum
    l, covering both orders.

    The angular half does not depend on the order: the harmonic is the same
    polynomial of the vector between the atoms either way. Only the prefactor
    differs, carrying (a / p) raised to the power l when the harmonic sits on the
    ket side and -(b / p) raised to it when it sits on the bra side, so the order is
    selected once per pair of primitives and not inside any loop.
    """
    tag = "SL" + LETTER[l]
    fn = tag.lower()
    ncomp = 2 * l + 1
    ms = list(range(-l, l + 1))

    per = {m: terms(0, l, 0, m) for m in ms}
    for m in ms:
        assert set(per[m]) == {0}, "a single S type function must leave one term"

    hpp = LIC_HPP + f"#ifndef SimdOverlapRec{tag}_hpp\n#define SimdOverlapRec{tag}_hpp\n\n"
    hpp += "#include <cstddef>\n\n"
    hpp += '#include "BasisFunction.hpp"\n#include "SimdMatrix.hpp"\n\n'
    hpp += "namespace simdovl {  // simdovl namespace\n\n"
    hpp += ("/// @brief Computes the overlap integrals of a combination of one basis function\n"
            f"/// of zero angular momentum and one of angular momentum {MOMENT[l]}, in either order.\n"
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
            "/// @note The angular half is the same for both orders, as the harmonic is the same\n"
            "/// polynomial of the vector between the atoms either way. The orders differ only in\n"
            "/// the prefactor, which is selected once for the whole combination.\n"
            "/// @note One term survives the integration over the Gaussian product center, so the\n"
            "/// buffer holds a single accumulator and the integrals of the angular components are\n"
            "/// formed straight into the values.\n")
    hpp += (f"auto compute_{fn}_overlap(double               *values,\n"
            "                         const size_t          nvalues,\n"
            "                         const CBasisFunction &bra,\n"
            "                         const CBasisFunction &ket,\n"
            "                         const CSimdMatrix    &coordinates,\n"
            "                         const double          threshold) -> void;\n\n")
    hpp += "}  // namespace simdovl\n\n"
    hpp += f"#endif /* SimdOverlapRec{tag}_hpp */\n"

    computed = [(m, i) for i, m in enumerate(ms)]
    groups, cur, load = [], [], 0
    for m, idx in computed:
        n = len(sp.Poly(per[m][0], RX, RY, RZ).terms())
        if cur and (len(cur) >= CHUNK or load + n > BUDGET):
            groups.append(cur)
            cur, load = [], 0
        cur.append((m, idx))
        load += n
    if cur:
        groups.append(cur)

    ptrs = "\n".join(f"    auto *pc_{i} = values + {i} * nvalues;" for _, i in computed)

    loops = ""
    for grp in groups:
        used = set()
        for m, _ in grp:
            used |= per[m][0].free_symbols
        vs = [(s, n) for s, n in ((RX, "x"), (RY, "y"), (RZ, "z")) if s in used]
        aligned = ", ".join(["pe_0"] + ["ab_" + n for _, n in vs])
        loops += f"#pragma omp simd aligned({aligned} : simd::cache_line_size())\n"
        loops += "    for (size_t k = 0; k < nmax; k++)\n    {\n"
        for _, n in vs:
            loops += f"        const auto {n} = ab_{n}[k];\n"
        loops += "\n        const auto e_0 = pe_0[k];\n\n"
        for m, i in grp:
            loops += f"        pc_{i}[k] = e_0 * ({expression(per[m][0])});\n\n"
        loops = loops.rstrip("\n") + "\n    }\n\n"

    cpp = LICENSE + f'#include "SimdOverlapRec{tag}.hpp"\n\n'
    cpp += "#include <algorithm>\n#include <cmath>\n#include <cstddef>\n#include <ranges>\n#include <string>\n\n"
    cpp += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "ScreeningFunc.hpp"\n'
            '#include "SimdAlign.hpp"\n#include "SimdDimensions.hpp"\n#include "SimdPrimitives.hpp"\n\n')
    cpp += "namespace simdovl {  // simdovl namespace\n\n"
    cpp += (f"auto\ncompute_{fn}_overlap(double               *values,\n"
            "                    const size_t          nvalues,\n"
            "                    const CBasisFunction &bra,\n"
            "                    const CBasisFunction &ket,\n"
            "                    const CSimdMatrix    &coordinates,\n"
            "                    const double          threshold) -> void\n{\n")
    cpp += (f"    const auto lbra = bra.get_angular_momentum();\n\n"
            f"    const auto lket = ket.get_angular_momentum();\n\n"
            f"    if (!(((lbra == 0) && (lket == {l})) || ((lbra == {l}) && (lket == 0))))\n"
            "    {\n        errors::assertMsgCritical(\n"
            f'            false, std::string("SimdOverlapRec{tag}.compute_{fn}_overlap: Basis functions must be of angular momenta zero and {MOMENT[l]}"));\n'
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
    cpp += ("    // NOTE: the buffer holds the contracted prefactor alone, as the harmonic\n"
            "    // factors out of the sum over the pairs of primitives and the integrals of the\n"
            "    // angular components are formed straight into the values.\n\n")
    cpp += "    auto buffer = simdfunc::make_primitive_buffer(dimensions, 1);\n\n"
    cpp += ("    if (buffer.number_of_columns() == 0)\n    {\n"
            f"        std::fill(values, values + {ncomp} * nvalues, 0.0);\n\n        return;\n    }}\n\n")
    cpp += "    const auto nmax = buffer.number_of_columns();\n\n    auto *pe_0 = buffer.data(0);\n\n"
    cpp += ("    // NOTE: the components of the vector between the atoms and its squared length\n"
            "    // are carried by the coordinates, so the angular half below reads rows which\n"
            "    // are already in place.\n\n")
    cpp += ("    const auto *ab_x = coordinates.data(6);\n    const auto *ab_y = coordinates.data(7);\n"
            "    const auto *ab_z = coordinates.data(8);\n\n    const auto *ab_2 = coordinates.data(9);\n\n")
    cpp += "    constexpr auto fpi = mathconst::pi_value();\n\n"
    cpp += ("    // NOTE: the harmonic sits on whichever side carries the angular momentum, and\n"
            "    // the Gaussian product center is displaced from it by (a / p) times the vector\n"
            "    // between the atoms when that is the ket side and by -(b / p) when it is the bra\n"
            "    // side. The order is therefore settled once here and not inside any loop.\n\n")
    cpp += "    const auto on_ket = (lbra == 0);\n\n"
    cpp += "    // accumulate the prefactor of each pair of primitives\n\n"
    cpp += ("    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {\n"
            "        const auto ncols = pair.ncols;\n\n"
            "        const auto fexp = pair.aexp + pair.bexp;\n\n"
            "        const auto fmu = pair.aexp * pair.bexp / fexp;\n\n"
            "        const auto fovl = fpi / fexp;\n\n"
            "        const auto fbase = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl);\n\n"
            "        const auto fr = on_ket ? (pair.aexp / fexp) : (-pair.bexp / fexp);\n\n"
            "        const auto ffact = fbase" + " * fr" * l + ";\n\n"
            "#pragma omp simd aligned(pe_0, ab_2 : simd::cache_line_size())\n"
            "        for (size_t k = 0; k < ncols; k++)\n        {\n"
            "            pe_0[k] += ffact * std::exp(-fmu * ab_2[k]);\n        }\n    });\n\n")
    cpp += ("    // NOTE: the rows of the values are not aligned, as they start at the offset of\n"
            "    // this combination of basis functions in the values block, so they are kept out\n"
            "    // of the aligned clauses below.\n\n")
    cpp += ptrs + "\n\n"
    if len(groups) > 1:
        cpp += (f"    // NOTE: the components are formed in {len(groups)} loops, as the vectorizer runs out\n"
                "    // of registers with all of them in one. Only the prefactor and the vector\n"
                "    // between the atoms are loaded by more than one loop.\n\n")
    cpp += loops
    cpp += ("    // NOTE: the atom pairs beyond the reach of every pair of primitives have no\n"
            "    // contribution and are set to zero.\n\n")
    cpp += (f"    for (size_t m = 0; m < {ncomp}; m++)\n    {{\n"
            "        auto *pv = values + m * nvalues;\n\n"
            "        std::fill(pv + nmax, pv + nvalues, 0.0);\n    }\n")
    cpp += "}\n\n}  // namespace simdovl\n"

    io.open(SRC + f"SimdOverlapRec{tag}.hpp", "w").write(hpp)
    io.open(SRC + f"SimdOverlapRec{tag}.cpp", "w").write(cpp)
    return tag
