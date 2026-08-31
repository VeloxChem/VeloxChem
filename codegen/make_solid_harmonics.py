"""Emits the solid harmonics recursion of one angular momentum, in the shape of
the generated files which already exist for the orders below it."""
import io
import os
from fractions import Fraction as F

# NOTE: the root of the checkout is taken from the location of this file, so the
# generators run wherever the repository is placed.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src", "simd_func") + "/"
LICENSE = io.open(SRC + 'SimdHarmonicsRecQ.cpp').read().split('#include "SimdHarmonicsRecQ.hpp"')[0]
LIC_HPP = io.open(SRC + 'SimdHarmonicsRecQ.hpp').read().split('#ifndef')[0]

LETTER = {1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i', 7: 'k', 8: 'l',
          9: 'm', 10: 'n', 11: 'o', 12: 'q', 13: 'r', 14: 't'}
WORD = {10: 'ten', 11: 'eleven', 12: 'twelve', 13: 'thirteen', 14: 'fourteen'}


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
    """A fraction as the files write it."""
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


def emit(l):
    prev, prev2 = LETTER[l - 1], LETTER[l - 2]
    me = LETTER[l]
    nrows = 2 * l + 1
    fact = F(2 * l - 1, 2 * l)
    fz_edge = F(2 * l - 1)

    # the loops, as the order below splits them: the extremes and the negative
    # orders, then the middle, then the positive orders
    loop1 = [l, -l, -(l - 1)] + list(range(-(l - 2), -3))
    loop2 = list(range(-3, 5))
    loop3 = list(range(5, l))
    loops = [loop1, loop2, loop3]
    nloops = len(loops)

    def gen_factors():
        out = [f"    const auto fact = {factor(fact)};",
               f"    const auto fz_m{l-1} = {factor(fz_edge)};"]
        for m in range(-(l - 2), l - 1):
            fz = F((2 * l - 1) ** 2, (l + m) * (l - m))
            fr = F((l - 1 + m) * (l - 1 - m), (l + m) * (l - m))
            out.append(f"    const auto fz_{name(m)} = {factor(fz)};")
            out.append(f"    const auto fr_{name(m)} = {factor(fr)};")
        out.append(f"    const auto fz_p{l-1} = {factor(fz_edge)};")
        return "\n".join(out)

    def body(ms, first):
        """One vectorized loop over the atom pairs."""
        needs_xy = (l in ms) or (-l in ms)
        t_in, u_in = [], []
        for m in ms:
            if abs(m) == l:
                t_in += [l - 1, -(l - 1)]
            elif abs(m) == l - 1:
                t_in += [m]
            else:
                t_in += [m]; u_in += [m]
        t_in = sorted(set(t_in), key=lambda v: (v, ))
        u_in = sorted(set(u_in))
        aligned = (["pr_x", "pr_y"] if needs_xy else []) + ["pr_z", "pr_2"]
        aligned += [f"pt_{name(m)}" for m in t_in] + [f"pu_{name(m)}" for m in u_in]
        aligned += [f"ps_{name(m)}" for m in sorted(ms)]
        s = ""
        if not first:
            s += "\n"
        s += ("    // NOTE: the rows are formed in three loops, as the vectorizer runs out of\n"
              f"    // registers with all {nrows} of them in one. Only the components and the squared\n"
              "    // distance are loaded by more than one loop, every other value once.\n\n")
        s += f"#pragma omp simd aligned({', '.join(aligned)} : simd::cache_line_size())\n"
        s += "    for (size_t i = 0; i < ncols; i++)\n    {\n"
        if needs_xy:
            s += "        const auto x = pr_x[i];\n        const auto y = pr_y[i];\n"
        s += "        const auto z = pr_z[i];\n\n        const auto r_2 = pr_2[i];\n\n"
        for m in t_in:
            s += f"        const auto t_{name(m)} = pt_{name(m)}[i];\n"
        if u_in:
            s += "\n"
        for m in u_in:
            s += f"        const auto u_{name(m)} = pu_{name(m)}[i];\n"
        s += "\n"
        for m in ms:
            if m == l:
                s += f"        ps_p{l}[i] = fact * (x * t_p{l-1} - y * t_m{l-1});\n\n"
            elif m == -l:
                s += f"        ps_m{l}[i] = fact * (y * t_p{l-1} + x * t_m{l-1});\n\n"
            elif abs(m) == l - 1:
                s += f"        ps_{name(m)}[i] = fz_{name(m)} * z * t_{name(m)};\n\n"
            else:
                s += (f"        ps_{name(m)}[i] = fz_{name(m)} * z * t_{name(m)}"
                      f" - fr_{name(m)} * r_2 * u_{name(m)};\n\n")
        return s.rstrip("\n") + "\n    }\n"

    args = (f"const CSimdMatrix &{prev}_harmonics, const CSimdMatrix &{prev2}_harmonics, "
            f"const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates")
    sig = f"make_{me}_solid_harmonics({args}) -> CSimdMatrix"

    cpp = LICENSE + f'#include "SimdHarmonicsRec{me.upper()}.hpp"\n\n'
    cpp += "#include <cmath>\n#include <cstddef>\n#include <string>\n\n"
    cpp += '#include "ErrorHandler.hpp"\n#include "SimdAlign.hpp"\n\n'
    cpp += "namespace simdfunc {  // simdfunc namespace\n\nauto\n" + sig + "\n{\n"
    cpp += (f'    errors::assertMsgCritical({prev}_harmonics.number_of_rows() == {2*(l-1)+1},\n'
            f'                              std::string("SimdHarmonics.make_{me}_solid_harmonics: Harmonics of angular momentum {WORD[l-1]} must have {2*(l-1)+1} rows"));\n\n')
    cpp += (f'    errors::assertMsgCritical({prev2}_harmonics.number_of_rows() == {2*(l-2)+1},\n'
            f'                              std::string("SimdHarmonics.make_{me}_solid_harmonics: Harmonics of angular momentum {WORD[l-2]} must have {2*(l-2)+1} rows"));\n\n')
    cpp += ('    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,\n'
            f'                              std::string("SimdHarmonics.make_{me}_solid_harmonics: Harmonics of angular momentum one must have three rows"));\n\n')
    cpp += ('    errors::assertMsgCritical(coordinates.number_of_rows() == 7,\n'
            f'                              std::string("SimdHarmonics.make_{me}_solid_harmonics: Coordinates must have seven rows"));\n\n')
    cpp += (f'    errors::assertMsgCritical({prev}_harmonics.number_of_columns() == coordinates.number_of_columns(),\n'
            f'                              std::string("SimdHarmonics.make_{me}_solid_harmonics: Harmonics and coordinates must have the same number of columns"));\n\n')
    cpp += "    const auto ncols = coordinates.number_of_columns();\n\n"
    cpp += f"    auto matrix = CSimdMatrix({nrows}, ncols);\n\n"
    cpp += ("    // NOTE: the factors of the recursion depend on the angular momentum alone, so\n"
            "    // they are formed once for the whole matrix instead of once for every atom pair.\n\n")
    cpp += gen_factors() + "\n\n"
    cpp += ("    // NOTE: the pointers to the rows are taken once, as the accessor of a row is\n"
            "    // bounds checked and would otherwise be called for every atom pair.\n\n")
    cpp += ("    const auto *pr_y = p_harmonics.data(0);\n    const auto *pr_z = p_harmonics.data(1);\n"
            "    const auto *pr_x = p_harmonics.data(2);\n\n    const auto *pr_2 = coordinates.data(6);\n\n")
    for i, m in enumerate(range(-(l - 1), l)):
        cpp += f"    const auto *pt_{name(m)} = {prev}_harmonics.data({i});\n"
    cpp += "\n"
    for i, m in enumerate(range(-(l - 2), l - 1)):
        cpp += f"    const auto *pu_{name(m)} = {prev2}_harmonics.data({i});\n"
    cpp += "\n"
    for i, m in enumerate(range(-l, l + 1)):
        cpp += f"    auto *ps_{name(m)} = matrix.data({i});\n"
    cpp += "\n"
    if nloops == 1:
        cpp += ("    // NOTE: every value entering the recursion is loaded before any result is\n"
                "    // stored, so that the rows of the previous orders are read once and the\n"
                "    // stores cannot force the loads to be repeated.\n\n")
    for n, ms in enumerate(loops):
        cpp += body(ms, n == 0)
    cpp += "\n    return matrix;\n}\n\n}  // namespace simdfunc\n"

    hpp = LIC_HPP + f"#ifndef SimdHarmonicsRec{me.upper()}_hpp\n#define SimdHarmonicsRec{me.upper()}_hpp\n\n"
    hpp += '#include "SimdMatrix.hpp"\n\nnamespace simdfunc {  // simdfunc namespace\n\n'
    hpp += (f"/// @brief Creates the real solid harmonics of angular momentum {WORD[l]} of the\n"
            "/// vectors between the atoms of the atom pairs.\n"
            f"/// @param {prev}_harmonics The solid harmonics of angular momentum {WORD[l-1]}.\n"
            f"/// @param {prev2}_harmonics The solid harmonics of angular momentum {WORD[l-2]}.\n"
            "/// @param p_harmonics The solid harmonics of angular momentum one, which supply\n"
            "/// the components of the vector between the atoms.\n"
            "/// @param coordinates The coordinates of the atom pairs, as seven rows, whose last\n"
            "/// row supplies the squared distance of the atom pair.\n"
            f"/// @return The matrix of {nrows} rows and as many columns as the coordinates, with\n"
            f"/// the row of index m + {l} holding the solid harmonic of order m.\n"
            "/// @note The harmonics follow Eqs. (A4) to (A6) of the supporting information of\n"
            "/// J. Chem. Theory Comput. 2020, 16, 2570.\n")
    hpp += f"auto {sig};\n\n}}  // namespace simdfunc\n\n#endif /* SimdHarmonicsRec{me.upper()}_hpp */\n"
    return cpp, hpp


for l in (13, 14):
    cpp, hpp = emit(l)
    u = LETTER[l].upper()
    io.open(SRC + f'SimdHarmonicsRec{u}.cpp', 'w').write(cpp)
    io.open(SRC + f'SimdHarmonicsRec{u}.hpp', 'w').write(hpp)
    print(f"wrote SimdHarmonicsRec{u}.cpp ({len(cpp.splitlines())} lines) and .hpp")
