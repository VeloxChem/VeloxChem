"""Emits the three-center electron repulsion kernels of two S type functions on
the a and b sides and one of angular momentum l on the c side, for l of two to
six.

    (ss|J|l)_m = (ss|J|s)^(l) sum_{l1} (alpha^l1 p^l2 / q^l)
                              sum_{m1,m2} C^{l,m}_{l1 m1, l2 m2} S_{l1,m1}(AB) S_{l2,m2}(BC)

with l2 = l - l1 and the auxiliary integral of order n the integral of three S
type functions with the Boys function of order n,

    (ss|J|s)^(n) = 2 pi^(5/2) N_a N_b N_c exp(-mu AB^2) F_n(rho PC^2) / (p gamma sqrt(q))

The coefficients come from make_addition_table.py; see the section on it in
README.md, and set VLX_HARM_PROBE so that its convention check runs.

Run as `python codegen/make_coulomb_addition_kernels.py <l> [<l> ...]`.
"""
import io
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from make_addition_table import table  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src", "simd_t3c_electron_repulsion") + "/"

LICENSE = io.open(SRC + 'SimdThreeCenterElectronRepulsionRecSSS.cpp').read().split('#include')[0]
LIC_HPP = io.open(SRC + 'SimdThreeCenterElectronRepulsionRecSSS.hpp').read().split('#ifndef')[0]

LETTER = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G', 5: 'H', 6: 'I'}
WORD = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six'}


def name(m):
    return "0" if m == 0 else (f"p{m}" if m > 0 else f"m{-m}")


def _isqrt_exact(n):
    r = int(n ** 0.5)
    for c in (r - 1, r, r + 1):
        if c >= 0 and c * c == n:
            return c
    return None


def _terminating(fr):
    d = fr.q
    while d % 2 == 0:
        d //= 2
    return d == 1


def _literal(fr):
    if _terminating(fr):
        s = repr(float(fr))
        return s if '.' in s or 'e' in s else s + '.0'
    return f"{fr.p}.0 / {fr.q}.0"


def coefficient(c):
    """A coefficient of the table as the kernels write it. Every one of them is a
    signed square root of a rational, so it is exact whenever that root is."""
    c = sp.nsimplify(sp.simplify(c))
    sign = '-' if c < 0 else ''
    sq = sp.Rational(sp.simplify(c ** 2))
    p, q = _isqrt_exact(sq.p), _isqrt_exact(sq.q)
    if p is not None and q is not None:
        return sign + _literal(sp.Rational(p, q))
    return sign + f"std::sqrt({_literal(sq)})"


def _terms(coeffs, m, rows):
    """The bracket of one angular component, as a sum over the products."""
    out = []
    for key in rows:
        c = coeffs.get(key, {}).get(m)
        if c is None:
            continue
        s = coefficient(c)
        var = rows[key]
        if s == '1.0':
            out.append(f"{var}[k]")
        elif s == '-1.0':
            out.append(f"-{var}[k]")
        else:
            out.append(f"{s} * {var}[k]")
    return out


def emit(l):
    me = LETTER[l]
    ncomps = 2 * l + 1
    coeffs = table(l)

    hpp = LIC_HPP + f"\n#ifndef SimdThreeCenterElectronRepulsionRecSS{me}_hpp\n"
    hpp += f"#define SimdThreeCenterElectronRepulsionRecSS{me}_hpp\n\n"
    hpp += "#include <cstddef>\n#include <vector>\n\n"
    hpp += '#include "BasisFunction.hpp"\n#include "SimdMatrix.hpp"\n\n'
    hpp += "namespace simdt3ceri {  // simdt3ceri namespace\n\n"
    hpp += ("/// @brief Computes the three-center electron repulsion integrals of two basis\n"
            "/// functions of zero angular momentum on a and b sides and one of angular\n"
            f"/// momentum {WORD[l]} on c side, over the atom pairs of a block and for one atom on\n"
            "/// c side.\n"
            "/// @param values The values of the combination of basis functions, whose slices\n"
            "/// of the atom on c side this kernel writes.\n"
            "/// @param npairs The number of surviving atom pairs of the combination.\n"
            "/// @param natoms The number of atoms on c side of the block.\n"
            "/// @param iatom The index of the atom on c side among the atoms of the block.\n"
            "/// @param a_function The basis function on a side.\n"
            "/// @param b_function The basis function on b side.\n"
            "/// @param c_function The basis function on c side.\n"
            "/// @param ab_harmonics The solid harmonics of the vectors between the atoms of\n"
            f"/// the atom pairs, reaching at least angular momentum {WORD[l]}.\n"
            "/// @param bc_harmonics The solid harmonics of the vectors from the atoms on b\n"
            f"/// side to the atom on c side, reaching at least angular momentum {WORD[l]}.\n"
            "/// @param ab_coordinates The coordinates of the atom pairs, whose rows zero to\n"
            "/// two carry the atoms on a side and whose row six carries their squared\n"
            "/// distances.\n"
            "/// @param bc_coordinates The coordinates of the atoms on b side in rows zero to\n"
            "/// two and of the atom on c side in rows three to five.\n"
            "/// @param threshold The screening threshold of the integrals.\n"
            f"/// @note The addition theorem splits the harmonic of the vector from the product\n"
            f"/// center of an atom pair to the atom on c side into {l + 1} bidegrees. The first and\n"
            "/// the last are the harmonics of one side alone, and the ones between them couple\n"
            "/// the orders of both, with the coefficients make_addition_table.py produces.\n")
    hpp += f"auto compute_ss{me.lower()}_electron_repulsion(double                         *values,\n"
    pad = " " * len(f"auto compute_ss{me.lower()}_electron_repulsion(")
    for a in ["const size_t                    npairs,", "const size_t                    natoms,",
              "const size_t                    iatom,", "const CBasisFunction           &a_function,",
              "const CBasisFunction           &b_function,", "const CBasisFunction           &c_function,",
              "const std::vector<CSimdMatrix> &ab_harmonics,", "const std::vector<CSimdMatrix> &bc_harmonics,",
              "const CSimdMatrix              &ab_coordinates,", "const CSimdMatrix              &bc_coordinates,",
              "const double                    threshold) -> void;"]:
        hpp += pad + a + "\n"
    hpp += f"\n}}  // namespace simdt3ceri\n\n#endif /* SimdThreeCenterElectronRepulsionRecSS{me}_hpp */\n"

    c = LICENSE + f'#include "SimdThreeCenterElectronRepulsionRecSS{me}.hpp"\n\n'
    c += "#include <algorithm>\n#include <cmath>\n#include <ranges>\n#include <string>\n#include <vector>\n\n"
    c += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "ScreeningFunc.hpp"\n'
          '#include "SimdAlign.hpp"\n#include "SimdBoysFunc.hpp"\n#include "SimdDimensions.hpp"\n'
          '#include "SimdVariableMatrix.hpp"\n\n')
    c += "namespace simdt3ceri {  // simdt3ceri namespace\n\nauto\n"
    c += f"compute_ss{me.lower()}_electron_repulsion(double                         *values,\n"
    pad = " " * len(f"compute_ss{me.lower()}_electron_repulsion(")
    for a in ["const size_t                    npairs,", "const size_t                    natoms,",
              "const size_t                    iatom,", "const CBasisFunction           &a_function,",
              "const CBasisFunction           &b_function,", "const CBasisFunction           &c_function,",
              "const std::vector<CSimdMatrix> &ab_harmonics,", "const std::vector<CSimdMatrix> &bc_harmonics,",
              "const CSimdMatrix              &ab_coordinates,", "const CSimdMatrix              &bc_coordinates,",
              "const double                    threshold) -> void"]:
        c += pad + a + "\n"
    c += "{\n"
    n = f"SimdThreeCenterElectronRepulsionRecSS{me}.compute_ss{me.lower()}_electron_repulsion"
    c += (f"    if ((a_function.get_angular_momentum() != 0) || (b_function.get_angular_momentum() != 0) ||\n"
          f"        (c_function.get_angular_momentum() != {l}))\n    {{\n"
          f"        errors::assertMsgCritical(\n            false,\n"
          f'            std::string("{n}: Basis functions must be of angular momenta zero, zero and {WORD[l]}"));\n    }}\n\n')
    c += (f"    if ((ab_harmonics.size() < {l}) || (bc_harmonics.size() < {l}))\n    {{\n"
          f"        errors::assertMsgCritical(\n            false, std::string(\"{n}: Harmonics must reach angular momentum {WORD[l]}\"));\n    }}\n\n")
    c += (f"    if (npairs > ab_coordinates.number_of_columns())\n    {{\n"
          f"        errors::assertMsgCritical(\n            false, std::string(\"{n}: Number of atom pairs exceeds coordinates\"));\n    }}\n\n")
    c += (f"    if (iatom >= natoms)\n    {{\n"
          f"        errors::assertMsgCritical(\n            false, std::string(\"{n}: Index of atom on c side is out of range\"));\n    }}\n\n")
    c += "    if (npairs == 0) return;\n\n"
    for s in ["a_exps = a_function.exponents()", "b_exps = b_function.exponents()", "c_exps = c_function.exponents()",
              "a_norms = a_function.normalization_factors()", "b_norms = b_function.normalization_factors()",
              "c_norms = c_function.normalization_factors()"]:
        c += f"    const auto &{s};\n\n"
    c += ("    const auto nprim_a = a_exps.size();\n\n    const auto nprim_b = b_exps.size();\n\n"
          "    const auto nprim_c = c_exps.size();\n\n    const auto nprims = nprim_a * nprim_b * nprim_c;\n\n")
    c += ("    // NOTE: the triples of primitives are screened with the threshold of the\n"
          "    // integrals divided by their number, as their contributions accumulate into\n"
          "    // a single value and the error of the sum is bounded by the number of terms.\n\n")
    c += ("    const auto dimensions = simdfunc::make_column_dimensions(a_function,\n"
          "                                                             b_function,\n"
          "                                                             c_function,\n"
          "                                                             npairs,\n"
          "                                                             ab_coordinates,\n"
          "                                                             screenfunc::three_center_electron_repulsion_primitive_bound,\n"
          "                                                             threshold / static_cast<double>(nprims));\n\n")
    c += "    const auto nmax = *std::ranges::max_element(dimensions);\n\n"
    c += ("    // NOTE: the values of one atom on c side are contiguous over the atom pairs,\n"
          "    // the atoms are npairs apart and the angular components of the atoms on c\n"
          "    // side are npairs times the number of those atoms apart, as the sparsity\n"
          "    // pattern lays them out.\n\n")
    c += "    const auto stride = natoms * npairs;\n\n"
    c += f"    double *slices[{ncomps}];\n\n"
    c += f"    for (size_t m = 0; m < {ncomps}; m++) slices[m] = values + m * stride + iatom * npairs;\n\n"
    c += ("    if (nmax == 0)\n    {\n"
          f"        for (size_t m = 0; m < {ncomps}; m++) std::fill(slices[m], slices[m] + npairs, 0.0);\n\n"
          "        return;\n    }\n\n")
    c += ("    const auto *a_x = ab_coordinates.data(0);\n    const auto *a_y = ab_coordinates.data(1);\n"
          "    const auto *a_z = ab_coordinates.data(2);\n\n    const auto *ab_2 = ab_coordinates.data(6);\n\n"
          "    const auto *b_x = bc_coordinates.data(0);\n    const auto *b_y = bc_coordinates.data(1);\n"
          "    const auto *b_z = bc_coordinates.data(2);\n\n"
          "    const auto *c_x = bc_coordinates.data(3);\n    const auto *c_y = bc_coordinates.data(4);\n"
          "    const auto *c_z = bc_coordinates.data(5);\n\n")
    c += ("    auto factors = CSimdMatrix(2, nmax);\n\n    auto *e_ab = factors.data(0);\n\n"
          "    auto *pc_2 = factors.data(1);\n\n")
    c += ("    // NOTE: one row accumulates for each bidegree of the addition theorem, as\n"
          "    // they carry different powers of the exponents and cannot share an\n"
          "    // accumulator. The row of index l1 multiplies the harmonics of degree l1 of\n"
          f"    // the atom pairs and of degree {l} less l1 of the atoms on c side.\n\n")
    c += f"    auto buffer = CSimdMatrix({l + 1}, nmax);\n\n    buffer.zero();\n\n"
    for l1 in range(l + 1):
        c += f"    auto *acc_{l1} = buffer.data({l1});\n\n"
    c += ("    constexpr auto fpi = mathconst::pi_value();\n\n"
          "    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);\n\n")
    c += ("    for (size_t i = 0; i < nprim_a; i++)\n    {\n        const auto aexp = a_exps[i];\n\n"
          "        const auto anorm = a_norms[i];\n\n"
          "        for (size_t j = 0; j < nprim_b; j++)\n        {\n            const auto bexp = b_exps[j];\n\n")
    c += ("            // NOTE: the widest of the triples of this pair of primitives is\n"
          "            // searched for rather than assumed to be the last, as the bound of a\n"
          "            // triple carries its prefactor as well as its decay.\n\n"
          "            const auto first = dimensions.begin() + static_cast<long>((i * nprim_b + j) * nprim_c);\n\n"
          "            const auto npair_max = *std::ranges::max_element(first, first + static_cast<long>(nprim_c));\n\n"
          "            if (npair_max == 0) continue;\n\n"
          "            const auto pexp = aexp + bexp;\n\n            const auto fmu = aexp * bexp / pexp;\n\n"
          "            const auto frp = 1.0 / pexp;\n\n")
    c += ("#pragma omp simd aligned(e_ab, pc_2, ab_2, a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z : simd::cache_line_size())\n"
          "            for (size_t k = 0; k < npair_max; k++)\n            {\n"
          "                e_ab[k] = std::exp(-fmu * ab_2[k]);\n\n"
          "                const auto p_x = frp * (aexp * (a_x[k] - c_x[k]) + bexp * (b_x[k] - c_x[k]));\n\n"
          "                const auto p_y = frp * (aexp * (a_y[k] - c_y[k]) + bexp * (b_y[k] - c_y[k]));\n\n"
          "                const auto p_z = frp * (aexp * (a_z[k] - c_z[k]) + bexp * (b_z[k] - c_z[k]));\n\n"
          "                pc_2[k] = p_x * p_x + p_y * p_y + p_z * p_z;\n            }\n\n")
    c += (f"            // NOTE: the Boys function of every primitive on c side of this pair\n"
          f"            // is computed by one call, which fills the orders zero to {WORD[l]} of\n"
          f"            // every row. The integrals need the order {WORD[l]} alone, and the lower\n"
          f"            // orders are formed on the way to it by the recursion.\n\n")
    c += (f"            auto boys = CSimdVariableMatrix(std::vector<size_t>(first, first + static_cast<long>(nprim_c)), {l + 2});\n\n")
    c += ("            for (size_t k = 0; k < nprim_c; k++)\n            {\n"
          "                const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];\n\n"
          "                if (ncols == 0) continue;\n\n"
          "                const auto frho = pexp * c_exps[k] / (pexp + c_exps[k]);\n\n"
          "                auto *bargs = boys.data(0, k);\n\n"
          "#pragma omp simd aligned(bargs, pc_2 : simd::cache_line_size())\n"
          "                for (size_t l = 0; l < ncols; l++)\n                {\n"
          "                    bargs[l] = frho * pc_2[l];\n                }\n            }\n\n"
          "            simdfunc::compute_boys_function(boys);\n\n")
    c += ("            for (size_t k = 0; k < nprim_c; k++)\n            {\n"
          "                const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];\n\n"
          "                if (ncols == 0) continue;\n\n"
          "                const auto cexp = c_exps[k];\n\n"
          "                const auto qexp = pexp + cexp;\n\n")
    c += ("                // NOTE: the total exponent is raised to the angular momentum by\n"
          "                // repeated multiplication, as the angular momentum is small.\n\n"
          "                auto faux = fcoul * anorm * b_norms[j] * c_norms[k] / (pexp * cexp * std::sqrt(qexp));\n\n"
          "                const auto frq = 1.0 / qexp;\n\n")
    c += "".join("                faux *= frq;\n\n" for _ in range(l))
    for l1 in range(l + 1):
        terms = ["faux"] + ["aexp"] * l1 + ["pexp"] * (l - l1)
        c += f"                const auto f_{l1} = {' * '.join(terms)};\n\n"
    c += f"                const auto *bvals = boys.data({l + 1}, k);\n\n"
    accs = ", ".join(f"acc_{l1}" for l1 in range(l + 1))
    c += (f"#pragma omp simd aligned({accs}, e_ab, bvals : simd::cache_line_size())\n"
          "                for (size_t l = 0; l < ncols; l++)\n                {\n"
          "                    const auto fval = e_ab[l] * bvals[l];\n\n")
    for l1 in range(l + 1):
        c += f"                    acc_{l1}[l] += f_{l1} * fval;\n\n"
    c = c.rstrip("\n") + "\n                }\n            }\n        }\n    }\n\n"

    c += ("    // NOTE: the bidegrees are accumulated into the angular components one at a\n"
          "    // time, so that no loop holds the harmonics of every degree at once and the\n"
          "    // vectorizer keeps its registers.\n\n")
    c += f"    auto components = CSimdMatrix({ncomps}, nmax);\n\n    components.zero();\n\n"
    for m in range(-l, l + 1):
        c += f"    auto *out_{name(m)} = components.data({m + l});\n"
    c += "\n"
    outs = ", ".join(f"out_{name(m)}" for m in range(-l, l + 1))

    for l1 in range(l + 1):
        l2 = l - l1
        if l1 == 0 or l2 == 0:
            side, arr = ("bc_harmonics", "atoms on c side") if l1 == 0 else ("ab_harmonics", "atom pairs")
            deg = l2 if l1 == 0 else l1
            c += (f"    // the bidegree of degree {WORD[l1]} on the atom pairs and {WORD[l2]} on the atoms\n"
                  f"    // on c side, whose coefficients are one: the harmonic of the other side is\n"
                  f"    // of degree zero and is one for every atom pair\n\n    {{\n")
            for m in range(-l, l + 1):
                c += f"        const auto *h_{name(m)} = {side}[{deg - 1}].data({m + l});\n"
            hs = ", ".join(f"h_{name(m)}" for m in range(-l, l + 1))
            c += (f"\n#pragma omp simd aligned({outs}, acc_{l1}, {hs} : simd::cache_line_size())\n"
                  "        for (size_t k = 0; k < nmax; k++)\n        {\n"
                  f"            const auto f = acc_{l1}[k];\n\n")
            for m in range(-l, l + 1):
                c += f"            out_{name(m)}[k] += f * h_{name(m)}[k];\n"
            c += "        }\n    }\n\n"
            continue

        rows, keys = {}, []
        for m1 in range(-l1, l1 + 1):
            for m2 in range(-l2, l2 + 1):
                if (l1, m1, m2) in coeffs:
                    keys.append((l1, m1, m2))
                    rows[(l1, m1, m2)] = f"r_{name(m1)}_{name(m2)}"
        c += (f"    // the bidegree of degree {WORD[l1]} on the atom pairs and {WORD[l2]} on the atoms\n"
              f"    // on c side, whose {len(keys)} products of harmonics are formed once and read by\n"
              f"    // the angular components which carry them\n\n    {{\n")
        c += f"        auto products = CSimdMatrix({len(keys)}, nmax);\n\n"
        for m1 in range(-l1, l1 + 1):
            c += f"        const auto *p_{name(m1)} = ab_harmonics[{l1 - 1}].data({m1 + l1});\n"
        c += "\n"
        for m2 in range(-l2, l2 + 1):
            c += f"        const auto *q_{name(m2)} = bc_harmonics[{l2 - 1}].data({m2 + l2});\n"
        c += "\n"
        for idx, key in enumerate(keys):
            c += f"        auto *{rows[key]} = products.data({idx});\n"
        ps = ", ".join(f"p_{name(m1)}" for m1 in range(-l1, l1 + 1))
        qs = ", ".join(f"q_{name(m2)}" for m2 in range(-l2, l2 + 1))
        rs = ", ".join(rows[k] for k in keys)
        c += (f"\n#pragma omp simd aligned({rs}, {ps}, {qs} : simd::cache_line_size())\n"
              "        for (size_t k = 0; k < nmax; k++)\n        {\n")
        for key in keys:
            _, m1, m2 = key
            c += f"            {rows[key]}[k] = p_{name(m1)}[k] * q_{name(m2)}[k];\n"
        c += "        }\n"
        c += ("\n        // NOTE: an angular component is accumulated by a loop of its own, as\n"
              "        // it reads only the products which carry it and the vectorizer would\n"
              "        // otherwise hold every product of the bidegree at once.\n")
        for m in range(-l, l + 1):
            ts = _terms(coeffs, m, rows)
            if not ts:
                continue
            used = [rows[k] for k in keys if coeffs.get(k, {}).get(m) is not None]
            c += (f"\n#pragma omp simd aligned(out_{name(m)}, acc_{l1}, {', '.join(used)} : simd::cache_line_size())\n"
                  "        for (size_t k = 0; k < nmax; k++)\n        {\n"
                  f"            out_{name(m)}[k] += acc_{l1}[k] * ({' + '.join(ts).replace('+ -', '- ')});\n"
                  "        }\n")
        c += "    }\n\n"

    c += ("    // NOTE: the atom pairs beyond the reach of every triple of primitives have no\n"
          "    // contribution and are set to zero.\n\n")
    c += (f"    for (size_t m = 0; m < {ncomps}; m++)\n    {{\n"
          "        const auto *row = components.data(m);\n\n"
          "        std::copy(row, row + nmax, slices[m]);\n\n"
          "        std::fill(slices[m] + nmax, slices[m] + npairs, 0.0);\n    }\n}\n\n"
          "}  // namespace simdt3ceri\n")
    return c, hpp


for l in [int(a) for a in sys.argv[1:]] or [3, 4, 5, 6]:
    cpp, hpp = emit(l)
    u = LETTER[l]
    io.open(SRC + f'SimdThreeCenterElectronRepulsionRecSS{u}.cpp', 'w').write(cpp)
    io.open(SRC + f'SimdThreeCenterElectronRepulsionRecSS{u}.hpp', 'w').write(hpp)
    print(f"wrote SimdThreeCenterElectronRepulsionRecSS{u}.cpp ({len(cpp.splitlines())} lines) and .hpp")
