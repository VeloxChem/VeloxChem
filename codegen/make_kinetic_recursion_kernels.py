"""Emits a kinetic energy kernel of two basis functions of non-zero angular
momentum from the recursion data."""
import io
import os
from fractions import Fraction as F

# NOTE: the root of the checkout is taken from the location of this file, so the
# generators run wherever the repository is placed.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src") + "/"

DST = os.path.join(SRC, "simd_t2c_kinetic_energy") + "/"
LIC = io.open(SRC + 'simd_t2c_overlap/SimdOverlapRecSS.cpp').read().split('#include "SimdOverlapRecSS.hpp"')[0]
LIC_H = io.open(SRC + 'simd_t2c_overlap/SimdOverlapRecSS.hpp').read().split('#ifndef')[0]

LET = {0:'s',1:'p',2:'d',3:'f',4:'g',5:'h',6:'i'}
WORD = {0:'zero',1:'one',2:'two',3:'three',4:'four',5:'five',6:'six'}
# NOTE: a loop is grouped so that the values it loads, the prefactors and the
# harmonics and the powers of the squared distance together, stay within the
# thirty two vector registers of the machine. Grouping by a fixed number of
# angular components instead lets the large kernels load sixty values in one
# loop, and the vectorizer then gives up and emits scalar code for all of them.
LOAD_BUDGET = 18


def parse(path):
    la = lb = None; sym = False
    E, G, rows = [], [], []
    cur = None
    for line in io.open(path):
        t = line.split()
        if not t or t[0].startswith('#'): continue
        if t[0] == 'SHELLS': la, lb = int(t[1]), int(t[2])
        elif t[0] == 'E': E.append(tuple(int(v) for v in t[2:6]))       # a b mu p
        elif t[0] == 'G': G.append(tuple(int(v) for v in t[2:5]))       # 2k L M
        elif t[0] == 'TABLE': sym = (len(t) > 4 and t[4] == 'SYMMETRIC')
        elif t[0] == 'ROW':
            cur = (int(t[1]), int(t[2]), [])
            rows.append(cur)
        elif t[0] == 'T':
            cur[2].append((int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5])))  # e num den rad g
    return la, lb, sym, E, G, rows


def mname(m):
    return "0" if m == 0 else (f"p{m}" if m > 0 else f"m{-m}")


def cname(num, den, rad):
    """The name of a coefficient, its sign carried by the expression."""
    n = abs(num)
    return f"f_{n}_{den}" if rad == 1 else f"fs_{n}_{den}_{rad}"


def clit(num, den, rad):
    n = abs(num)
    if rad == 1:
        if den == 1: return f"{float(n)}"
        return f"{n}.0 / {den}.0"
    if den == 1: return f"std::sqrt({n * n * rad}.0)"
    return f"std::sqrt({F(n * n * rad, den * den).numerator}.0 / {F(n * n * rad, den * den).denominator}.0)"


def emit(path):
    la, lb, sym, E, G, rows = parse(path)
    name = LET[la] + LET[lb]; up = name.upper()
    na, nb = 2 * la + 1, 2 * lb + 1
    nslots = na * nb
    lmax = max((g[1] for g in G), default=0)

    slot = lambda m, mp: (m + la) * nb + (mp + lb)
    stored = {slot(m, mp): i for i, (m, mp, _) in enumerate(rows)}

    # the harmonics and the powers of the squared distance which the terms use
    used_h = sorted({(g[1], g[2]) for g in G if g[1] > 0})
    coeffs = {}
    for _, _, terms in rows:
        for e, num, den, rad, g in terms:
            coeffs[cname(num, den, rad)] = clit(num, den, rad)

    ind = ' ' * (len(name) + 25)
    sig = (f"compute_{name}_kinetic_energy(double                         *values,\n"
           f"{ind}const size_t                    nvalues,\n"
           f"{ind}const CBasisFunction           &bra,\n"
           f"{ind}const CBasisFunction           &ket,\n"
           f"{ind}const std::vector<CSimdMatrix> &harmonics,\n"
           f"{ind}const CSimdMatrix              &coordinates,\n"
           f"{ind}const double                    threshold) -> void")

    s = LIC + f'#include "SimdKineticEnergyRec{up}.hpp"\n\n'
    s += "#include <algorithm>\n#include <ranges>\n#include <cmath>\n#include <string>\n\n"
    s += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "ScreeningFunc.hpp"\n'
          '#include "SimdAlign.hpp"\n#include "SimdDimensions.hpp"\n\n')
    s += "namespace simdkin {  // simdkin namespace\n\nauto\n" + sig + "\n{\n"
    s += (f'    if ((bra.get_angular_momentum() != {la}) || (ket.get_angular_momentum() != {lb}))\n    {{\n'
          '        errors::assertMsgCritical(\n'
          f'            false, std::string("SimdKineticEnergyFunc.compute_{name}_kinetic_energy: '
          f'Basis functions must be of angular momenta {WORD[la]} and {WORD[lb]}"));\n    }}\n\n')
    s += (f'    if (harmonics.size() < {lmax})\n    {{\n'
          '        errors::assertMsgCritical(\n'
          f'            false, std::string("SimdKineticEnergyFunc.compute_{name}_kinetic_energy: '
          f'Harmonics must reach angular momentum {WORD.get(lmax, str(lmax))}"));\n    }}\n\n')
    s += ('    if (nvalues > coordinates.number_of_columns())\n    {\n'
          '        errors::assertMsgCritical(\n'
          f'            false, std::string("SimdKineticEnergyFunc.compute_{name}_kinetic_energy: '
          'Number of values exceeds number of atom pairs"));\n    }\n\n')
    s += "    if (nvalues == 0) return;\n\n"
    s += ("    const auto &a_exps = bra.exponents();\n\n    const auto &b_exps = ket.exponents();\n\n"
          "    const auto &a_norms = bra.normalization_factors();\n\n"
          "    const auto &b_norms = ket.normalization_factors();\n\n"
          "    const auto nprim_a = a_exps.size();\n\n    const auto nprim_b = b_exps.size();\n\n"
          "    const auto nprims = nprim_a * nprim_b;\n\n")
    s += ("    // NOTE: the pairs of primitives are screened with the threshold of the\n"
          "    // integrals divided by their number, as their contributions accumulate into\n"
          "    // a single value and the error of the sum is bounded by the number of terms.\n\n"
          "    const auto dimensions = simdfunc::make_column_dimensions(\n"
          "        bra, ket, nvalues, coordinates, screenfunc::two_center_kinetic_energy_primitive_bound,\n"
          "        threshold / static_cast<double>(nprims));\n\n")
    s += ("    // NOTE: the buffer spans the atom pairs reached by the pair of primitives\n"
          "    // reaching furthest, which is searched for rather than assumed. The\n"
          "    // primitives are sorted by descending exponent, but the bound of a pair of\n"
          "    // primitives carries their prefactor as well as their decay, so a tighter\n"
          "    // pair with a larger prefactor reaches further than a more diffuse pair with\n"
          "    // a smaller one, and the last pair is not always the furthest reaching.\n\n"
          "    const auto nmax = *std::ranges::max_element(dimensions);\n\n")
    s += (f"    if (nmax == 0)\n    {{\n        std::fill(values, values + {nslots} * nvalues, 0.0);\n\n"
          "        return;\n    }\n\n")
    s += ("    // NOTE: the buffer holds the contracted prefactor of each exponent factor of\n"
          "    // the terms, as the integrals of the angular components are formed straight\n"
          "    // into the values and are not written a second time.\n\n"
          f"    auto buffer = CSimdMatrix({len(E)}, nmax);\n\n")
    for i in range(len(E)):
        s += f"    auto *pe_{i} = buffer.data({i});\n"
    s += "\n"
    for i in range(len(E)):
        s += f"    std::fill(pe_{i}, pe_{i} + nmax, 0.0);\n"
    s += "\n    const auto *ab_2 = coordinates.data(6);\n\n    constexpr auto fpi = mathconst::pi_value();\n\n"
    s += "    // accumulate the prefactor of each exponent factor over the pairs of primitives\n\n"
    s += ("    for (size_t i = 0; i < nprim_a; i++)\n    {\n        const auto aexp = a_exps[i];\n\n"
          "        const auto anorm = a_norms[i];\n\n"
          "        for (size_t j = 0; j < nprim_b; j++)\n        {\n"
          "            const auto ncols = dimensions[i * nprim_b + j];\n\n"
          "            if (ncols == 0) continue;\n\n"
          "            const auto bexp = b_exps[j];\n\n"
          "            const auto fexp = aexp + bexp;\n\n"
          "            const auto fmu = aexp * bexp / fexp;\n\n"
          "            const auto fovl = fpi / fexp;\n\n"
          "            const auto fbase = anorm * b_norms[j] * fovl * std::sqrt(fovl);\n\n")
    for i, (pa, pb, pm, pp) in enumerate(E):
        parts = ["fbase"]
        parts += ["aexp"] * pa + ["bexp"] * pb + ["fmu"] * pm
        expr = " * ".join(parts)
        if pp < 0:
            expr += " / (" + " * ".join(["fexp"] * (-pp)) + ")"
        s += f"            const auto ff_{i} = {expr};\n\n"
    al = ", ".join([f"pe_{i}" for i in range(len(E))] + ["ab_2"])
    s += (f"#pragma omp simd aligned({al} : simd::cache_line_size())\n"
          "            for (size_t k = 0; k < ncols; k++)\n            {\n"
          "                const auto fterm = std::exp(-fmu * ab_2[k]);\n\n")
    for i in range(len(E)):
        s += f"                pe_{i}[k] += ff_{i} * fterm;\n"
    s += "            }\n        }\n    }\n\n"

    # the harmonics rows and the value rows
    for (L, M) in used_h:
        s += f"    const auto *ph{L}_{mname(M)} = harmonics[{L - 1}].data({M + L});\n"
    s += "\n"
    for si, (idx, r) in enumerate(sorted(stored.items())):
        s += f"    auto *pc_{r} = values + {idx} * nvalues;\n"
    s += "\n"
    s += ("    // NOTE: the factors of the terms depend on the angular momenta alone, so they\n"
          "    // are formed once for the whole matrix instead of once for every atom pair.\n\n")
    for k in sorted(coeffs):
        s += f"    const auto {k} = {coeffs[k]};\n"
    s += "\n"

    maxpow = max((g[0] // 2 for g in G), default=0)
    def loads_of(rs):
        h = {(G[t[4]][1], G[t[4]][2]) for _, _, ts in rs for t in ts if G[t[4]][1] > 0}
        e = {t[0] for _, _, ts in rs for t in ts}
        pw = max((G[t[4]][0] // 2 for _, _, ts in rs for t in ts), default=0)
        return len(h) + len(e) + (1 if pw else 0)

    groups, cur = [], []
    for r in rows:
        if cur and loads_of(cur + [r]) > LOAD_BUDGET:
            groups.append(cur); cur = [r]
        else:
            cur.append(r)
    if cur:
        groups.append(cur)
    for gi, grp in enumerate(groups):
        es = sorted({t[0] for _, _, ts in grp for t in ts})
        hs = sorted({(G[t[4]][1], G[t[4]][2]) for _, _, ts in grp for t in ts if G[t[4]][1] > 0})
        pw = max((G[t[4]][0] // 2 for _, _, ts in grp for t in ts), default=0)
        outs = [rows.index(r) for r in grp]
        al = ([f"pe_{i}" for i in es] + [f"ph{L}_{mname(M)}" for L, M in hs] +
              (["ab_2"] if pw > 0 else []) + [f"pc_{o}" for o in outs])
        if gi:
            s += "\n"
        if len(groups) > 1:
            s += (f"    // NOTE: the angular components are formed in {len(groups)} loops, grouped so that\n"
                  "    // the values a loop loads stay within the vector registers of the machine.\n"
                  "    // Only the prefactors, the harmonics and the squared distance are loaded by\n"
                  "    // more than one loop, every other value once.\n\n")
        s += f"#pragma omp simd aligned({', '.join(al)} : simd::cache_line_size())\n"
        s += "    for (size_t k = 0; k < nmax; k++)\n    {\n"
        for i in es:
            s += f"        const auto e_{i} = pe_{i}[k];\n"
        if pw > 0:
            s += "\n        const auto r_2 = ab_2[k];\n"
            for j in range(2, pw + 1):
                s += f"        const auto r_{2 * j} = r_2 * r_{2 * (j - 1)};\n"
        if hs:
            s += "\n"
        for L, M in hs:
            s += f"        const auto h{L}_{mname(M)} = ph{L}_{mname(M)}[k];\n"
        s += "\n"
        for r in grp:
            o = rows.index(r)
            parts = []
            for e, num, den, rad, g in r[2]:
                k2, L, M = G[g]
                fac = [f"e_{e}"]
                if not (num == den == rad == 1):
                    fac.append(cname(num, den, rad))
                if k2 > 0:
                    fac.append("r_2" if k2 == 2 else f"r_{k2}")
                if L > 0:
                    fac.append(f"h{L}_{mname(M)}")
                parts.append(("-" if num < 0 else "+") + " " + " * ".join(fac))
            expr = " ".join(parts).lstrip("+ ").strip()
            s += f"        pc_{o}[k] = {expr};\n\n"
        s = s.rstrip("\n") + "\n    }\n"

    s += ("\n    // NOTE: the values of a combination of angular components are stored as one\n"
          "    // row of nvalues columns, with the component on bra side running slowest. The\n"
          "    // rows which the symmetry relates to an already formed one are copied from it,\n"
          "    // and the atom pairs beyond the reach of every pair of primitives are set to\n"
          "    // zero.\n\n")
    src = []
    for m in range(-la, la + 1):
        for mp in range(-lb, lb + 1):
            i = slot(m, mp)
            src.append(i if i in stored else slot(mp, m))
    s += f"    const size_t sources[{nslots}] = {{{', '.join(str(v) for v in src)}}};\n\n"
    s += (f"    for (size_t n = 0; n < {nslots}; n++)\n    {{\n"
          "        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);\n\n"
          "        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);\n    }\n")
    s += "}\n\n}  // namespace simdkin\n"

    h = LIC_H + f"#ifndef SimdKineticEnergyRec{up}_hpp\n#define SimdKineticEnergyRec{up}_hpp\n\n"
    h += "#include <cstddef>\n#include <vector>\n\n"
    h += '#include "BasisFunction.hpp"\n#include "SimdMatrix.hpp"\n\nnamespace simdkin {  // simdkin namespace\n\n'
    h += (f"/// @brief Computes the kinetic energy integrals of the basis functions of angular\n"
          f"/// momenta {WORD[la]} and {WORD[lb]} over the atom pairs of a block.\n"
          "/// @param values The values of the combination of basis functions to compute.\n"
          "/// @param nvalues The number of atom pairs the combination reaches.\n"
          "/// @param bra The basis function on bra side.\n"
          "/// @param ket The basis function on ket side.\n"
          f"/// @param harmonics The solid harmonics of the vectors between the atoms, to angular momentum {lmax}.\n"
          "/// @param coordinates The coordinates of the atom pairs.\n"
          "/// @param threshold The screening threshold of the integrals.\n")
    h += f"auto {sig};\n\n}}  // namespace simdkin\n\n#endif /* SimdKineticEnergyRec{up}_hpp */\n"
    return up, s, h


if __name__ == "__main__":
    import sys

    # the directory holding the kinetic_harmonics_*.txt recursion data
    DATA = os.environ.get("KINETIC_DATA", os.path.expanduser("~/Downloads"))
    for tag in sys.argv[1:]:
        up, cpp, hpp = emit(os.path.join(DATA, f"kinetic_harmonics_{tag}.txt"))
        io.open(DST + f'SimdKineticEnergyRec{up}.cpp', 'w').write(cpp)
        io.open(DST + f'SimdKineticEnergyRec{up}.hpp', 'w').write(hpp)
        print(f"wrote SimdKineticEnergyRec{up}.cpp, {len(cpp.splitlines())} lines")
