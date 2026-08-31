"""Emits an overlap, kinetic energy or two-center Coulomb kernel of two basis
functions of non-zero angular momentum from the recursion data.

The three kernels share the format of the data and the shape of the generated
code. The overlap and the kinetic energy differ only in the name they are given,
the namespace they sit in and the bound they screen the pairs of primitives with.
The Coulomb kernel is not screened at all and its source integral carries an
order, so its data has an order on every term and its code replaces the
exponential of a pair of primitives with the Boys function of that order."""
import io
import os
from fractions import Fraction as F

# NOTE: the root of the checkout is taken from the location of this file, so the
# generators run wherever the repository is placed.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src") + "/"

# NOTE: the destination and the file the licence is taken from follow the kernel.
LIC = io.open(SRC + 'simd_t2c_overlap/SimdOverlapRecSS.cpp').read().split('#include "SimdOverlapRecSS.hpp"')[0]
LIC_H = io.open(SRC + 'simd_t2c_overlap/SimdOverlapRecSS.hpp').read().split('#ifndef')[0]

LET = {0:'s',1:'p',2:'d',3:'f',4:'g',5:'h',6:'i'}

# NOTE: what distinguishes the two kernels. The data of both is read the same way
# and the code of both has the same shape.
KERNELS = {
    "S": dict(word="overlap", ns="simdovl", dir="simd_t2c_overlap", stem="SimdOverlapRec",
              func="{name}_overlap", header="SimdOverlapFunc",
              bound="screenfunc::two_center_overlap_primitive_bound"),
    "T": dict(word="kinetic energy", ns="simdkin", dir="simd_t2c_kinetic_energy",
              stem="SimdKineticEnergyRec", func="{name}_kinetic_energy",
              header="SimdKineticEnergyFunc",
              bound="screenfunc::two_center_kinetic_energy_primitive_bound"),
    "J": dict(word="two-center electron repulsion", ns="simdt2ceri",
              dir="simd_t2c_electron_repulsion", stem="SimdTwoCenterElectronRepulsionRec",
              func="{name}_electron_repulsion", header="SimdTwoCenterElectronRepulsionFunc",
              bound=None),
}
WORD = {0:'zero',1:'one',2:'two',3:'three',4:'four',5:'five',6:'six'}
# NOTE: a loop is grouped so that the values it loads, the prefactors and the
# harmonics and the powers of the squared distance together, stay within the
# thirty two vector registers of the machine. Grouping by a fixed number of
# angular components instead lets the large kernels load sixty values in one
# loop, and the vectorizer then gives up and emits scalar code for all of them.
LOAD_BUDGET = 18
# NOTE: the Coulomb kernels count the rows a loop writes as well as the values it
# loads, and their budget is the two together. Their loops hold more angular
# components than the kinetic energy loops of the same angular momenta, as one
# exponent factor of theirs serves one order of the Boys function where the
# kinetic energy has one per power of the reduced exponent, and counting the
# loads alone lets (d|J|d) and (p|J|f) fall back to scalar code. Twenty is the
# largest budget at which no kernel of the set falls back, and the vector
# operations of the largest kernels are within a percent of what any smaller
# budget reaches.
LOAD_BUDGET_J = 20


def parse(path):
    la = lb = None; sym = False; kernel = 'T'; orders = None
    E, G, rows = [], [], []
    cur = None
    for line in io.open(path):
        t = line.split()
        if not t or t[0].startswith('#'): continue
        if t[0] == 'KERNEL': kernel = t[1]
        elif t[0] == 'SHELLS': la, lb = int(t[1]), int(t[2])
        elif t[0] == 'ORDERS': orders = (int(t[1]), int(t[2]))
        elif t[0] == 'E': E.append(tuple(int(v) for v in t[2:6]))       # a b mu p
        elif t[0] == 'G': G.append(tuple(int(v) for v in t[2:5]))       # 2k L M
        elif t[0] == 'TABLE': sym = (len(t) > 4 and t[4] == 'SYMMETRIC')
        elif t[0] == 'ROW':
            cur = (int(t[1]), int(t[2]), [])
            rows.append(cur)
        elif t[0] == 'T':
            # NOTE: the Coulomb data carries the order of the Boys function of a
            # term as a sixth field. The overlap and the kinetic energy have no
            # order, and zero stands in for it there and is never read.
            n = int(t[6]) if len(t) > 6 else 0
            cur[2].append((int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), n))  # e num den rad g n
    return la, lb, sym, E, G, rows, kernel, orders


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
    la, lb, sym, E, G, rows, kernel, orders = parse(path)
    K = KERNELS[kernel]
    # NOTE: the Coulomb source integral carries an order, and every exponent
    # factor of the data is used with one order alone, so an accumulator per
    # exponent factor still suffices and the order follows from the factor.
    boys = kernel == 'J'
    border = {}
    for _, _, terms in rows:
        for e, num, den, rad, g, n in terms:
            border.setdefault(e, n)
    maxorder = max(border.values()) if boys else 0
    name = LET[la] + LET[lb]; up = name.upper()
    DST = os.path.join(SRC, K["dir"]) + "/"
    fn = K["func"].format(name=name)
    na, nb = 2 * la + 1, 2 * lb + 1
    nslots = na * nb
    lmax = max((g[1] for g in G), default=0)

    slot = lambda m, mp: (m + la) * nb + (mp + lb)
    stored = {slot(m, mp): i for i, (m, mp, _) in enumerate(rows)}

    # the harmonics and the powers of the squared distance which the terms use
    used_h = sorted({(g[1], g[2]) for g in G if g[1] > 0})
    coeffs = {}
    for _, _, terms in rows:
        for e, num, den, rad, g, n in terms:
            # NOTE: a coefficient of one is not named, as the sign is carried by
            # the term and multiplying by one is not written out.
            if abs(num) == den == rad == 1: continue
            coeffs[cname(num, den, rad)] = clit(num, den, rad)

    # NOTE: the indent is the one the kinetic energy kernels in the repository
    # carry, so that regenerating them reproduces them exactly. The overlap
    # kernels came from another generator and are indented one space less.
    ind = ' ' * (len(fn) + (9 if boys else 10))
    sig = (f"compute_{fn}(double                         *values,\n"
           f"{ind}const size_t                    nvalues,\n"
           f"{ind}const CBasisFunction           &bra,\n"
           f"{ind}const CBasisFunction           &ket,\n"
           f"{ind}const std::vector<CSimdMatrix> &harmonics,\n"
           f"{ind}const CSimdMatrix              &coordinates")
    sig += ") -> void" if boys else (f",\n{ind}const double                    threshold) -> void")

    s = LIC + f'#include "{K["stem"]}{up}.hpp"\n\n'
    if boys:
        s += "#include <algorithm>\n#include <cmath>\n#include <string>\n#include <vector>\n\n"
        s += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "SimdAlign.hpp"\n'
              '#include "SimdBoysFunc.hpp"\n#include "SimdVariableMatrix.hpp"\n\n')
    else:
        s += "#include <algorithm>\n#include <ranges>\n#include <cmath>\n#include <string>\n\n"
        s += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "ScreeningFunc.hpp"\n'
              '#include "SimdAlign.hpp"\n#include "SimdDimensions.hpp"\n\n')
    s += f'namespace {K["ns"]} {{  // {K["ns"]} namespace\n\nauto\n' + sig + "\n{\n"
    s += (f'    if ((bra.get_angular_momentum() != {la}) || (ket.get_angular_momentum() != {lb}))\n    {{\n'
          '        errors::assertMsgCritical(\n'
          f'            false, std::string("{K["header"]}.compute_{fn}: '
          f'Basis functions must be of angular momenta {WORD[la]} and {WORD[lb]}"));\n    }}\n\n')
    s += (f'    if (harmonics.size() < {lmax})\n    {{\n'
          '        errors::assertMsgCritical(\n'
          f'            false, std::string("{K["header"]}.compute_{fn}: '
          f'Harmonics must reach angular momentum {WORD.get(lmax, str(lmax))}"));\n    }}\n\n')
    s += ('    if (nvalues > coordinates.number_of_columns())\n    {\n'
          '        errors::assertMsgCritical(\n'
          f'            false, std::string("{K["header"]}.compute_{fn}: '
          'Number of values exceeds number of atom pairs"));\n    }\n\n')
    s += "    if (nvalues == 0) return;\n\n"
    s += ("    const auto &a_exps = bra.exponents();\n\n    const auto &b_exps = ket.exponents();\n\n"
          "    const auto &a_norms = bra.normalization_factors();\n\n"
          "    const auto &b_norms = ket.normalization_factors();\n\n"
          "    const auto nprim_a = a_exps.size();\n\n    const auto nprim_b = b_exps.size();\n\n"
          "    const auto nprims = nprim_a * nprim_b;\n\n")
    if boys:
        s += ("    const auto *ab_2 = coordinates.data(6);\n\n")
        s += ("    // NOTE: the pairs of primitives are not screened and every one of them\n"
              "    // reaches every atom pair, as the Coulomb operator decays as the inverse of\n"
              "    // the interatomic distance, so every row spans the atom pairs of the block.\n\n"
              f"    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), {maxorder + 2});\n\n")
        s += ("    // set up the arguments of the Boys function of each pair of primitives\n\n"
              "    for (size_t i = 0; i < nprim_a; i++)\n    {\n        const auto aexp = a_exps[i];\n\n"
              "        for (size_t j = 0; j < nprim_b; j++)\n        {\n"
              "            const auto fmu = aexp * b_exps[j] / (aexp + b_exps[j]);\n\n"
              "            auto *bargs = boys.data(0, i * nprim_b + j);\n\n"
              "#pragma omp simd aligned(bargs, ab_2 : simd::cache_line_size())\n"
              "            for (size_t k = 0; k < nvalues; k++)\n            {\n"
              "                bargs[k] = fmu * ab_2[k];\n            }\n        }\n    }\n\n")
        s += ("    // NOTE: the Boys function of every pair of primitives is computed by one\n"
              f"    // call, which fills the orders 0 to {maxorder} of every row. The terms read the\n"
              f"    // orders {min(border.values())} to {maxorder} alone, and the orders below them are formed on the\n"
              "    // way to them by the recursion the Boys function is evaluated with.\n\n"
              "    simdfunc::compute_boys_function(boys);\n\n")
        s += ("    // NOTE: the buffer holds the contracted prefactor of each exponent factor of\n"
              "    // the terms, as the integrals of the angular components are formed straight\n"
              "    // into the values and are not written a second time. Every exponent factor is\n"
              "    // used with one order of the Boys function alone, so the order follows from\n"
              "    // the factor and one accumulator per factor suffices.\n\n"
              f"    auto buffer = CSimdMatrix({len(E)}, nvalues);\n\n")
        for i in range(len(E)):
            s += f"    auto *pe_{i} = buffer.data({i});\n"
        s += "\n"
        for i in range(len(E)):
            s += f"    std::fill(pe_{i}, pe_{i} + nvalues, 0.0);\n"
        s += "\n    constexpr auto fpi = mathconst::pi_value();\n\n"
        s += ("    // NOTE: the two-center repulsion of a pair of s type primitives is two pi to\n"
              "    // the five halves over the exponents times the square root of their sum,\n"
              "    // times the Boys function of the order the terms ask for.\n\n"
              "    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);\n\n")
        s += "    // accumulate the prefactor of each exponent factor over the pairs of primitives\n\n"
        s += ("    for (size_t i = 0; i < nprim_a; i++)\n    {\n        const auto aexp = a_exps[i];\n\n"
              "        const auto anorm = a_norms[i];\n\n"
              "        for (size_t j = 0; j < nprim_b; j++)\n        {\n"
              "            const auto bexp = b_exps[j];\n\n"
              "            const auto fexp = aexp + bexp;\n\n"
              "            const auto fmu = aexp * bexp / fexp;\n\n"
              "            const auto fbase = fcoul * anorm * b_norms[j] / (aexp * bexp * std::sqrt(fexp));\n\n")
        for i, (pa, pb, pm, pp) in enumerate(E):
            parts = ["fbase"] + ["aexp"] * pa + ["bexp"] * pb + ["fmu"] * pm
            expr = " * ".join(parts)
            if pp < 0:
                expr += " / (" + " * ".join(["fexp"] * (-pp)) + ")"
            s += f"            const auto ff_{i} = {expr};\n\n"
        for i in range(len(E)):
            s += f"            const auto *bv_{i} = boys.data({border[i] + 1}, i * nprim_b + j);\n\n"
        al = ", ".join([f"pe_{i}" for i in range(len(E))] + [f"bv_{i}" for i in range(len(E))])
        s += (f"#pragma omp simd aligned({al} : simd::cache_line_size())\n"
              "            for (size_t k = 0; k < nvalues; k++)\n            {\n")
        for i in range(len(E)):
            s += f"                pe_{i}[k] += ff_{i} * bv_{i}[k];\n"
        s += "            }\n        }\n    }\n\n"
    else:
        s += ("    // NOTE: the pairs of primitives are screened with the threshold of the\n"
              "    // integrals divided by their number, as their contributions accumulate into\n"
              "    // a single value and the error of the sum is bounded by the number of terms.\n\n"
              "    const auto dimensions = simdfunc::make_column_dimensions(\n"
              f'        bra, ket, nvalues, coordinates, {K["bound"]},\n'
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
        # NOTE: the rows a loop writes are counted for the Coulomb kernels and not
        # for the other two, whose grouping is the one their kernels in the
        # repository were emitted with.
        return len(h) + len(e) + (1 if pw else 0) + (len(rs) if boys else 0)

    groups, cur = [], []
    for r in rows:
        if cur and loads_of(cur + [r]) > (LOAD_BUDGET_J if boys else LOAD_BUDGET):
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
        s += f"    for (size_t k = 0; k < {'nvalues' if boys else 'nmax'}; k++)\n    {{\n"
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
            for e, num, den, rad, g, n in r[2]:
                k2, L, M = G[g]
                fac = [f"e_{e}"]
                if not (abs(num) == den == rad == 1):
                    fac.append(cname(num, den, rad))
                if k2 > 0:
                    fac.append("r_2" if k2 == 2 else f"r_{k2}")
                if L > 0:
                    fac.append(f"h{L}_{mname(M)}")
                parts.append(("-" if num < 0 else "+") + " " + " * ".join(fac))
            expr = " ".join(parts).lstrip("+ ").strip()
            s += f"        pc_{o}[k] = {expr};\n\n"
        s = s.rstrip("\n") + "\n    }\n"

    if boys:
        s += ("\n    // NOTE: the values of a combination of angular components are stored as one\n"
              "    // row of nvalues columns, with the component on bra side running slowest. The\n"
              "    // rows which the symmetry relates to an already formed one are copied from it.\n\n")
    else:
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
    if boys:
        s += (f"    for (size_t n = 0; n < {nslots}; n++)\n    {{\n"
              "        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);\n    }\n")
    else:
        s += (f"    for (size_t n = 0; n < {nslots}; n++)\n    {{\n"
              "        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);\n\n"
              "        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);\n    }\n")
    s += f'}}\n\n}}  // namespace {K["ns"]}\n'

    h = LIC_H + f'#ifndef {K["stem"]}{up}_hpp\n#define {K["stem"]}{up}_hpp\n\n'
    h += "#include <cstddef>\n#include <vector>\n\n"
    h += f'#include "BasisFunction.hpp"\n#include "SimdMatrix.hpp"\n\nnamespace {K["ns"]} {{  // {K["ns"]} namespace\n\n'
    h += (f'/// @brief Computes the {K["word"]} integrals of the basis functions of angular\n'
          f"/// momenta {WORD[la]} and {WORD[lb]} over the atom pairs of a block.\n"
          "/// @param values The values of the combination of basis functions to compute.\n"
          f"/// @param nvalues The number of atom pairs {'of the block' if boys else 'the combination reaches'}.\n"
          "/// @param bra The basis function on bra side.\n"
          "/// @param ket The basis function on ket side.\n"
          f"/// @param harmonics The solid harmonics of the vectors between the atoms, to angular momentum {lmax}.\n"
          "/// @param coordinates The coordinates of the atom pairs.\n")
    h += ("/// @note There is no screening threshold, as the integrals are not screened: the\n"
          "/// Coulomb operator decays as the inverse of the interatomic distance and neither\n"
          "/// an atom pair nor a pair of primitives falls below a threshold.\n"
          if boys else "/// @param threshold The screening threshold of the integrals.\n")
    h += f'auto {sig};\n\n}}  // namespace {K["ns"]}\n\n#endif /* {K["stem"]}{up}_hpp */\n'
    return up, s, h, K


if __name__ == "__main__":
    import sys

    # the directory holding the kinetic_harmonics_*.txt recursion data
    # the recursion data ships with the repository
    DATA = os.environ.get("RECURSION_DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))
    KIND = os.environ.get("RECURSION_KIND", "kinetic")
    for tag in sys.argv[1:]:
        path = os.path.join(DATA, f"{KIND}_harmonics_{tag}.txt")
        up, cpp, hpp, K = emit(path)
        dst = os.path.join(SRC, K["dir"]) + "/"
        io.open(dst + f'{K["stem"]}{up}.cpp', 'w').write(cpp)
        io.open(dst + f'{K["stem"]}{up}.hpp', 'w').write(hpp)
        print(f'wrote {K["stem"]}{up}.cpp, {len(cpp.splitlines())} lines')
