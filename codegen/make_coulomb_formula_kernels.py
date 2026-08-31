"""Emits the two-center electron repulsion kernels which follow a closed
formula: (s|J|s), and (s|J|l) and (l|J|s) for l from one to six."""
import io
import os

# NOTE: the root of the checkout is taken from the location of this file, so the
# generators run wherever the repository is placed.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src") + "/"

DST = os.path.join(SRC, "simd_t2c_electron_repulsion") + "/"
LIC = io.open(SRC + 'simd_t2c_overlap/SimdOverlapRecSS.cpp').read().split('#include "SimdOverlapRecSS.hpp"')[0]
LIC_H = io.open(SRC + 'simd_t2c_overlap/SimdOverlapRecSS.hpp').read().split('#ifndef')[0]

LET = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i'}
WORD = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six'}


def comps(l):
    """The names of the angular components, from m = -l to m = +l."""
    return [("m%d" % -m if m < 0 else ("0" if m == 0 else "p%d" % m)) for m in range(-l, l + 1)]


def emit(la, lb):
    """One kernel. Exactly one of the two angular momenta is zero."""
    l = la + lb
    name = LET[la] + LET[lb]
    up = name.upper()
    ncomp = 2 * l + 1
    on_bra = la > 0
    pad = ' ' * (len(name) + 28)

    sig = (f"compute_{name}_electron_repulsion(double               *values,\n"
           f"{pad}const size_t          nvalues,\n"
           f"{pad}const CBasisFunction &bra,\n"
           f"{pad}const CBasisFunction &ket,\n")
    if l > 0:
        sig += f"{pad}const CSimdMatrix    &harmonics,\n"
    sig += f"{pad}const CSimdMatrix    &coordinates) -> void"

    s = LIC + f'#include "SimdTwoCenterElectronRepulsionRec{up}.hpp"\n\n'
    s += "#include <algorithm>\n#include <cmath>\n#include <string>\n#include <vector>\n\n"
    s += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "SimdAlign.hpp"\n'
          '#include "SimdBoysFunc.hpp"\n#include "SimdVariableMatrix.hpp"\n\n')
    s += "namespace simdt2ceri {  // simdt2ceri namespace\n\nauto\n" + sig + "\n{\n"

    err = f"SimdTwoCenterElectronRepulsionRec{up}.compute_{name}_electron_repulsion"

    if l == 0:
        s += ('    if ((bra.get_angular_momentum() != 0) || (ket.get_angular_momentum() != 0))\n    {\n'
              '        errors::assertMsgCritical(\n'
              f'            false, std::string("{err}: '
              'Basis functions must be of zero angular momentum"));\n    }\n\n')
    else:
        s += (f'    if ((bra.get_angular_momentum() != {la}) || (ket.get_angular_momentum() != {lb}))\n    {{\n'
              '        errors::assertMsgCritical(\n'
              f'            false, std::string("{err}: '
              f'Basis functions must be of angular momenta {WORD[la]} and {WORD[lb]}"));\n    }}\n\n')
        s += (f'    if (harmonics.number_of_rows() != {ncomp})\n    {{\n'
              '        errors::assertMsgCritical(\n'
              f'            false, std::string("{err}: '
              f'Harmonics must have {ncomp} rows"));\n    }}\n\n')

    s += ('    if (nvalues > coordinates.number_of_columns())\n    {\n'
          '        errors::assertMsgCritical(\n'
          f'            false, std::string("{err}: '
          'Number of values exceeds number of atom pairs"));\n    }\n\n')
    s += "    if (nvalues == 0) return;\n\n"
    s += ("    const auto &a_exps = bra.exponents();\n\n    const auto &b_exps = ket.exponents();\n\n"
          "    const auto &a_norms = bra.normalization_factors();\n\n"
          "    const auto &b_norms = ket.normalization_factors();\n\n"
          "    const auto nprim_a = a_exps.size();\n\n    const auto nprim_b = b_exps.size();\n\n"
          "    const auto nprims = nprim_a * nprim_b;\n\n")
    s += ("    // NOTE: the squared distances of the atom pairs are carried by the\n"
          "    // coordinates, so that they are formed once for the whole block instead of\n"
          "    // once for every combination of basis functions.\n\n"
          "    const auto *ab_2 = coordinates.data(6);\n\n")

    # the arguments of the Boys function, one row per pair of primitives

    s += ("    // NOTE: the pairs of primitives are not screened and every one of them\n"
          "    // reaches every atom pair, as the Coulomb operator decays as the inverse of\n"
          "    // the interatomic distance, so every row spans the atom pairs of the block.\n\n"
          f"    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), {l + 2});\n\n")
    s += ("    // set up the arguments of the Boys function of each pair of primitives\n\n"
          "    for (size_t i = 0; i < nprim_a; i++)\n    {\n        const auto aexp = a_exps[i];\n\n"
          "        for (size_t j = 0; j < nprim_b; j++)\n        {\n"
          "            const auto fmu = aexp * b_exps[j] / (aexp + b_exps[j]);\n\n"
          "            auto *bargs = boys.data(0, i * nprim_b + j);\n\n"
          "#pragma omp simd aligned(bargs, ab_2 : simd::cache_line_size())\n"
          "            for (size_t k = 0; k < nvalues; k++)\n            {\n"
          "                bargs[k] = fmu * ab_2[k];\n            }\n        }\n    }\n\n")
    s += ("    // NOTE: the Boys function of every pair of primitives is computed by one\n"
          f"    // call, which fills the orders zero to {WORD[l]} of every row. The integrals need\n"
          f"    // the order {WORD[l]} alone, and the lower orders are formed on the way to it by\n"
          "    // the recursion the Boys function is evaluated with.\n\n"
          "    simdfunc::compute_boys_function(boys);\n\n")

    if l == 0:
        s += ("    // NOTE: the integrals of all pairs of primitives are accumulated in a single\n"
              "    // row, which starts at a cache line boundary.\n\n"
              "    auto buffer = CSimdMatrix(1, nvalues);\n\n    buffer.zero();\n\n    auto *prim = buffer.data(0);\n\n")
    else:
        s += ("    // NOTE: the first row accumulates the factor which the angular components\n"
              "    // share, and the remaining rows hold the integrals of the angular components.\n\n"
              f"    auto buffer = CSimdMatrix({ncomp + 1}, nvalues);\n\n"
              "    auto *prim = buffer.data(0);\n\n    std::fill(prim, prim + nvalues, 0.0);\n\n")

    s += ("    constexpr auto fpi = mathconst::pi_value();\n\n"
          "    // NOTE: the two-center repulsion of a pair of s type primitives is two pi to\n"
          "    // the five halves over the exponents times the square root of their sum,\n"
          "    // times the Boys function of the order the angular momenta ask for.\n\n"
          "    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);\n\n")
    s += "    // accumulate the integrals of each pair of primitives\n\n"
    s += ("    for (size_t i = 0; i < nprim_a; i++)\n    {\n        const auto aexp = a_exps[i];\n\n"
          "        const auto anorm = a_norms[i];\n\n"
          "        for (size_t j = 0; j < nprim_b; j++)\n        {\n"
          "            const auto fexp = aexp + b_exps[j];\n\n"
          "            auto ffact = fcoul * anorm * b_norms[j] / (aexp * b_exps[j] * std::sqrt(fexp));\n\n")
    if l > 0:
        # NOTE: the ratio belongs to the side carrying no angular momentum: the
        # bra exponent where the momentum is on ket, the negated ket exponent
        # where it is on bra, as the overlap kernels of the same combinations do.
        ratio = "-b_exps[j] / fexp" if on_bra else "aexp / fexp"
        s += ("            // NOTE: the ratio is raised to the angular momentum by repeated\n"
              "            // multiplication, as the angular momentum is small.\n\n"
              f"            const auto frat = {ratio};\n\n")
        for _ in range(l):
            s += "            ffact *= frat;\n\n"
    s += (f"            const auto *bvals = boys.data({l + 1}, i * nprim_b + j);\n\n"
          "            // NOTE: the row of the buffer and the row of the Boys function start at\n"
          "            // a cache line boundary, so the loop is vectorized with aligned loads\n"
          "            // and stores.\n\n"
          "#pragma omp simd aligned(prim, bvals : simd::cache_line_size())\n"
          "            for (size_t k = 0; k < nvalues; k++)\n            {\n"
          "                prim[k] += ffact * bvals[k];\n"
          "            }\n        }\n    }\n\n")

    if l == 0:
        s += "    std::copy(prim, prim + nvalues, values);\n"
    else:
        cs = comps(l)
        s += ("    // NOTE: the integral of an angular component is the accumulated factor\n"
              "    // times the solid harmonic of that component.\n\n")
        for i, c in enumerate(cs):
            s += f"    const auto *ph_{c} = harmonics.data({i});\n"
        s += "\n"
        for i, c in enumerate(cs):
            s += f"    auto *pb_{c} = buffer.data({i + 1});\n"
        s += "\n"
        al = ["prim"] + [f"ph_{c}" for c in cs] + [f"pb_{c}" for c in cs]
        s += f"#pragma omp simd aligned({', '.join(al)} : simd::cache_line_size())\n"
        s += "    for (size_t k = 0; k < nvalues; k++)\n    {\n        const auto fval = prim[k];\n\n"
        for c in cs:
            s += f"        pb_{c}[k] = fval * ph_{c}[k];\n"
        s += "    }\n\n"
        s += ("    // NOTE: the values of an angular component are stored as one row of nvalues\n"
              "    // columns.\n\n"
              f"    for (size_t m = 0; m < {ncomp}; m++)\n    {{\n"
              "        const auto *pb = buffer.data(m + 1);\n\n"
              "        std::copy(pb, pb + nvalues, values + m * nvalues);\n    }\n")
    s += "}\n\n}  // namespace simdt2ceri\n"

    h = LIC_H + f"#ifndef SimdTwoCenterElectronRepulsionRec{up}_hpp\n#define SimdTwoCenterElectronRepulsionRec{up}_hpp\n\n"
    h += "#include <cstddef>\n\n"
    h += '#include "BasisFunction.hpp"\n#include "SimdMatrix.hpp"\n\nnamespace simdt2ceri {  // simdt2ceri namespace\n\n'
    h += ("/// @brief Computes the two-center electron repulsion integrals of the basis\n"
          f"/// functions of angular momenta {WORD[la]} and {WORD[lb]} over the atom pairs of a block.\n"
          "/// @param values The values of the combination of basis functions to compute.\n"
          "/// @param nvalues The number of atom pairs of the block.\n"
          "/// @param bra The basis function on bra side.\n"
          "/// @param ket The basis function on ket side.\n")
    if l > 0:
        h += f"/// @param harmonics The solid harmonics of angular momentum {WORD[l]} of the vectors between the atoms.\n"
    h += "/// @param coordinates The coordinates of the atom pairs.\n"
    h += ("/// @note There is no screening threshold, as the integrals are not screened: the\n"
          "/// Coulomb operator decays as the inverse of the interatomic distance and neither\n"
          "/// an atom pair nor a pair of primitives falls below a threshold.\n")
    h += f"auto {sig};\n\n}}  // namespace simdt2ceri\n\n#endif /* SimdTwoCenterElectronRepulsionRec{up}_hpp */\n"
    return s, h


made = []
for la, lb in [(0, 0)] + [(0, l) for l in range(1, 7)] + [(l, 0) for l in range(1, 7)]:
    cpp, hpp = emit(la, lb)
    up = (LET[la] + LET[lb]).upper()
    io.open(DST + f'SimdTwoCenterElectronRepulsionRec{up}.cpp', 'w').write(cpp)
    io.open(DST + f'SimdTwoCenterElectronRepulsionRec{up}.hpp', 'w').write(hpp)
    made.append(up)
print("wrote:", " ".join(made))
