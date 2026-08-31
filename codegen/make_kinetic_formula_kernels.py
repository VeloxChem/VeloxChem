"""Emits the kinetic energy kernels which follow a closed formula: (s|T|s), and
(s|T|l) and (l|T|s) for l from one to six."""
import io
import os

# NOTE: the root of the checkout is taken from the location of this file, so the
# generators run wherever the repository is placed.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src") + "/"

DST = os.path.join(SRC, "simd_t2c_kinetic_energy") + "/"
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

    sig = (f"compute_{name}_kinetic_energy(double               *values,\n"
           f"{' ' * (len(name) + 25)}const size_t          nvalues,\n"
           f"{' ' * (len(name) + 25)}const CBasisFunction &bra,\n"
           f"{' ' * (len(name) + 25)}const CBasisFunction &ket,\n")
    if l > 0:
        sig += f"{' ' * (len(name) + 25)}const CSimdMatrix    &harmonics,\n"
    sig += (f"{' ' * (len(name) + 25)}const CSimdMatrix    &coordinates,\n"
            f"{' ' * (len(name) + 25)}const double          threshold) -> void")

    s = LIC + f'#include "SimdKineticEnergyRec{up}.hpp"\n\n'
    s += "#include <algorithm>\n#include <ranges>\n#include <cmath>\n#include <string>\n\n"
    s += ('#include "ErrorHandler.hpp"\n#include "MathConst.hpp"\n#include "ScreeningFunc.hpp"\n'
          '#include "SimdAlign.hpp"\n#include "SimdDimensions.hpp"\n\n')
    s += "namespace simdkin {  // simdkin namespace\n\nauto\n" + sig + "\n{\n"

    if l == 0:
        s += ('    if ((bra.get_angular_momentum() != 0) || (ket.get_angular_momentum() != 0))\n    {\n'
              '        errors::assertMsgCritical(\n'
              f'            false, std::string("SimdKineticEnergyFunc.compute_{name}_kinetic_energy: '
              'Basis functions must be of zero angular momentum"));\n    }\n\n')
    else:
        first, second = (la, lb) if on_bra else (la, lb)
        s += (f'    if ((bra.get_angular_momentum() != {la}) || (ket.get_angular_momentum() != {lb}))\n    {{\n'
              '        errors::assertMsgCritical(\n'
              f'            false, std::string("SimdKineticEnergyFunc.compute_{name}_kinetic_energy: '
              f'Basis functions must be of angular momenta {WORD[la]} and {WORD[lb]}"));\n    }}\n\n')
        s += (f'    if (harmonics.number_of_rows() != {ncomp})\n    {{\n'
              '        errors::assertMsgCritical(\n'
              f'            false, std::string("SimdKineticEnergyFunc.compute_{name}_kinetic_energy: '
              f'Harmonics must have {ncomp} rows"));\n    }}\n\n')

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
    s += (f"    if (nmax == 0)\n    {{\n        std::fill(values, values + {ncomp} * nvalues, 0.0);\n\n"
          "        return;\n    }\n\n")

    if l == 0:
        s += ("    // NOTE: the integrals of all pairs of primitives are accumulated in a single\n"
              "    // row, which starts at a cache line boundary and spans only the atom pairs\n"
              "    // reached by the furthest reaching pair of primitives.\n\n"
              "    auto buffer = CSimdMatrix(1, nmax);\n\n    buffer.zero();\n\n    auto *prim = buffer.data(0);\n\n")
    else:
        s += ("    // NOTE: the first row accumulates the factor which the angular components\n"
              "    // share, and the remaining rows hold the integrals of the angular components.\n\n"
              f"    auto buffer = CSimdMatrix({ncomp + 1}, nmax);\n\n"
              "    auto *prim = buffer.data(0);\n\n    std::fill(prim, prim + nmax, 0.0);\n\n")

    s += ("    // NOTE: the squared distances of the atom pairs are carried by the\n"
          "    // coordinates, so that they are formed once for the whole block instead of\n"
          "    // once for every combination of basis functions.\n\n"
          "    const auto *ab_2 = coordinates.data(6);\n\n"
          "    constexpr auto fpi = mathconst::pi_value();\n\n")
    s += "    // accumulate the integrals of each pair of primitives\n\n"
    s += ("    for (size_t i = 0; i < nprim_a; i++)\n    {\n        const auto aexp = a_exps[i];\n\n"
          "        const auto anorm = a_norms[i];\n\n"
          "        for (size_t j = 0; j < nprim_b; j++)\n        {\n"
          "            const auto ncols = dimensions[i * nprim_b + j];\n\n"
          "            if (ncols == 0) continue;\n\n"
          "            const auto fexp = aexp + b_exps[j];\n\n"
          "            const auto fmu = aexp * b_exps[j] / fexp;\n\n"
          "            const auto fovl = fpi / fexp;\n\n"
          "            auto ffact = anorm * b_norms[j] * fovl * std::sqrt(fovl);\n\n")
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
    s += ("            // NOTE: the kinetic energy of a pair of primitives is the overlap of the\n"
          f"            // pair times {2 * l + 3} mu less twice mu squared times the squared distance,\n"
          "            // so the factor of the loop carries the two terms of that bracket.\n\n"
          f"            const auto fkin = {'3.0' if l == 0 else f'{2 * l + 3}.0'} * fmu * ffact;\n\n"
          "            const auto fsqd = 2.0 * fmu * fmu * ffact;\n\n")
    s += ("            // NOTE: the row of the buffer and the row of the coordinates start at\n"
          "            // a cache line boundary, so the loop is vectorized with aligned loads\n"
          "            // and stores.\n\n"
          "#pragma omp simd aligned(prim, ab_2 : simd::cache_line_size())\n"
          "            for (size_t k = 0; k < ncols; k++)\n            {\n"
          "                prim[k] += (fkin - fsqd * ab_2[k]) * std::exp(-fmu * ab_2[k]);\n"
          "            }\n        }\n    }\n\n")

    if l == 0:
        s += ("    // NOTE: the atom pairs beyond the reach of every pair of primitives have no\n"
              "    // contribution and are set to zero.\n\n"
              "    std::copy(prim, prim + nmax, values);\n\n"
              "    std::fill(values + nmax, values + nvalues, 0.0);\n")
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
        s += "    for (size_t k = 0; k < nmax; k++)\n    {\n        const auto fval = prim[k];\n\n"
        for c in cs:
            s += f"        pb_{c}[k] = fval * ph_{c}[k];\n"
        s += "    }\n\n"
        s += ("    // NOTE: the values of an angular component are stored as one row of nvalues\n"
              "    // columns, and the atom pairs beyond the reach of every pair of primitives\n"
              "    // have no contribution and are set to zero.\n\n"
              f"    for (size_t m = 0; m < {ncomp}; m++)\n    {{\n"
              "        const auto *pb = buffer.data(m + 1);\n\n"
              "        std::copy(pb, pb + nmax, values + m * nvalues);\n\n"
              "        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);\n    }\n")
    s += "}\n\n}  // namespace simdkin\n"

    h = LIC_H + f"#ifndef SimdKineticEnergyRec{up}_hpp\n#define SimdKineticEnergyRec{up}_hpp\n\n"
    h += "#include <cstddef>\n\n"
    h += '#include "BasisFunction.hpp"\n#include "SimdMatrix.hpp"\n\nnamespace simdkin {  // simdkin namespace\n\n'
    h += (f"/// @brief Computes the kinetic energy integrals of the basis functions of angular\n"
          f"/// momenta {WORD[la]} and {WORD[lb]} over the atom pairs of a block.\n"
          "/// @param values The values of the combination of basis functions to compute.\n"
          "/// @param nvalues The number of atom pairs the combination reaches.\n"
          "/// @param bra The basis function on bra side.\n"
          "/// @param ket The basis function on ket side.\n")
    if l > 0:
        h += f"/// @param harmonics The solid harmonics of angular momentum {WORD[l]} of the vectors between the atoms.\n"
    h += ("/// @param coordinates The coordinates of the atom pairs.\n"
          "/// @param threshold The screening threshold of the integrals.\n")
    h += f"auto {sig};\n\n}}  // namespace simdkin\n\n#endif /* SimdKineticEnergyRec{up}_hpp */\n"
    return s, h


made = []
for la, lb in [(0, 0)] + [(0, l) for l in range(1, 7)] + [(l, 0) for l in range(1, 7)]:
    cpp, hpp = emit(la, lb)
    up = (LET[la] + LET[lb]).upper()
    io.open(DST + f'SimdKineticEnergyRec{up}.cpp', 'w').write(cpp)
    io.open(DST + f'SimdKineticEnergyRec{up}.hpp', 'w').write(hpp)
    made.append(up)
print("wrote:", " ".join(made))
