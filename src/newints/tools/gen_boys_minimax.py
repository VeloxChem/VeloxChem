#!/usr/bin/env python3
"""Generate src/newints/BoysMinimaxCoefficients.hpp from the published
rational-minimax coefficients of the Boys functions.

Source: R. Vikhamar-Sandberg and M. Repisky, "Efficient Boys function
evaluation using minimax approximation", arXiv:2512.10059v3 (2025),
ancillary files anc/coefficients.txt (region-A tables F0..F32) and
anc/minimax.c / anc/constants.h (region-B table and region boundaries).

This script is the single source of truth for the generated header -- do NOT
hand-edit BoysMinimaxCoefficients.hpp; change this script (or the input data)
and regenerate:

    python3 src/newints/tools/gen_boys_minimax.py

The coefficients are published numerical constants; the algorithm they feed is
implemented in src/newints/BoysFunction.hpp.
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, "boys_minimax_coefficients.txt")
OUTPUT = os.path.normpath(os.path.join(HERE, "..", "BoysMinimaxCoefficients.hpp"))

MAX_ORDER = 32  # F0 .. F32

# Region boundaries and sqrt(pi)/2 from anc/constants.h.
X0 = "11.899848152108484"
X1 = "28.989337738820740"
HALF_SQRT_PI = "0.886226925452758014"

# Region-B rational minimax (F0 on [x0, x1]) from anc/minimax.c (c*/d*),
# stored ascending in degree; the denominator is monic so a trailing 1.0 is
# appended to match the region-A storage convention.
RB_NUM = [
    "5.74537531702047552E+7",
    "2.73330925890901898E+6",
    "7.52922255805293133E+04",
    "2.33846894861346960E+05",
    "8.34841284469484906E+03",
    "3.90892739018191431E+01",
]
RB_DEN = [
    "4.79893571439451030E+07",
    "3.04808499107506708E+07",
    "-1.66693114610725015E+06",
    "5.63505368535215625E+05",
    "6.39702496081641495E+04",
    "8.53693546919731980E+02",
    "1.00000000000000000E+00",
]


def parse_region_a(path):
    """Return {order: (num_list, den_list)} of coefficient strings, ascending
    in polynomial degree (constant term first, monic leading term last).

    The published file mislabels the numerator/denominator boundary for some
    orders (the denominator constant term is printed under the numerator
    label). We therefore ignore the labels and recover the true split from the
    exact value F_k(0) = num[0]/den[0] = 1/(2k+1): the denominator constant
    term is the coefficient closest to (2k+1)*num[0]."""
    blocks = {}  # order -> list of (raw_string, float_value)
    order = None
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*F(\d+)\s+(?:numerator|denominator)", line)
            if m:
                order = int(m.group(1))
                blocks.setdefault(order, [])
                continue
            s = line.strip()
            if s and order is not None:
                blocks[order].append((s, float(s)))

    tables = {}
    for k, items in blocks.items():
        raws = [r for r, _ in items]
        vals = [v for _, v in items]
        assert abs(vals[-1] - 1.0) < 1e-12, f"F{k}: denominator not monic ({vals[-1]})"
        target = (2 * k + 1) * vals[0]
        split = min(range(1, len(vals) - 1), key=lambda s: abs(vals[s] / target - 1.0))
        num, den = raws[:split], raws[split:]
        r0 = vals[0] / vals[split]
        assert abs(r0 * (2 * k + 1) - 1.0) < 1e-9, f"F{k}: bad split, r(0)={r0}"
        tables[k] = (num, den)
    return tables


def fmt_array(name, coeffs):
    body = ",\n    ".join(coeffs)
    return (
        f"    static constexpr std::array<double, {len(coeffs)}> {name} = {{\n"
        f"    {body}}};"
    )


def main():
    tables = parse_region_a(INPUT)
    assert sorted(tables) == list(range(MAX_ORDER + 1)), "missing orders"

    out = []
    out.append("""//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
""")
    out.append("// GENERATED FILE -- do not edit by hand.")
    out.append("// Regenerate with: python3 src/newints/tools/gen_boys_minimax.py")
    out.append("//")
    out.append("// Rational-minimax coefficients of the Boys functions from")
    out.append("// R. Vikhamar-Sandberg and M. Repisky, arXiv:2512.10059v3 (2025),")
    out.append("// ancillary files anc/coefficients.txt, anc/minimax.c, anc/constants.h.")
    out.append("//")
    out.append("// Each rational approximation is r(x) = num(x) / den(x) with the")
    out.append("// coefficients stored ascending in degree (constant term first, monic")
    out.append("// leading term last), to be evaluated by Horner's scheme.")
    out.append("")
    out.append("#ifndef newints_BoysMinimaxCoefficients_hpp")
    out.append("#define newints_BoysMinimaxCoefficients_hpp")
    out.append("")
    out.append("#include <array>")
    out.append("#include <cstddef>")
    out.append("")
    out.append("namespace newints::boys_detail {")
    out.append("")
    out.append("/// @brief Maximum supported Boys function order (F0 .. F_max_order).")
    out.append(f"inline constexpr int max_order = {MAX_ORDER};")
    out.append("")
    out.append("/// @brief Start of the middle region B; below this, region A (downward recursion).")
    out.append(f"inline constexpr double region_b_start = {X0};")
    out.append("")
    out.append("/// @brief Start of the asymptotic region C; above this, the asymptotic form.")
    out.append(f"inline constexpr double region_c_start = {X1};")
    out.append("")
    out.append("/// @brief sqrt(pi) / 2, the leading factor of the asymptotic F0.")
    out.append(f"inline constexpr double half_sqrt_pi = {HALF_SQRT_PI};")
    out.append("")
    out.append("/// @brief Region-B rational minimax for F0 on [region_b_start, region_c_start).")
    out.append(fmt_array("region_b_num", RB_NUM))
    out.append(fmt_array("region_b_den", RB_DEN))
    out.append("")
    out.append("/// @brief Region-A rational minimax tables for F_K on [0, region_b_start).")
    out.append("/// Selected at compile time by the top Boys order K.")
    out.append("template <int K>")
    out.append("struct region_a;")
    out.append("")
    for k in range(MAX_ORDER + 1):
        num, den = tables[k]
        out.append(f"template <>")
        out.append(f"struct region_a<{k}>")
        out.append("{")
        out.append(fmt_array("num", num))
        out.append(fmt_array("den", den))
        out.append("};")
        out.append("")
    out.append("}  // namespace newints::boys_detail")
    out.append("")
    out.append("#endif /* newints_BoysMinimaxCoefficients_hpp */")

    with open(OUTPUT, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"wrote {OUTPUT}")
    print(f"orders F0..F{MAX_ORDER}: "
          + ", ".join(f"F{k}({len(tables[k][0])}/{len(tables[k][1])})" for k in (0, 1, 12, 32)))


if __name__ == "__main__":
    main()
