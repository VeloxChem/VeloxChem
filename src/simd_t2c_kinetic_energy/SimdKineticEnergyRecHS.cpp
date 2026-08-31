//
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



#include "SimdKineticEnergyRecHS.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hs_kinetic_energy(double               *values,
                           const size_t          nvalues,
                           const CBasisFunction &bra,
                           const CBasisFunction &ket,
                           const CSimdMatrix    &harmonics,
                           const CSimdMatrix    &coordinates,
                           const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 0))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hs_kinetic_energy: Basis functions must be of angular momenta five and zero"));
    }

    if (harmonics.number_of_rows() != 11)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hs_kinetic_energy: Harmonics must have 11 rows"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_hs_kinetic_energy: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_kinetic_energy_primitive_bound,
        threshold / static_cast<double>(nprims));

    // NOTE: the buffer spans the atom pairs reached by the pair of primitives
    // reaching furthest, which is searched for rather than assumed. The
    // primitives are sorted by descending exponent, but the bound of a pair of
    // primitives carries their prefactor as well as their decay, so a tighter
    // pair with a larger prefactor reaches further than a more diffuse pair with
    // a smaller one, and the last pair is not always the furthest reaching.

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 11 * nvalues, 0.0);

        return;
    }

    // NOTE: the first row accumulates the factor which the angular components
    // share, and the remaining rows hold the integrals of the angular components.

    auto buffer = CSimdMatrix(12, nmax);

    auto *prim = buffer.data(0);

    std::fill(prim, prim + nmax, 0.0);

    // NOTE: the squared distances of the atom pairs are carried by the
    // coordinates, so that they are formed once for the whole block instead of
    // once for every combination of basis functions.

    const auto *ab_2 = coordinates.data(6);

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the integrals of each pair of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto fexp = aexp + b_exps[j];

            const auto fmu = aexp * b_exps[j] / fexp;

            const auto fovl = fpi / fexp;

            auto ffact = anorm * b_norms[j] * fovl * std::sqrt(fovl);

            // NOTE: the ratio is raised to the angular momentum by repeated
            // multiplication, as the angular momentum is small.

            const auto frat = -b_exps[j] / fexp;

            ffact *= frat;

            ffact *= frat;

            ffact *= frat;

            ffact *= frat;

            ffact *= frat;

            // NOTE: the kinetic energy of a pair of primitives is the overlap of the
            // pair times 13 mu less twice mu squared times the squared distance,
            // so the factor of the loop carries the two terms of that bracket.

            const auto fkin = 13.0 * fmu * ffact;

            const auto fsqd = 2.0 * fmu * fmu * ffact;

            // NOTE: the row of the buffer and the row of the coordinates start at
            // a cache line boundary, so the loop is vectorized with aligned loads
            // and stores.

#pragma omp simd aligned(prim, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                prim[k] += (fkin - fsqd * ab_2[k]) * std::exp(-fmu * ab_2[k]);
            }
        }
    }

    // NOTE: the integral of an angular component is the accumulated factor
    // times the solid harmonic of that component.

    const auto *ph_m5 = harmonics.data(0);
    const auto *ph_m4 = harmonics.data(1);
    const auto *ph_m3 = harmonics.data(2);
    const auto *ph_m2 = harmonics.data(3);
    const auto *ph_m1 = harmonics.data(4);
    const auto *ph_0 = harmonics.data(5);
    const auto *ph_p1 = harmonics.data(6);
    const auto *ph_p2 = harmonics.data(7);
    const auto *ph_p3 = harmonics.data(8);
    const auto *ph_p4 = harmonics.data(9);
    const auto *ph_p5 = harmonics.data(10);

    auto *pb_m5 = buffer.data(1);
    auto *pb_m4 = buffer.data(2);
    auto *pb_m3 = buffer.data(3);
    auto *pb_m2 = buffer.data(4);
    auto *pb_m1 = buffer.data(5);
    auto *pb_0 = buffer.data(6);
    auto *pb_p1 = buffer.data(7);
    auto *pb_p2 = buffer.data(8);
    auto *pb_p3 = buffer.data(9);
    auto *pb_p4 = buffer.data(10);
    auto *pb_p5 = buffer.data(11);

#pragma omp simd aligned(prim, ph_m5, ph_m4, ph_m3, ph_m2, ph_m1, ph_0, ph_p1, ph_p2, ph_p3, ph_p4, ph_p5, pb_m5, pb_m4, pb_m3, pb_m2, pb_m1, pb_0, pb_p1, pb_p2, pb_p3, pb_p4, pb_p5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto fval = prim[k];

        pb_m5[k] = fval * ph_m5[k];
        pb_m4[k] = fval * ph_m4[k];
        pb_m3[k] = fval * ph_m3[k];
        pb_m2[k] = fval * ph_m2[k];
        pb_m1[k] = fval * ph_m1[k];
        pb_0[k] = fval * ph_0[k];
        pb_p1[k] = fval * ph_p1[k];
        pb_p2[k] = fval * ph_p2[k];
        pb_p3[k] = fval * ph_p3[k];
        pb_p4[k] = fval * ph_p4[k];
        pb_p5[k] = fval * ph_p5[k];
    }

    // NOTE: the values of an angular component are stored as one row of nvalues
    // columns, and the atom pairs beyond the reach of every pair of primitives
    // have no contribution and are set to zero.

    for (size_t m = 0; m < 11; m++)
    {
        const auto *pb = buffer.data(m + 1);

        std::copy(pb, pb + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
