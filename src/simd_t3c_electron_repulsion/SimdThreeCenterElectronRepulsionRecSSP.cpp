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



#include "SimdThreeCenterElectronRepulsionRecSSP.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdBoysFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdVariableMatrix.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ssp_electron_repulsion(double               *values,
                               const size_t          npairs,
                               const size_t          natoms,
                               const size_t          iatom,
                               const CBasisFunction &a_function,
                               const CBasisFunction &b_function,
                               const CBasisFunction &c_function,
                               const CSimdMatrix    &ab_harmonics,
                               const CSimdMatrix    &bc_harmonics,
                               const CSimdMatrix    &ab_coordinates,
                               const CSimdMatrix    &bc_coordinates,
                               const double          threshold) -> void
{
    if ((a_function.get_angular_momentum() != 0) || (b_function.get_angular_momentum() != 0) ||
        (c_function.get_angular_momentum() != 1))
    {
        errors::assertMsgCritical(
            false,
            std::string("SimdThreeCenterElectronRepulsionRecSSP.compute_ssp_electron_repulsion: Basis functions must be of angular momenta zero, zero and one"));
    }

    if ((ab_harmonics.number_of_rows() != 3) || (bc_harmonics.number_of_rows() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdThreeCenterElectronRepulsionRecSSP.compute_ssp_electron_repulsion: Harmonics must have 3 rows"));
    }

    if (npairs > ab_coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdThreeCenterElectronRepulsionRecSSP.compute_ssp_electron_repulsion: Number of atom pairs exceeds coordinates"));
    }

    if (iatom >= natoms)
    {
        errors::assertMsgCritical(
            false, std::string("SimdThreeCenterElectronRepulsionRecSSP.compute_ssp_electron_repulsion: Index of atom on c side is out of range"));
    }

    if (npairs == 0) return;

    const auto &a_exps = a_function.exponents();

    const auto &b_exps = b_function.exponents();

    const auto &c_exps = c_function.exponents();

    const auto &a_norms = a_function.normalization_factors();

    const auto &b_norms = b_function.normalization_factors();

    const auto &c_norms = c_function.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprim_c = c_exps.size();

    const auto nprims = nprim_a * nprim_b * nprim_c;

    // NOTE: the triples of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(a_function,
                                                             b_function,
                                                             c_function,
                                                             npairs,
                                                             ab_coordinates,
                                                             screenfunc::three_center_electron_repulsion_primitive_bound,
                                                             threshold / static_cast<double>(nprims));

    const auto nmax = *std::ranges::max_element(dimensions);

    // NOTE: the values of one atom on c side are contiguous over the atom pairs,
    // the atoms are npairs apart and the angular components of the atoms on c
    // side are npairs times the number of those atoms apart, as the sparsity
    // pattern lays them out.

    const auto stride = natoms * npairs;

    auto *slice_m1 = values + iatom * npairs;
    auto *slice_z0 = slice_m1 + stride;
    auto *slice_p1 = slice_z0 + stride;

    if (nmax == 0)
    {
        std::fill(slice_m1, slice_m1 + npairs, 0.0);
        std::fill(slice_z0, slice_z0 + npairs, 0.0);
        std::fill(slice_p1, slice_p1 + npairs, 0.0);

        return;
    }

    const auto *a_x = ab_coordinates.data(0);
    const auto *a_y = ab_coordinates.data(1);
    const auto *a_z = ab_coordinates.data(2);

    const auto *ab_2 = ab_coordinates.data(6);

    const auto *b_x = bc_coordinates.data(0);
    const auto *b_y = bc_coordinates.data(1);
    const auto *b_z = bc_coordinates.data(2);

    const auto *c_x = bc_coordinates.data(3);
    const auto *c_y = bc_coordinates.data(4);
    const auto *c_z = bc_coordinates.data(5);

    // NOTE: the exponential of the pair on a and b sides and the squared distance
    // from their product center to the atom on c side depend on the exponents of
    // that pair alone, the ket exponent entering neither, so both are formed once
    // for a pair of primitives and read by every primitive on c side.

    auto factors = CSimdMatrix(2, nmax);

    auto *e_ab = factors.data(0);

    auto *pc_2 = factors.data(1);

    // NOTE: one row accumulates the factor which multiplies the harmonics of the
    // atom pairs and the other the factor which multiplies the harmonics of the
    // atoms on c side, as the two bidegrees of the addition theorem carry
    // different powers of the exponents and cannot share an accumulator.

    auto buffer = CSimdMatrix(2, nmax);

    buffer.zero();

    auto *acc_ab = buffer.data(0);

    auto *acc_bc = buffer.data(1);

    constexpr auto fpi = mathconst::pi_value();

    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto bexp = b_exps[j];

            // NOTE: the widest of the triples of this pair of primitives is
            // searched for rather than assumed to be the last, as the bound of a
            // triple carries its prefactor as well as its decay.

            const auto first = dimensions.begin() + static_cast<long>((i * nprim_b + j) * nprim_c);

            const auto npair_max = *std::ranges::max_element(first, first + static_cast<long>(nprim_c));

            if (npair_max == 0) continue;

            const auto pexp = aexp + bexp;

            const auto fmu = aexp * bexp / pexp;

            const auto frp = 1.0 / pexp;

#pragma omp simd aligned(e_ab, pc_2, ab_2, a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z : simd::cache_line_size())
            for (size_t k = 0; k < npair_max; k++)
            {
                e_ab[k] = std::exp(-fmu * ab_2[k]);

                const auto p_x = frp * (aexp * (a_x[k] - c_x[k]) + bexp * (b_x[k] - c_x[k]));

                const auto p_y = frp * (aexp * (a_y[k] - c_y[k]) + bexp * (b_y[k] - c_y[k]));

                const auto p_z = frp * (aexp * (a_z[k] - c_z[k]) + bexp * (b_z[k] - c_z[k]));

                pc_2[k] = p_x * p_x + p_y * p_y + p_z * p_z;
            }

            // NOTE: the Boys function of every primitive on c side of this pair
            // is computed by one call, which fills the orders zero and one of
            // every row. The integrals need the order one alone, and the order
            // zero is formed on the way to it by the recursion.

            auto boys = CSimdVariableMatrix(std::vector<size_t>(first, first + static_cast<long>(nprim_c)), 3);

            for (size_t k = 0; k < nprim_c; k++)
            {
                const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                if (ncols == 0) continue;

                const auto frho = pexp * c_exps[k] / (pexp + c_exps[k]);

                auto *bargs = boys.data(0, k);

#pragma omp simd aligned(bargs, pc_2 : simd::cache_line_size())
                for (size_t l = 0; l < ncols; l++)
                {
                    bargs[l] = frho * pc_2[l];
                }
            }

            simdfunc::compute_boys_function(boys);

            for (size_t k = 0; k < nprim_c; k++)
            {
                const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                if (ncols == 0) continue;

                const auto cexp = c_exps[k];

                const auto qexp = pexp + cexp;

                // NOTE: the auxiliary integral of order one is the integral of
                // three S type primitives with the Boys function of order one in
                // place of the order zero. The bidegree of the addition theorem
                // weighs the harmonics of the atom pairs by the exponent on a
                // side and those of the atoms on c side by the exponent of the
                // pair, both over the total exponent.

                const auto faux = fcoul * anorm * b_norms[j] * c_norms[k] / (pexp * cexp * std::sqrt(qexp) * qexp);

                const auto fab = faux * aexp;

                const auto fbc = faux * pexp;

                const auto *bvals = boys.data(2, k);

#pragma omp simd aligned(acc_ab, acc_bc, e_ab, bvals : simd::cache_line_size())
                for (size_t l = 0; l < ncols; l++)
                {
                    const auto fval = e_ab[l] * bvals[l];

                    acc_ab[l] += fab * fval;

                    acc_bc[l] += fbc * fval;
                }
            }
        }
    }

    // NOTE: the harmonics of angular momentum one are the components of the
    // vector taken in the order of the spherical components, so the row of index
    // m + 1 holds the harmonic of order m on both sides.

    const auto *hab_m1 = ab_harmonics.data(0);
    const auto *hab_z0 = ab_harmonics.data(1);
    const auto *hab_p1 = ab_harmonics.data(2);

    const auto *hbc_m1 = bc_harmonics.data(0);
    const auto *hbc_z0 = bc_harmonics.data(1);
    const auto *hbc_p1 = bc_harmonics.data(2);

    auto components = CSimdMatrix(3, nmax);

    auto *out_m1 = components.data(0);
    auto *out_z0 = components.data(1);
    auto *out_p1 = components.data(2);

#pragma omp simd aligned(out_m1, out_z0, out_p1, acc_ab, acc_bc, hab_m1, hab_z0, hab_p1, hbc_m1, hbc_z0, hbc_p1 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        out_m1[k] = acc_ab[k] * hab_m1[k] + acc_bc[k] * hbc_m1[k];

        out_z0[k] = acc_ab[k] * hab_z0[k] + acc_bc[k] * hbc_z0[k];

        out_p1[k] = acc_ab[k] * hab_p1[k] + acc_bc[k] * hbc_p1[k];
    }

    // NOTE: the atom pairs beyond the reach of every triple of primitives have no
    // contribution and are set to zero.

    std::copy(out_m1, out_m1 + nmax, slice_m1);
    std::copy(out_z0, out_z0 + nmax, slice_z0);
    std::copy(out_p1, out_p1 + nmax, slice_p1);

    std::fill(slice_m1 + nmax, slice_m1 + npairs, 0.0);
    std::fill(slice_z0 + nmax, slice_z0 + npairs, 0.0);
    std::fill(slice_p1 + nmax, slice_p1 + npairs, 0.0);
}

}  // namespace simdt3ceri
