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


#include "SimdNuclearPotentialRecPG.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdNuclearPotentialVrrRecSD.hpp"
#include "SimdNuclearPotentialVrrRecSF.hpp"
#include "SimdNuclearPotentialVrrRecSG.hpp"
#include "SimdNuclearPotentialVrrRecSH.hpp"
#include "SimdNuclearPotentialVrrRecSP.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_pg_nuclear_potential(double                    *values,
                             const size_t               nvalues,
                             const CBasisFunction      &bra,
                             const CBasisFunction      &ket,
                             const CSimdMatrix         &coordinates,
                             const std::vector<double> &charges,
                             const std::vector<double> &points,
                             CSimdMatrix               &buffer,
                             const double               threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_pg_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_pg_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 27 * nvalues, 0.0);

        return;
    }

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number and by the number of charges, as every
    // integral is a sum over both and the error of a sum is bounded by the
    // number of its terms.

    const auto terms = static_cast<double>(nprims * charges.size());

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_nuclear_potential_primitive_bound, threshold / terms);

    const auto nmax = simdfunc::prepare_buffer(buffer, 242, 134, 36, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 27 * nvalues, 0.0);

        return;
    }

    errors::assertMsgCritical(dimensions.size() == nprim_a * nprim_b,
                              std::string("Dimensions do not match the pairs of primitives"));

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto fnpot = 2.0 * mathconst::pi_value() / p * a_norms[i] * b_norms[j];

            const auto fb = a_exps[i] / p;

            const auto fc = b_exps[j] / p;

            simdfunc::compute_pb(buffer, coordinates, 0, ncols, fb);

            simdfunc::compute_pair_exponent(buffer, coordinates, 6, ncols, mu);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 5, ncols,
                                                          fz, 6, p);

                compute_prim_sp_nuclear_potential_0(buffer, 14, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 17, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 20, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 23, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 26, 0, 3, 12, 13, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 29, 0, 3, 8, 14, 9, 17, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 35, 0, 3, 9, 17, 10, 20, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 41, 0, 3, 10, 20, 11, 23, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 47, 0, 3, 11, 23, 12, 26, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 53, 0, 3, 14, 29, 17, 35, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 63, 0, 3, 17, 35, 20, 41, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 73, 0, 3, 20, 41, 23, 47, ncols, p);

                compute_prim_sg_nuclear_potential_0(buffer, 83, 0, 3, 29, 53, 35, 63, ncols, p);

                compute_prim_sg_nuclear_potential_0(buffer, 98, 0, 3, 35, 63, 41, 73, ncols, p);

                compute_prim_sh_nuclear_potential_0(buffer, 113, 0, 3, 53, 83, 63, 98, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 134, 83, 15, ncols);

                simdfunc::contract_primitives(buffer, 149, 113, 21, ncols);
            }
        }
    }

    simdtrf::compute_hrr_pg(buffer, coordinates, 170, 134, 149, 1, nmax);

    simdtrf::transform_g_inner(buffer, 215, 170, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 215, 9, nmax);

    for (size_t m = 0; m < 27; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
