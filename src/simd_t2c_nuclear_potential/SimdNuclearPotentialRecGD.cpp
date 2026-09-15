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


#include "SimdNuclearPotentialRecGD.hpp"

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

#include "SimdNuclearPotentialVrrRecDS.hpp"
#include "SimdNuclearPotentialVrrRecFS.hpp"
#include "SimdNuclearPotentialVrrRecGS.hpp"
#include "SimdNuclearPotentialVrrRecHS.hpp"
#include "SimdNuclearPotentialVrrRecIS.hpp"
#include "SimdNuclearPotentialVrrRecPS.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_gd_nuclear_potential(double                    *values,
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
            false, std::string("compute_gd_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_gd_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 45 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 555, 218, 64, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 45 * nvalues, 0.0);

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

            const auto fa = -b_exps[j] / p;

            const auto fc = b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pair_exponent(buffer, coordinates, 6, ncols, mu);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 6, ncols,
                                                          fz, 6, p);

                compute_prim_ps_nuclear_potential_0(buffer, 15, 0, 3, 8, 9, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 18, 0, 3, 9, 10, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 21, 0, 3, 10, 11, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 24, 0, 3, 11, 12, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 27, 0, 3, 12, 13, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 30, 0, 3, 13, 14, ncols);

                compute_prim_ds_nuclear_potential_0(buffer, 33, 0, 3, 8, 9, 15, 18, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 39, 0, 3, 9, 10, 18, 21, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 45, 0, 3, 10, 11, 21, 24, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 51, 0, 3, 11, 12, 24, 27, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 57, 0, 3, 12, 13, 27, 30, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 63, 0, 3, 15, 18, 33, 39, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 73, 0, 3, 18, 21, 39, 45, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 83, 0, 3, 21, 24, 45, 51, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 93, 0, 3, 24, 27, 51, 57, ncols, p);

                compute_prim_gs_nuclear_potential_0(buffer, 103, 0, 3, 33, 39, 63, 73, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 118, 0, 3, 39, 45, 73, 83, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 133, 0, 3, 45, 51, 83, 93, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 148, 0, 3, 63, 73, 103, 118, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 169, 0, 3, 73, 83, 118, 133, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 190, 0, 3, 103, 118, 148, 169, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 218, 103, 15, ncols);

                simdfunc::contract_primitives(buffer, 233, 148, 21, ncols);

                simdfunc::contract_primitives(buffer, 254, 190, 28, ncols);
            }
        }
    }

    simdtrf::compute_hrr_gp(buffer, coordinates, 282, 218, 233, 1, nmax);

    simdtrf::compute_hrr_hp(buffer, coordinates, 327, 233, 254, 1, nmax);

    simdtrf::compute_hrr_gd(buffer, coordinates, 390, 282, 327, 1, nmax);

    simdtrf::transform_d_inner(buffer, 480, 390, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 480, 5, nmax);

    for (size_t m = 0; m < 45; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
