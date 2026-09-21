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


#include "SimdNuclearPotentialRecHH.hpp"

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
#include "SimdNuclearPotentialVrrRecSI.hpp"
#include "SimdNuclearPotentialVrrRecSK.hpp"
#include "SimdNuclearPotentialVrrRecSL.hpp"
#include "SimdNuclearPotentialVrrRecSM.hpp"
#include "SimdNuclearPotentialVrrRecSN.hpp"
#include "SimdNuclearPotentialVrrRecSP.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferGH.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferHH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransformH.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_hh_nuclear_potential(double                    *values,
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
            false, std::string("compute_hh_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_hh_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 121 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 4852, 1009, 251, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 121 * nvalues, 0.0);

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

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 10, ncols,
                                                          fz, 6, p);

                compute_prim_sp_nuclear_potential_0(buffer, 19, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 22, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 25, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 28, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 31, 0, 3, 12, 13, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 34, 0, 3, 13, 14, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 37, 0, 3, 14, 15, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 40, 0, 3, 15, 16, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 43, 0, 3, 16, 17, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 46, 0, 3, 17, 18, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 49, 0, 3, 8, 19, 9, 22, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 55, 0, 3, 9, 22, 10, 25, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 61, 0, 3, 10, 25, 11, 28, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 67, 0, 3, 11, 28, 12, 31, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 73, 0, 3, 12, 31, 13, 34, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 79, 0, 3, 13, 34, 14, 37, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 85, 0, 3, 14, 37, 15, 40, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 91, 0, 3, 15, 40, 16, 43, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 97, 0, 3, 16, 43, 17, 46, ncols, p);

                compute_prim_sf_nuclear_potential_0(buffer, 103, 0, 3, 19, 49, 22, 55, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 113, 0, 3, 22, 55, 25, 61, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 123, 0, 3, 25, 61, 28, 67, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 133, 0, 3, 28, 67, 31, 73, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 143, 0, 3, 31, 73, 34, 79, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 153, 0, 3, 34, 79, 37, 85, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 163, 0, 3, 37, 85, 40, 91, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 173, 0, 3, 40, 91, 43, 97, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 183, 0, 3, 49, 103, 55, 113, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 198, 0, 3, 55, 113, 61, 123, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 213, 0, 3, 61, 123, 67, 133, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 228, 0, 3, 67, 133, 73, 143, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 243, 0, 3, 73, 143, 79, 153, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 258, 0, 3, 79, 153, 85, 163, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 273, 0, 3, 85, 163, 91, 173, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 288, 0, 3, 103, 183, 113, 198, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 309, 0, 3, 113, 198, 123, 213, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 330, 0, 3, 123, 213, 133, 228, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 351, 0, 3, 133, 228, 143, 243, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 372, 0, 3, 143, 243, 153, 258, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 393, 0, 3, 153, 258, 163, 273, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 414, 0, 3, 183, 288, 198, 309, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 442, 0, 3, 198, 309, 213, 330, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 470, 0, 3, 213, 330, 228, 351, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 498, 0, 3, 228, 351, 243, 372, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 526, 0, 3, 243, 372, 258, 393, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 554, 0, 3, 288, 414, 309, 442, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 590, 0, 3, 309, 442, 330, 470, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 626, 0, 3, 330, 470, 351, 498, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 662, 0, 3, 351, 498, 372, 526, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 698, 0, 3, 414, 554, 442, 590, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 743, 0, 3, 442, 590, 470, 626, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 788, 0, 3, 470, 626, 498, 662, ncols,
                                                    p);

                compute_prim_sm_nuclear_potential_0(buffer, 833, 0, 3, 554, 698, 590, 743, ncols,
                                                    p);

                compute_prim_sm_nuclear_potential_0(buffer, 888, 0, 3, 590, 743, 626, 788, ncols,
                                                    p);

                compute_prim_sn_nuclear_potential_0(buffer, 943, 0, 3, 698, 833, 743, 888, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 1009, 288, 21, ncols);

                simdfunc::contract_primitives(buffer, 1030, 414, 28, ncols);

                simdfunc::contract_primitives(buffer, 1058, 554, 36, ncols);

                simdfunc::contract_primitives(buffer, 1094, 698, 45, ncols);

                simdfunc::contract_primitives(buffer, 1139, 833, 55, ncols);

                simdfunc::contract_primitives(buffer, 1194, 943, 66, ncols);
            }
        }
    }

    simdtrf::compute_hrr_ph(buffer, coordinates, 1260, 1009, 1030, 1, nmax);

    simdtrf::compute_hrr_pi(buffer, coordinates, 1323, 1030, 1058, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 1407, 1058, 1094, 1, nmax);

    simdtrf::compute_hrr_pl(buffer, coordinates, 1515, 1094, 1139, 1, nmax);

    simdtrf::compute_hrr_pm(buffer, coordinates, 1650, 1139, 1194, 1, nmax);

    simdtrf::compute_hrr_dh(buffer, coordinates, 1815, 1260, 1323, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 1941, 1323, 1407, 1, nmax);

    simdtrf::compute_hrr_dk(buffer, coordinates, 2109, 1407, 1515, 1, nmax);

    simdtrf::compute_hrr_dl(buffer, coordinates, 2325, 1515, 1650, 1, nmax);

    simdtrf::compute_hrr_fh(buffer, coordinates, 2595, 1815, 1941, 1, nmax);

    simdtrf::compute_hrr_fi(buffer, coordinates, 2805, 1941, 2109, 1, nmax);

    simdtrf::compute_hrr_fk(buffer, coordinates, 3085, 2109, 2325, 1, nmax);

    simdtrf::compute_hrr_gh(buffer, coordinates, 3445, 2595, 2805, 1, nmax);

    simdtrf::compute_hrr_gi(buffer, coordinates, 3760, 2805, 3085, 1, nmax);

    simdtrf::compute_hrr_hh(buffer, coordinates, 4180, 3445, 3760, 1, nmax);

    simdtrf::transform_h_inner(buffer, 4621, 4180, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 4621, 11, nmax);

    for (size_t m = 0; m < 121; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
