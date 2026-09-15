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


#include "SimdNuclearPotentialRecHI.hpp"

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
#include "SimdNuclearPotentialVrrRecSO.hpp"
#include "SimdNuclearPotentialVrrRecSP.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferDM.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferFL.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferGK.hpp"
#include "SimdTransferHI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransferPN.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_hi_nuclear_potential(double                    *values,
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
            false, std::string("compute_hi_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_hi_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 143 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 6266, 1373, 308, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 143 * nvalues, 0.0);

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

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 11, ncols,
                                                          fz, 6, p);

                compute_prim_sp_nuclear_potential_0(buffer, 20, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 23, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 26, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 29, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 32, 0, 3, 12, 13, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 35, 0, 3, 13, 14, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 38, 0, 3, 14, 15, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 41, 0, 3, 15, 16, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 44, 0, 3, 16, 17, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 47, 0, 3, 17, 18, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 50, 0, 3, 18, 19, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 53, 0, 3, 8, 20, 9, 23, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 59, 0, 3, 9, 23, 10, 26, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 65, 0, 3, 10, 26, 11, 29, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 71, 0, 3, 11, 29, 12, 32, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 77, 0, 3, 12, 32, 13, 35, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 83, 0, 3, 13, 35, 14, 38, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 89, 0, 3, 14, 38, 15, 41, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 95, 0, 3, 15, 41, 16, 44, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 101, 0, 3, 16, 44, 17, 47, ncols,
                                                    p);

                compute_prim_sd_nuclear_potential_0(buffer, 107, 0, 3, 17, 47, 18, 50, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 113, 0, 3, 20, 53, 23, 59, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 123, 0, 3, 23, 59, 26, 65, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 133, 0, 3, 26, 65, 29, 71, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 143, 0, 3, 29, 71, 32, 77, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 153, 0, 3, 32, 77, 35, 83, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 163, 0, 3, 35, 83, 38, 89, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 173, 0, 3, 38, 89, 41, 95, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 183, 0, 3, 41, 95, 44, 101, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 193, 0, 3, 44, 101, 47, 107, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 203, 0, 3, 53, 113, 59, 123, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 218, 0, 3, 59, 123, 65, 133, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 233, 0, 3, 65, 133, 71, 143, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 248, 0, 3, 71, 143, 77, 153, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 263, 0, 3, 77, 153, 83, 163, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 278, 0, 3, 83, 163, 89, 173, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 293, 0, 3, 89, 173, 95, 183, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 308, 0, 3, 95, 183, 101, 193, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 323, 0, 3, 113, 203, 123, 218, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 344, 0, 3, 123, 218, 133, 233, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 365, 0, 3, 133, 233, 143, 248, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 386, 0, 3, 143, 248, 153, 263, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 407, 0, 3, 153, 263, 163, 278, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 428, 0, 3, 163, 278, 173, 293, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 449, 0, 3, 173, 293, 183, 308, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 470, 0, 3, 203, 323, 218, 344, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 498, 0, 3, 218, 344, 233, 365, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 526, 0, 3, 233, 365, 248, 386, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 554, 0, 3, 248, 386, 263, 407, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 582, 0, 3, 263, 407, 278, 428, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 610, 0, 3, 278, 428, 293, 449, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 638, 0, 3, 323, 470, 344, 498, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 674, 0, 3, 344, 498, 365, 526, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 710, 0, 3, 365, 526, 386, 554, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 746, 0, 3, 386, 554, 407, 582, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 782, 0, 3, 407, 582, 428, 610, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 818, 0, 3, 470, 638, 498, 674, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 863, 0, 3, 498, 674, 526, 710, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 908, 0, 3, 526, 710, 554, 746, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 953, 0, 3, 554, 746, 582, 782, ncols,
                                                    p);

                compute_prim_sm_nuclear_potential_0(buffer, 998, 0, 3, 638, 818, 674, 863, ncols,
                                                    p);

                compute_prim_sm_nuclear_potential_0(buffer, 1053, 0, 3, 674, 863, 710, 908,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1108, 0, 3, 710, 908, 746, 953,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1163, 0, 3, 818, 998, 863, 1053,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1229, 0, 3, 863, 1053, 908, 1108,
                                                    ncols, p);

                compute_prim_so_nuclear_potential_0(buffer, 1295, 0, 3, 998, 1163, 1053, 1229,
                                                    ncols, p);

                simdfunc::contract_primitives(buffer, 1373, 470, 28, ncols);

                simdfunc::contract_primitives(buffer, 1401, 638, 36, ncols);

                simdfunc::contract_primitives(buffer, 1437, 818, 45, ncols);

                simdfunc::contract_primitives(buffer, 1482, 998, 55, ncols);

                simdfunc::contract_primitives(buffer, 1537, 1163, 66, ncols);

                simdfunc::contract_primitives(buffer, 1603, 1295, 78, ncols);
            }
        }
    }

    simdtrf::compute_hrr_pi(buffer, coordinates, 1681, 1373, 1401, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 1765, 1401, 1437, 1, nmax);

    simdtrf::compute_hrr_pl(buffer, coordinates, 1873, 1437, 1482, 1, nmax);

    simdtrf::compute_hrr_pm(buffer, coordinates, 2008, 1482, 1537, 1, nmax);

    simdtrf::compute_hrr_pn(buffer, coordinates, 2173, 1537, 1603, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 2371, 1681, 1765, 1, nmax);

    simdtrf::compute_hrr_dk(buffer, coordinates, 2539, 1765, 1873, 1, nmax);

    simdtrf::compute_hrr_dl(buffer, coordinates, 2755, 1873, 2008, 1, nmax);

    simdtrf::compute_hrr_dm(buffer, coordinates, 3025, 2008, 2173, 1, nmax);

    simdtrf::compute_hrr_fi(buffer, coordinates, 3355, 2371, 2539, 1, nmax);

    simdtrf::compute_hrr_fk(buffer, coordinates, 3635, 2539, 2755, 1, nmax);

    simdtrf::compute_hrr_fl(buffer, coordinates, 3995, 2755, 3025, 1, nmax);

    simdtrf::compute_hrr_gi(buffer, coordinates, 4445, 3355, 3635, 1, nmax);

    simdtrf::compute_hrr_gk(buffer, coordinates, 4865, 3635, 3995, 1, nmax);

    simdtrf::compute_hrr_hi(buffer, coordinates, 5405, 4445, 4865, 1, nmax);

    simdtrf::transform_i_inner(buffer, 5993, 5405, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 5993, 13, nmax);

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
