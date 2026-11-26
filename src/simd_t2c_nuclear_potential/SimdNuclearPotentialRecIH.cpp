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


#include "SimdNuclearPotentialRecIH.hpp"

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
#include "SimdNuclearPotentialVrrRecKS.hpp"
#include "SimdNuclearPotentialVrrRecLS.hpp"
#include "SimdNuclearPotentialVrrRecMS.hpp"
#include "SimdNuclearPotentialVrrRecNS.hpp"
#include "SimdNuclearPotentialVrrRecOS.hpp"
#include "SimdNuclearPotentialVrrRecPS.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIH.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_ih_nuclear_potential(double                    *values,
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
            false, std::string("compute_ih_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_ih_nuclear_potential: Expecting three coordinates for each charge"));

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 6301, 1373, 308, dimensions);

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

            const auto fa = -b_exps[j] / p;

            const auto fc = b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pair_exponent(buffer, coordinates, 6, ncols, mu);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 11, ncols,
                                                          fz, 6, p);

                compute_prim_ps_nuclear_potential_0(buffer, 20, 0, 3, 8, 9, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 23, 0, 3, 9, 10, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 26, 0, 3, 10, 11, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 29, 0, 3, 11, 12, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 32, 0, 3, 12, 13, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 35, 0, 3, 13, 14, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 38, 0, 3, 14, 15, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 41, 0, 3, 15, 16, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 44, 0, 3, 16, 17, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 47, 0, 3, 17, 18, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 50, 0, 3, 18, 19, ncols);

                compute_prim_ds_nuclear_potential_0(buffer, 53, 0, 3, 8, 9, 20, 23, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 59, 0, 3, 9, 10, 23, 26, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 65, 0, 3, 10, 11, 26, 29, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 71, 0, 3, 11, 12, 29, 32, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 77, 0, 3, 12, 13, 32, 35, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 83, 0, 3, 13, 14, 35, 38, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 89, 0, 3, 14, 15, 38, 41, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 95, 0, 3, 15, 16, 41, 44, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 101, 0, 3, 16, 17, 44, 47, ncols,
                                                    p);

                compute_prim_ds_nuclear_potential_0(buffer, 107, 0, 3, 17, 18, 47, 50, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 113, 0, 3, 20, 23, 53, 59, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 123, 0, 3, 23, 26, 59, 65, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 133, 0, 3, 26, 29, 65, 71, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 143, 0, 3, 29, 32, 71, 77, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 153, 0, 3, 32, 35, 77, 83, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 163, 0, 3, 35, 38, 83, 89, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 173, 0, 3, 38, 41, 89, 95, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 183, 0, 3, 41, 44, 95, 101, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 193, 0, 3, 44, 47, 101, 107, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 203, 0, 3, 53, 59, 113, 123, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 218, 0, 3, 59, 65, 123, 133, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 233, 0, 3, 65, 71, 133, 143, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 248, 0, 3, 71, 77, 143, 153, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 263, 0, 3, 77, 83, 153, 163, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 278, 0, 3, 83, 89, 163, 173, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 293, 0, 3, 89, 95, 173, 183, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 308, 0, 3, 95, 101, 183, 193, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 323, 0, 3, 113, 123, 203, 218, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 344, 0, 3, 123, 133, 218, 233, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 365, 0, 3, 133, 143, 233, 248, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 386, 0, 3, 143, 153, 248, 263, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 407, 0, 3, 153, 163, 263, 278, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 428, 0, 3, 163, 173, 278, 293, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 449, 0, 3, 173, 183, 293, 308, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 470, 0, 3, 203, 218, 323, 344, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 498, 0, 3, 218, 233, 344, 365, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 526, 0, 3, 233, 248, 365, 386, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 554, 0, 3, 248, 263, 386, 407, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 582, 0, 3, 263, 278, 407, 428, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 610, 0, 3, 278, 293, 428, 449, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 638, 0, 3, 323, 344, 470, 498, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 674, 0, 3, 344, 365, 498, 526, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 710, 0, 3, 365, 386, 526, 554, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 746, 0, 3, 386, 407, 554, 582, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 782, 0, 3, 407, 428, 582, 610, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 818, 0, 3, 470, 498, 638, 674, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 863, 0, 3, 498, 526, 674, 710, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 908, 0, 3, 526, 554, 710, 746, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 953, 0, 3, 554, 582, 746, 782, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 998, 0, 3, 638, 674, 818, 863, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 1053, 0, 3, 674, 710, 863, 908,
                                                    ncols, p);

                compute_prim_ms_nuclear_potential_0(buffer, 1108, 0, 3, 710, 746, 908, 953,
                                                    ncols, p);

                compute_prim_ns_nuclear_potential_0(buffer, 1163, 0, 3, 818, 863, 998, 1053,
                                                    ncols, p);

                compute_prim_ns_nuclear_potential_0(buffer, 1229, 0, 3, 863, 908, 1053, 1108,
                                                    ncols, p);

                compute_prim_os_nuclear_potential_0(buffer, 1295, 0, 3, 998, 1053, 1163, 1229,
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

    simdtrf::compute_hrr_ip(buffer, coordinates, 1681, 1373, 1401, 1, nmax);

    simdtrf::compute_hrr_kp(buffer, coordinates, 1765, 1401, 1437, 1, nmax);

    simdtrf::compute_hrr_lp(buffer, coordinates, 1873, 1437, 1482, 1, nmax);

    simdtrf::compute_hrr_mp(buffer, coordinates, 2008, 1482, 1537, 1, nmax);

    simdtrf::compute_hrr_np(buffer, coordinates, 2173, 1537, 1603, 1, nmax);

    simdtrf::compute_hrr_id(buffer, coordinates, 2371, 1681, 1765, 1, nmax);

    simdtrf::compute_hrr_kd(buffer, coordinates, 2539, 1765, 1873, 1, nmax);

    simdtrf::compute_hrr_ld(buffer, coordinates, 2755, 1873, 2008, 1, nmax);

    simdtrf::compute_hrr_md(buffer, coordinates, 3025, 2008, 2173, 1, nmax);

    simdtrf::compute_hrr_if(buffer, coordinates, 3355, 2371, 2539, 1, nmax);

    simdtrf::compute_hrr_kf(buffer, coordinates, 3635, 2539, 2755, 1, nmax);

    simdtrf::compute_hrr_lf(buffer, coordinates, 3995, 2755, 3025, 1, nmax);

    simdtrf::compute_hrr_ig(buffer, coordinates, 4445, 3355, 3635, 1, nmax);

    simdtrf::compute_hrr_kg(buffer, coordinates, 4865, 3635, 3995, 1, nmax);

    simdtrf::compute_hrr_ih(buffer, coordinates, 5405, 4445, 4865, 1, nmax);

    simdtrf::transform_h_inner(buffer, 5993, 5405, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 5993, 11, nmax);

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
