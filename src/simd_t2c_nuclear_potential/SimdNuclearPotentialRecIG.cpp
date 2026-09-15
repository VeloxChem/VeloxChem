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


#include "SimdNuclearPotentialRecIG.hpp"

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
#include "SimdNuclearPotentialVrrRecPS.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_ig_nuclear_potential(double                    *values,
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
            false, std::string("compute_ig_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_ig_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 117 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3697, 1009, 230, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 117 * nvalues, 0.0);

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

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 10, ncols,
                                                          fz, 6, p);

                compute_prim_ps_nuclear_potential_0(buffer, 19, 0, 3, 8, 9, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 22, 0, 3, 9, 10, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 25, 0, 3, 10, 11, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 28, 0, 3, 11, 12, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 31, 0, 3, 12, 13, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 34, 0, 3, 13, 14, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 37, 0, 3, 14, 15, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 40, 0, 3, 15, 16, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 43, 0, 3, 16, 17, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 46, 0, 3, 17, 18, ncols);

                compute_prim_ds_nuclear_potential_0(buffer, 49, 0, 3, 8, 9, 19, 22, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 55, 0, 3, 9, 10, 22, 25, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 61, 0, 3, 10, 11, 25, 28, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 67, 0, 3, 11, 12, 28, 31, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 73, 0, 3, 12, 13, 31, 34, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 79, 0, 3, 13, 14, 34, 37, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 85, 0, 3, 14, 15, 37, 40, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 91, 0, 3, 15, 16, 40, 43, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 97, 0, 3, 16, 17, 43, 46, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 103, 0, 3, 19, 22, 49, 55, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 113, 0, 3, 22, 25, 55, 61, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 123, 0, 3, 25, 28, 61, 67, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 133, 0, 3, 28, 31, 67, 73, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 143, 0, 3, 31, 34, 73, 79, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 153, 0, 3, 34, 37, 79, 85, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 163, 0, 3, 37, 40, 85, 91, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 173, 0, 3, 40, 43, 91, 97, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 183, 0, 3, 49, 55, 103, 113, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 198, 0, 3, 55, 61, 113, 123, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 213, 0, 3, 61, 67, 123, 133, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 228, 0, 3, 67, 73, 133, 143, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 243, 0, 3, 73, 79, 143, 153, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 258, 0, 3, 79, 85, 153, 163, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 273, 0, 3, 85, 91, 163, 173, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 288, 0, 3, 103, 113, 183, 198, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 309, 0, 3, 113, 123, 198, 213, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 330, 0, 3, 123, 133, 213, 228, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 351, 0, 3, 133, 143, 228, 243, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 372, 0, 3, 143, 153, 243, 258, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 393, 0, 3, 153, 163, 258, 273, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 414, 0, 3, 183, 198, 288, 309, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 442, 0, 3, 198, 213, 309, 330, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 470, 0, 3, 213, 228, 330, 351, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 498, 0, 3, 228, 243, 351, 372, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 526, 0, 3, 243, 258, 372, 393, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 554, 0, 3, 288, 309, 414, 442, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 590, 0, 3, 309, 330, 442, 470, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 626, 0, 3, 330, 351, 470, 498, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 662, 0, 3, 351, 372, 498, 526, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 698, 0, 3, 414, 442, 554, 590, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 743, 0, 3, 442, 470, 590, 626, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 788, 0, 3, 470, 498, 626, 662, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 833, 0, 3, 554, 590, 698, 743, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 888, 0, 3, 590, 626, 743, 788, ncols,
                                                    p);

                compute_prim_ns_nuclear_potential_0(buffer, 943, 0, 3, 698, 743, 833, 888, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 1009, 414, 28, ncols);

                simdfunc::contract_primitives(buffer, 1037, 554, 36, ncols);

                simdfunc::contract_primitives(buffer, 1073, 698, 45, ncols);

                simdfunc::contract_primitives(buffer, 1118, 833, 55, ncols);

                simdfunc::contract_primitives(buffer, 1173, 943, 66, ncols);
            }
        }
    }

    simdtrf::compute_hrr_ip(buffer, coordinates, 1239, 1009, 1037, 1, nmax);

    simdtrf::compute_hrr_kp(buffer, coordinates, 1323, 1037, 1073, 1, nmax);

    simdtrf::compute_hrr_lp(buffer, coordinates, 1431, 1073, 1118, 1, nmax);

    simdtrf::compute_hrr_mp(buffer, coordinates, 1566, 1118, 1173, 1, nmax);

    simdtrf::compute_hrr_id(buffer, coordinates, 1731, 1239, 1323, 1, nmax);

    simdtrf::compute_hrr_kd(buffer, coordinates, 1899, 1323, 1431, 1, nmax);

    simdtrf::compute_hrr_ld(buffer, coordinates, 2115, 1431, 1566, 1, nmax);

    simdtrf::compute_hrr_if(buffer, coordinates, 2385, 1731, 1899, 1, nmax);

    simdtrf::compute_hrr_kf(buffer, coordinates, 2665, 1899, 2115, 1, nmax);

    simdtrf::compute_hrr_ig(buffer, coordinates, 3025, 2385, 2665, 1, nmax);

    simdtrf::transform_g_inner(buffer, 3445, 3025, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 3445, 9, nmax);

    for (size_t m = 0; m < 117; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
