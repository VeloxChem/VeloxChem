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


#include "SimdNuclearPotentialRecHG.hpp"

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
#include "SimdNuclearPotentialVrrRecPS.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHG.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_hg_nuclear_potential(double                    *values,
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
            false, std::string("compute_hg_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_hg_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 99 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2801, 722, 185, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 99 * nvalues, 0.0);

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

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 6, 3, 9, ncols,
                                                          fz, mu, p);

                compute_prim_ps_nuclear_potential_0(buffer, 17, 0, 3, 7, 8, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 20, 0, 3, 8, 9, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 23, 0, 3, 9, 10, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 26, 0, 3, 10, 11, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 29, 0, 3, 11, 12, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 32, 0, 3, 12, 13, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 35, 0, 3, 13, 14, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 38, 0, 3, 14, 15, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 41, 0, 3, 15, 16, ncols);

                compute_prim_ds_nuclear_potential_0(buffer, 44, 0, 3, 7, 8, 17, 20, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 50, 0, 3, 8, 9, 20, 23, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 56, 0, 3, 9, 10, 23, 26, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 62, 0, 3, 10, 11, 26, 29, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 68, 0, 3, 11, 12, 29, 32, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 74, 0, 3, 12, 13, 32, 35, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 80, 0, 3, 13, 14, 35, 38, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 86, 0, 3, 14, 15, 38, 41, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 92, 0, 3, 17, 20, 44, 50, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 102, 0, 3, 20, 23, 50, 56, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 112, 0, 3, 23, 26, 56, 62, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 122, 0, 3, 26, 29, 62, 68, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 132, 0, 3, 29, 32, 68, 74, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 142, 0, 3, 32, 35, 74, 80, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 152, 0, 3, 35, 38, 80, 86, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 162, 0, 3, 44, 50, 92, 102, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 177, 0, 3, 50, 56, 102, 112, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 192, 0, 3, 56, 62, 112, 122, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 207, 0, 3, 62, 68, 122, 132, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 222, 0, 3, 68, 74, 132, 142, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 237, 0, 3, 74, 80, 142, 152, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 252, 0, 3, 92, 102, 162, 177, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 273, 0, 3, 102, 112, 177, 192, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 294, 0, 3, 112, 122, 192, 207, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 315, 0, 3, 122, 132, 207, 222, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 336, 0, 3, 132, 142, 222, 237, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 357, 0, 3, 162, 177, 252, 273, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 385, 0, 3, 177, 192, 273, 294, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 413, 0, 3, 192, 207, 294, 315, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 441, 0, 3, 207, 222, 315, 336, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 469, 0, 3, 252, 273, 357, 385, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 505, 0, 3, 273, 294, 385, 413, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 541, 0, 3, 294, 315, 413, 441, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 577, 0, 3, 357, 385, 469, 505, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 622, 0, 3, 385, 413, 505, 541, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 667, 0, 3, 469, 505, 577, 622, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 722, 252, 21, ncols);

                simdfunc::contract_primitives(buffer, 743, 357, 28, ncols);

                simdfunc::contract_primitives(buffer, 771, 469, 36, ncols);

                simdfunc::contract_primitives(buffer, 807, 577, 45, ncols);

                simdfunc::contract_primitives(buffer, 852, 667, 55, ncols);
            }
        }
    }

    simdtrf::compute_hrr_hp(buffer, coordinates, 907, 722, 743, 1, nmax);

    simdtrf::compute_hrr_ip(buffer, coordinates, 970, 743, 771, 1, nmax);

    simdtrf::compute_hrr_kp(buffer, coordinates, 1054, 771, 807, 1, nmax);

    simdtrf::compute_hrr_lp(buffer, coordinates, 1162, 807, 852, 1, nmax);

    simdtrf::compute_hrr_hd(buffer, coordinates, 1297, 907, 970, 1, nmax);

    simdtrf::compute_hrr_id(buffer, coordinates, 1423, 970, 1054, 1, nmax);

    simdtrf::compute_hrr_kd(buffer, coordinates, 1591, 1054, 1162, 1, nmax);

    simdtrf::compute_hrr_hf(buffer, coordinates, 1807, 1297, 1423, 1, nmax);

    simdtrf::compute_hrr_if(buffer, coordinates, 2017, 1423, 1591, 1, nmax);

    simdtrf::compute_hrr_hg(buffer, coordinates, 2297, 1807, 2017, 1, nmax);

    simdtrf::transform_g_inner(buffer, 2612, 2297, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 2612, 9, nmax);

    for (size_t m = 0; m < 99; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
