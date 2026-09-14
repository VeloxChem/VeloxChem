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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3696, 1008, 230, dimensions);

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

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 6, 3, 10, ncols,
                                                          fz, mu, p);

                compute_prim_ps_nuclear_potential_0(buffer, 18, 0, 3, 7, 8, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 21, 0, 3, 8, 9, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 24, 0, 3, 9, 10, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 27, 0, 3, 10, 11, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 30, 0, 3, 11, 12, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 33, 0, 3, 12, 13, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 36, 0, 3, 13, 14, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 39, 0, 3, 14, 15, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 42, 0, 3, 15, 16, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 45, 0, 3, 16, 17, ncols);

                compute_prim_ds_nuclear_potential_0(buffer, 48, 0, 3, 7, 8, 18, 21, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 54, 0, 3, 8, 9, 21, 24, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 60, 0, 3, 9, 10, 24, 27, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 66, 0, 3, 10, 11, 27, 30, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 72, 0, 3, 11, 12, 30, 33, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 78, 0, 3, 12, 13, 33, 36, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 84, 0, 3, 13, 14, 36, 39, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 90, 0, 3, 14, 15, 39, 42, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 96, 0, 3, 15, 16, 42, 45, ncols, p);

                compute_prim_fs_nuclear_potential_0(buffer, 102, 0, 3, 18, 21, 48, 54, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 112, 0, 3, 21, 24, 54, 60, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 122, 0, 3, 24, 27, 60, 66, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 132, 0, 3, 27, 30, 66, 72, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 142, 0, 3, 30, 33, 72, 78, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 152, 0, 3, 33, 36, 78, 84, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 162, 0, 3, 36, 39, 84, 90, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 172, 0, 3, 39, 42, 90, 96, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 182, 0, 3, 48, 54, 102, 112, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 197, 0, 3, 54, 60, 112, 122, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 212, 0, 3, 60, 66, 122, 132, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 227, 0, 3, 66, 72, 132, 142, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 242, 0, 3, 72, 78, 142, 152, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 257, 0, 3, 78, 84, 152, 162, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 272, 0, 3, 84, 90, 162, 172, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 287, 0, 3, 102, 112, 182, 197, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 308, 0, 3, 112, 122, 197, 212, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 329, 0, 3, 122, 132, 212, 227, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 350, 0, 3, 132, 142, 227, 242, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 371, 0, 3, 142, 152, 242, 257, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 392, 0, 3, 152, 162, 257, 272, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 413, 0, 3, 182, 197, 287, 308, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 441, 0, 3, 197, 212, 308, 329, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 469, 0, 3, 212, 227, 329, 350, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 497, 0, 3, 227, 242, 350, 371, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 525, 0, 3, 242, 257, 371, 392, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 553, 0, 3, 287, 308, 413, 441, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 589, 0, 3, 308, 329, 441, 469, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 625, 0, 3, 329, 350, 469, 497, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 661, 0, 3, 350, 371, 497, 525, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 697, 0, 3, 413, 441, 553, 589, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 742, 0, 3, 441, 469, 589, 625, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 787, 0, 3, 469, 497, 625, 661, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 832, 0, 3, 553, 589, 697, 742, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 887, 0, 3, 589, 625, 742, 787, ncols,
                                                    p);

                compute_prim_ns_nuclear_potential_0(buffer, 942, 0, 3, 697, 742, 832, 887, ncols,
                                                    p);

                simdfunc::contract_primitives(buffer, 1008, 413, 28, ncols);

                simdfunc::contract_primitives(buffer, 1036, 553, 36, ncols);

                simdfunc::contract_primitives(buffer, 1072, 697, 45, ncols);

                simdfunc::contract_primitives(buffer, 1117, 832, 55, ncols);

                simdfunc::contract_primitives(buffer, 1172, 942, 66, ncols);
            }
        }
    }

    simdtrf::compute_hrr_ip(buffer, coordinates, 1238, 1008, 1036, 1, nmax);

    simdtrf::compute_hrr_kp(buffer, coordinates, 1322, 1036, 1072, 1, nmax);

    simdtrf::compute_hrr_lp(buffer, coordinates, 1430, 1072, 1117, 1, nmax);

    simdtrf::compute_hrr_mp(buffer, coordinates, 1565, 1117, 1172, 1, nmax);

    simdtrf::compute_hrr_id(buffer, coordinates, 1730, 1238, 1322, 1, nmax);

    simdtrf::compute_hrr_kd(buffer, coordinates, 1898, 1322, 1430, 1, nmax);

    simdtrf::compute_hrr_ld(buffer, coordinates, 2114, 1430, 1565, 1, nmax);

    simdtrf::compute_hrr_if(buffer, coordinates, 2384, 1730, 1898, 1, nmax);

    simdtrf::compute_hrr_kf(buffer, coordinates, 2664, 1898, 2114, 1, nmax);

    simdtrf::compute_hrr_ig(buffer, coordinates, 3024, 2384, 2664, 1, nmax);

    simdtrf::transform_g_inner(buffer, 3444, 3024, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 3444, 9, nmax);

    for (size_t m = 0; m < 117; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
