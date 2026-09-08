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


#include "SimdElectronRepulsionRecHF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformHF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hf_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1121, nvalues);

    buffer.zero();

    const auto nmax = nvalues;

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = nvalues;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto pi = mathconst::pi_value();

            const auto fj = 2.0 * a_norms[i] * b_norms[j] * pi * pi * std::sqrt(pi)
                            / (a_exps[i] * b_exps[j] * std::sqrt(p));

            const auto alpha = a_exps[i];

            const auto beta = b_exps[j];

            const auto fb = a_exps[i] / p;

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pb(buffer, coordinates, 3, ncols, fb);

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 15, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 36, 0, 7, 8, 18, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 39, 0, 8, 9, 21, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 42, 0, 9, 10, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 45, 0, 10, 11, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 48, 0, 11, 12, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 12, 13, 33, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 54, 0, 15, 18, 39, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 60, 0, 18, 21, 42, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 66, 0, 21, 24, 45, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_3(buffer, 75, 0, 24, 27, 48, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 83, 0, 27, 30, 51, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 89, 0, 36, 39, 60, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_6(buffer, 95, 0, 39, 42, 66, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_7(buffer, 104, 0, 42, 45, 75, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_4(buffer, 116, 0, 45, 48, 83, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 127, 0, 54, 60, 95, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_2(buffer, 136, 0, 60, 66, 104, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_4(buffer, 145, 0, 66, 75, 116, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 170, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 173, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 176, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 179, 3, 12, 30, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 182, 3, 18, 39, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 185, 3, 21, 42, ncols, p);

            compute_prim_dp_electron_repulsion_8(buffer, 188, 0, 3, 24, 176, 45, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 198, 3, 27, 48, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 207, 3, 30, 51, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 210, 3, 36, 54, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 213, 3, 39, 60, ncols, p);

            compute_prim_fp_electron_repulsion_13(buffer, 216, 0, 3, 42, 188, 66, ncols, p);

            compute_prim_fp_electron_repulsion_9(buffer, 228, 0, 3, 45, 198, 75, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 240, 3, 48, 83, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 249, 3, 60, 95, ncols, p);

            compute_prim_gp_electron_repulsion_9(buffer, 261, 0, 3, 66, 228, 104, ncols, p);

            compute_prim_gp_electron_repulsion_10(buffer, 283, 0, 3, 75, 240, 116, ncols, p);

            compute_prim_hp_electron_repulsion_2(buffer, 300, 3, 89, 127, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 315, 3, 95, 136, ncols, p);

            compute_prim_hp_electron_repulsion_4(buffer, 330, 0, 3, 104, 283, 145, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 360, 3, 9, 10, 161, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 363, 3, 10, 11, 164, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 366, 3, 11, 12, 167, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 369, 0, 158, 360, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 372, 0, 161, 363, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 375, 0, 164, 366, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 378, 3, 170, 36, 39, 185, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_10(buffer, 381, 3, 173, 39, 42, 188, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 384, 3, 176, 42, 45, 198, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 393, 3, 179, 45, 48, 207, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_16(buffer, 396, 0, 3, 185, 381, 54, 60, 216, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_6(buffer, 405, 0, 3, 188, 384, 60, 66, 228, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_17(buffer, 428, 0, 3, 198, 393, 66, 75, 240, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_9(buffer, 443, 0, 3, 378, 381, 216, 405, 89, 95, 261, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_10(buffer, 491, 0, 3, 381, 384, 228, 428, 95, 104, 283, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_1(buffer, 524, 0, 3, 396, 405, 261, 491, 127, 136, 330, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 605, 3, 369, 182, 185, 381, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_18(buffer, 608, 3, 372, 185, 188, 384, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_19(buffer, 611, 3, 375, 188, 198, 393, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 614, 0, 3, 378, 605, 210, 213, 396, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_16(buffer, 623, 0, 3, 381, 608, 213, 216, 405, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_17(buffer, 632, 0, 3, 384, 611, 216, 228, 428, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_10(buffer, 647, 0, 3, 605, 608, 405, 632, 249, 261, 491, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 701, 0, 3, 614, 623, 443, 647, 300, 315, 524, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 911, 701, 210, ncols);
        }
    }

    simdtrf::transform_hf(values, nvalues, buffer, 911, nmax);
}

}  // namespace simdt2ceri
