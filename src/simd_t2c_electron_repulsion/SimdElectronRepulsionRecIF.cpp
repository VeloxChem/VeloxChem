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


#include "SimdElectronRepulsionRecIF.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformIF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_if_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_if_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1627, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 40, 0, 7, 8, 19, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 43, 0, 8, 9, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 46, 0, 9, 10, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 49, 0, 10, 11, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 11, 12, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 12, 13, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 58, 0, 13, 14, 37, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 61, 0, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 67, 0, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_3(buffer, 73, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_3(buffer, 81, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 89, 0, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 95, 0, 31, 34, 58, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 101, 0, 40, 43, 67, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_9(buffer, 110, 0, 43, 46, 73, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_7(buffer, 119, 0, 46, 49, 81, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_4(buffer, 131, 0, 49, 52, 89, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 142, 0, 52, 55, 95, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 151, 0, 61, 67, 110, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_11(buffer, 160, 0, 67, 73, 119, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_12(buffer, 172, 0, 73, 81, 131, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_13(buffer, 188, 0, 81, 89, 142, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 203, 0, 101, 110, 160, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_4(buffer, 215, 0, 110, 119, 172, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_5(buffer, 227, 0, 119, 131, 188, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 245, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 248, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 251, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 254, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 257, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 260, 3, 14, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 263, 3, 9, 22, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 266, 3, 10, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 269, 3, 11, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 272, 3, 12, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 275, 3, 13, 34, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 278, 3, 16, 40, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 281, 3, 19, 43, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 284, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 287, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 290, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 293, 3, 31, 55, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 296, 3, 34, 58, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 299, 3, 43, 67, ncols, p);

            compute_prim_fp_electron_repulsion_14(buffer, 302, 3, 46, 73, ncols, p);

            compute_prim_fp_electron_repulsion_15(buffer, 305, 0, 3, 49, 290, 81, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 317, 3, 52, 89, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 326, 3, 55, 95, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 329, 3, 61, 101, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 332, 3, 67, 110, ncols, p);

            compute_prim_gp_electron_repulsion_15(buffer, 335, 0, 3, 73, 305, 119, ncols, p);

            compute_prim_gp_electron_repulsion_10(buffer, 353, 0, 3, 81, 317, 131, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 370, 3, 89, 142, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 379, 3, 110, 160, ncols, p);

            compute_prim_hp_electron_repulsion_9(buffer, 394, 0, 3, 119, 353, 172, ncols, p);

            compute_prim_hp_electron_repulsion_10(buffer, 424, 0, 3, 131, 370, 188, ncols, p);

            compute_prim_ip_electron_repulsion_2(buffer, 446, 3, 151, 203, ncols, p);

            compute_prim_ip_electron_repulsion_3(buffer, 464, 3, 160, 215, ncols, p);

            compute_prim_ip_electron_repulsion_4(buffer, 482, 0, 3, 172, 424, 227, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 521, 3, 8, 9, 248, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 524, 3, 9, 10, 251, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 527, 3, 10, 11, 254, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 530, 3, 11, 12, 257, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 533, 3, 12, 13, 260, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 536, 0, 245, 521, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 539, 0, 248, 524, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 542, 0, 251, 527, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 545, 0, 254, 530, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 548, 0, 257, 533, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 551, 3, 263, 40, 43, 284, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 554, 3, 266, 43, 46, 287, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 557, 3, 269, 46, 49, 290, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 560, 3, 272, 49, 52, 293, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 563, 3, 275, 52, 55, 296, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 566, 0, 3, 284, 554, 61, 67, 302, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_22(buffer, 575, 0, 3, 287, 557, 67, 73, 305, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_23(buffer, 584, 0, 3, 290, 560, 73, 81, 317, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_24(buffer, 599, 0, 3, 293, 563, 81, 89, 326, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_23(buffer, 608, 0, 3, 551, 554, 302, 575, 101, 110, 335, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_24(buffer, 623, 0, 3, 554, 557, 305, 584, 110, 119, 353, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_25(buffer, 656, 0, 3, 557, 560, 317, 599, 119, 131, 370, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_10(buffer, 677, 0, 3, 566, 575, 335, 623, 151, 160, 394, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_11(buffer, 747, 0, 3, 575, 584, 353, 656, 160, 172, 424, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_1(buffer, 793, 0, 3, 608, 623, 394, 747, 203, 215, 482, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 902, 3, 536, 278, 281, 551, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 905, 3, 539, 281, 284, 554, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 908, 3, 542, 284, 287, 557, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 911, 3, 545, 287, 290, 560, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 914, 3, 548, 290, 293, 563, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 917, 0, 3, 554, 908, 299, 302, 575, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_25(buffer, 926, 0, 3, 557, 911, 302, 305, 584, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_26(buffer, 935, 0, 3, 560, 914, 305, 317, 599, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_23(buffer, 944, 0, 3, 902, 905, 566, 917, 329, 332, 608, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_24(buffer, 959, 0, 3, 905, 908, 575, 926, 332, 335, 623, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_25(buffer, 974, 0, 3, 908, 911, 584, 935, 335, 353, 656, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_10(buffer, 995, 0, 3, 917, 926, 623, 974, 379, 394, 747, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 1067, 0, 3, 944, 959, 677, 995, 446, 464, 793, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1347, 1067, 280, ncols);
        }
    }

    simdtrf::transform_if(values, nvalues, buffer, 1347, nmax);
}

}  // namespace simdt2ceri
