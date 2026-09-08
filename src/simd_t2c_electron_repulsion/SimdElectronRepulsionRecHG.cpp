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


#include "SimdElectronRepulsionRecHG.hpp"

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
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformHG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hg_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1890, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 43, 0, 7, 8, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 46, 0, 8, 9, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 49, 0, 9, 10, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 10, 11, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 11, 12, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 58, 0, 12, 13, 37, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 61, 0, 13, 14, 40, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 64, 0, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 70, 0, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 76, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 82, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 88, 0, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 94, 0, 31, 34, 58, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 100, 0, 34, 37, 61, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 106, 0, 43, 46, 76, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 115, 0, 46, 49, 82, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 124, 0, 49, 52, 88, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 133, 0, 52, 55, 94, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 142, 0, 55, 58, 100, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 151, 0, 64, 70, 106, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 160, 0, 70, 76, 115, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 169, 0, 76, 82, 124, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 178, 0, 82, 88, 133, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 187, 0, 88, 94, 142, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 196, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 199, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 202, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 205, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 208, 3, 14, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 211, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 214, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 217, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 220, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 223, 3, 13, 37, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 226, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 229, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 232, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 241, 3, 31, 55, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 250, 3, 34, 58, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 259, 3, 37, 61, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 262, 3, 46, 76, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 265, 3, 49, 82, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 274, 3, 52, 88, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 283, 3, 55, 94, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 292, 3, 58, 100, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 301, 3, 76, 115, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 313, 3, 82, 124, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 325, 3, 88, 133, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 337, 3, 94, 142, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 349, 3, 115, 169, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 364, 3, 124, 178, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 379, 3, 133, 187, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 394, 3, 9, 10, 199, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 397, 3, 10, 11, 202, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 400, 3, 11, 12, 205, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 403, 3, 12, 13, 208, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 406, 0, 196, 394, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 409, 0, 199, 397, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 412, 0, 202, 400, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 415, 0, 205, 403, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 418, 3, 211, 43, 46, 229, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 421, 3, 214, 46, 49, 232, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_7(buffer, 424, 0, 3, 217, 412, 49, 52, 241, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 436, 3, 220, 52, 55, 250, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 445, 3, 223, 55, 58, 259, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 448, 0, 3, 226, 418, 64, 70, 262, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 457, 0, 3, 229, 421, 70, 76, 265, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_8(buffer, 466, 0, 3, 232, 424, 76, 82, 274, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_9(buffer, 490, 0, 3, 241, 436, 82, 88, 283, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 511, 0, 3, 250, 445, 88, 94, 292, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_11(buffer, 526, 0, 3, 418, 421, 265, 466, 106, 115, 313, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_12(buffer, 550, 0, 3, 421, 424, 274, 490, 115, 124, 325, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_13(buffer, 586, 0, 3, 424, 436, 283, 511, 124, 133, 337, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_2(buffer, 616, 0, 3, 448, 457, 301, 526, 151, 160, 349, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_3(buffer, 643, 0, 3, 457, 466, 313, 550, 160, 169, 364, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_4(buffer, 670, 0, 3, 466, 490, 325, 586, 169, 178, 379, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 721, 3, 196, 199, 397, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 724, 3, 199, 202, 400, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 727, 3, 202, 205, 403, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 730, 0, 394, 721, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 733, 0, 397, 724, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 736, 0, 400, 727, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 739, 3, 406, 226, 229, 421, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_12(buffer, 742, 3, 409, 229, 232, 424, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_20(buffer, 745, 3, 412, 232, 241, 436, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_21(buffer, 754, 3, 415, 241, 250, 445, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_18(buffer, 757, 0, 3, 421, 742, 262, 265, 466, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_19(buffer, 766, 0, 3, 424, 745, 265, 274, 490, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_20(buffer, 807, 0, 3, 436, 754, 274, 283, 511, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_11(buffer, 831, 0, 3, 739, 742, 466, 766, 301, 313, 550, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_12(buffer, 928, 0, 3, 742, 745, 490, 807, 313, 325, 586, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_1(buffer, 987, 0, 3, 757, 766, 550, 928, 349, 364, 670, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 1128, 3, 730, 418, 421, 742, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_18(buffer, 1131, 3, 733, 421, 424, 745, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_19(buffer, 1134, 3, 736, 424, 436, 754, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_16(buffer, 1137, 0, 3, 739, 1128, 448, 457, 757, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_17(buffer, 1146, 0, 3, 742, 1131, 457, 466, 766, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_18(buffer, 1155, 0, 3, 745, 1134, 466, 490, 807, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_10(buffer, 1179, 0, 3, 1128, 1131, 766, 1155, 526, 550, 928, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 1260, 0, 3, 1137, 1146, 831, 1179, 616, 643, 987, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1575, 1260, 315, ncols);
        }
    }

    simdtrf::transform_hg(values, nvalues, buffer, 1575, nmax);
}

}  // namespace simdt2ceri
