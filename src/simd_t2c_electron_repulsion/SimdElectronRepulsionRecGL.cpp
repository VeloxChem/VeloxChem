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


#include "SimdElectronRepulsionRecGL.hpp"

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
#include "SimdElectronRepulsionVrrRecDH.hpp"
#include "SimdElectronRepulsionVrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDK.hpp"
#include "SimdElectronRepulsionVrrRecDL.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFL.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGL.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformGL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gl_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(5935, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 12, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 53, 0, 7, 8, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 56, 0, 8, 9, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 59, 0, 9, 10, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 62, 0, 10, 11, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 65, 0, 11, 12, 32, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 68, 0, 12, 13, 35, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 71, 0, 13, 14, 38, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 74, 0, 14, 15, 41, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 77, 0, 15, 16, 44, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 80, 0, 16, 17, 47, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 83, 0, 17, 18, 50, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 86, 0, 20, 23, 59, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 92, 0, 23, 26, 62, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 98, 0, 26, 29, 65, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 104, 0, 29, 32, 68, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 110, 0, 32, 35, 71, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 116, 0, 35, 38, 74, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 122, 0, 38, 41, 77, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 128, 0, 41, 44, 80, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 134, 0, 44, 47, 83, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 140, 0, 53, 56, 86, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 146, 0, 56, 59, 92, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 152, 0, 59, 62, 98, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 158, 0, 62, 65, 104, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 164, 0, 65, 68, 110, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 170, 0, 68, 71, 116, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 176, 0, 71, 74, 122, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 182, 0, 74, 77, 128, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 188, 0, 77, 80, 134, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 194, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 197, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 200, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 203, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 206, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 209, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 212, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 215, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 218, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 221, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 224, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 227, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 230, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 233, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 236, 3, 14, 38, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 239, 3, 15, 41, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 242, 3, 16, 44, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 245, 3, 17, 47, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 248, 3, 23, 59, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 251, 3, 26, 62, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 260, 3, 29, 65, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 269, 3, 32, 68, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 278, 3, 35, 71, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 287, 3, 38, 74, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 296, 3, 41, 77, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 305, 3, 44, 80, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 314, 3, 47, 83, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 323, 3, 59, 92, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 332, 3, 62, 98, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 341, 3, 65, 104, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 350, 3, 68, 110, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 359, 3, 71, 116, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 368, 3, 74, 122, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 377, 3, 77, 128, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 386, 3, 80, 134, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 395, 3, 92, 152, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 407, 3, 98, 158, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 419, 3, 104, 164, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 431, 3, 110, 170, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 443, 3, 116, 176, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 455, 3, 122, 182, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 467, 3, 128, 188, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 479, 3, 9, 10, 197, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 482, 3, 10, 11, 200, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 485, 3, 11, 12, 203, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 488, 3, 12, 13, 206, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 491, 3, 13, 14, 209, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 494, 3, 14, 15, 212, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 497, 3, 15, 16, 215, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 500, 3, 16, 17, 218, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 503, 0, 197, 482, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 506, 0, 200, 485, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 509, 0, 203, 488, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 512, 0, 206, 491, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 515, 0, 209, 494, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 518, 0, 212, 497, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 521, 0, 215, 500, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 524, 3, 221, 53, 56, 248, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 527, 3, 224, 56, 59, 251, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 530, 3, 227, 59, 62, 260, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 539, 3, 230, 62, 65, 269, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 548, 3, 233, 65, 68, 278, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 557, 3, 236, 68, 71, 287, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 566, 3, 239, 71, 74, 296, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 575, 3, 242, 74, 77, 305, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 584, 3, 245, 77, 80, 314, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 593, 0, 3, 251, 530, 86, 92, 332, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 608, 0, 3, 260, 539, 92, 98, 341, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 623, 0, 3, 269, 548, 98, 104, 350, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 638, 0, 3, 278, 557, 104, 110, 359, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 653, 0, 3, 287, 566, 110, 116, 368, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 668, 0, 3, 296, 575, 116, 122, 377, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 683, 0, 3, 305, 584, 122, 128, 386, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_2(buffer, 698, 0, 3, 524, 527, 323, 593, 140, 146, 395, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_5(buffer, 716, 0, 3, 527, 530, 332, 608, 146, 152, 407, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 734, 0, 3, 530, 539, 341, 623, 152, 158, 419, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 752, 0, 3, 539, 548, 350, 638, 158, 164, 431, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 770, 0, 3, 548, 557, 359, 653, 164, 170, 443, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 788, 0, 3, 557, 566, 368, 668, 170, 176, 455, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 806, 0, 3, 566, 575, 377, 683, 176, 182, 467, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 824, 3, 194, 197, 482, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 827, 3, 197, 200, 485, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 830, 3, 200, 203, 488, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 836, 3, 203, 206, 491, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 842, 3, 206, 209, 494, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 848, 3, 209, 212, 497, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 854, 3, 212, 215, 500, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 860, 0, 479, 824, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 863, 0, 482, 827, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 866, 0, 485, 830, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 869, 0, 488, 836, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 872, 0, 491, 842, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 875, 0, 494, 848, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 878, 0, 497, 854, ncols, p);

            compute_prim_df_electron_repulsion_14(buffer, 881, 3, 503, 248, 251, 530, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 884, 3, 506, 251, 260, 539, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 902, 3, 509, 260, 269, 548, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 920, 3, 512, 269, 278, 557, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 938, 3, 515, 278, 287, 566, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 956, 3, 518, 287, 296, 575, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 974, 3, 521, 296, 305, 584, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_10(buffer, 992, 0, 3, 530, 884, 323, 332, 608, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1016, 0, 3, 539, 902, 332, 341, 623, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1043, 0, 3, 548, 920, 341, 350, 638, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1070, 0, 3, 557, 938, 350, 359, 653, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1097, 0, 3, 566, 956, 359, 368, 668, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_13(buffer, 1124, 0, 3, 575, 974, 368, 377, 683, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_5(buffer, 1148, 0, 3, 881, 884, 608, 1016, 395, 407, 734, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_8(buffer, 1184, 0, 3, 884, 902, 623, 1043, 407, 419, 752, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_8(buffer, 1220, 0, 3, 902, 920, 638, 1070, 419, 431, 770, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_8(buffer, 1256, 0, 3, 920, 938, 653, 1097, 431, 443, 788, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_9(buffer, 1292, 0, 3, 938, 956, 668, 1124, 443, 455, 806, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 1328, 3, 479, 482, 827, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_8(buffer, 1331, 3, 482, 485, 830, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 1334, 3, 485, 488, 836, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 1343, 3, 488, 491, 842, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 1352, 3, 491, 494, 848, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 1361, 3, 494, 497, 854, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 1370, 0, 827, 1331, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 1373, 0, 830, 1334, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 1376, 0, 836, 1343, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 1379, 0, 842, 1352, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 1382, 0, 848, 1361, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 1385, 3, 860, 524, 527, 881, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_14(buffer, 1388, 3, 863, 527, 530, 884, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1391, 3, 866, 530, 539, 902, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1418, 3, 869, 539, 548, 920, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1445, 3, 872, 548, 557, 938, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1472, 3, 875, 557, 566, 956, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1499, 3, 878, 566, 575, 974, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_14(buffer, 1526, 0, 3, 884, 1391, 593, 608, 1016, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 1559, 0, 3, 902, 1418, 608, 623, 1043, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 1598, 0, 3, 920, 1445, 623, 638, 1070, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 1637, 0, 3, 938, 1472, 638, 653, 1097, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_13(buffer, 1676, 0, 3, 956, 1499, 653, 668, 1124, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_2(buffer, 1709, 0, 3, 1385, 1388, 992, 1526, 698, 716, 1148, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_5(buffer, 1763, 0, 3, 1388, 1391, 1016, 1559, 716, 734, 1184, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_8(buffer, 1817, 0, 3, 1391, 1418, 1043, 1598, 734, 752, 1220, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_8(buffer, 1871, 0, 3, 1418, 1445, 1070, 1637, 752, 770, 1256, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_9(buffer, 1925, 0, 3, 1445, 1472, 1097, 1676, 770, 788, 1292, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 1979, 3, 824, 827, 1331, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_18(buffer, 1982, 3, 827, 830, 1334, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_19(buffer, 1985, 3, 830, 836, 1343, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_17(buffer, 2000, 3, 836, 842, 1352, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_17(buffer, 2012, 3, 842, 848, 1361, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 2024, 0, 1328, 1979, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 2027, 0, 1331, 1982, ncols, p);

            compute_prim_ph_electron_repulsion_20(buffer, 2030, 0, 1334, 1985, ncols, p);

            compute_prim_ph_electron_repulsion_19(buffer, 2033, 0, 1343, 2000, ncols, p);

            compute_prim_ph_electron_repulsion_19(buffer, 2036, 0, 1352, 2012, ncols, p);

            compute_prim_dh_electron_repulsion_15(buffer, 2039, 3, 1370, 881, 884, 1391, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 2042, 3, 1373, 884, 902, 1418, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_18(buffer, 2081, 3, 1376, 902, 920, 1445, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 2123, 3, 1379, 920, 938, 1472, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 2162, 3, 1382, 938, 956, 1499, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_10(buffer, 2201, 0, 3, 1391, 2042, 992, 1016, 1559, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_11(buffer, 2252, 0, 3, 1418, 2081, 1016, 1043, 1598, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_12(buffer, 2306, 0, 3, 1445, 2123, 1043, 1070, 1637, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_13(buffer, 2363, 0, 3, 1472, 2162, 1070, 1097, 1676, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_5(buffer, 2408, 0, 3, 2039, 2042, 1559, 2252, 1148, 1184, 1817, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_6(buffer, 2486, 0, 3, 2042, 2081, 1598, 2306, 1184, 1220, 1871, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_7(buffer, 2564, 0, 3, 2081, 2123, 1637, 2363, 1220, 1256, 1925, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 2642, 3, 1328, 1331, 1982, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_14(buffer, 2645, 3, 1331, 1334, 1985, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_15(buffer, 2648, 3, 1334, 1343, 2000, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_16(buffer, 2663, 3, 1343, 1352, 2012, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 2678, 0, 1982, 2645, ncols, p);

            compute_prim_pi_electron_repulsion_16(buffer, 2681, 0, 3, 1985, 2648, 2033, ncols, p);

            compute_prim_pi_electron_repulsion_17(buffer, 2699, 0, 2000, 2663, ncols, p);

            compute_prim_di_electron_repulsion_9(buffer, 2702, 3, 2024, 1385, 1388, 2039, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_15(buffer, 2705, 3, 2027, 1388, 1391, 2042, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_16(buffer, 2708, 0, 3, 2030, 2681, 1391, 1418, 2081, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_17(buffer, 2774, 0, 3, 2033, 2699, 1418, 1445, 2123, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_2(buffer, 2834, 3, 2036, 1445, 1472, 2162, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_7(buffer, 2888, 0, 3, 2042, 2708, 1526, 1559, 2252, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_8(buffer, 2948, 0, 3, 2081, 2774, 1559, 1598, 2306, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_9(buffer, 3050, 0, 3, 2123, 2834, 1598, 1637, 2363, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_2(buffer, 3116, 0, 3, 2702, 2705, 2201, 2888, 1709, 1763, 2408, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_3(buffer, 3224, 0, 3, 2705, 2708, 2252, 2948, 1763, 1817, 2486, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_4(buffer, 3332, 0, 3, 2708, 2774, 2306, 3050, 1817, 1871, 2564, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 3458, 3, 1979, 1982, 2645, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_10(buffer, 3461, 3, 1982, 1985, 2648, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_11(buffer, 3464, 3, 1985, 2000, 2663, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_8(buffer, 3479, 0, 2642, 3458, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 3482, 0, 2645, 3461, ncols, p);

            compute_prim_pk_electron_repulsion_9(buffer, 3485, 0, 2648, 3464, ncols, p);

            compute_prim_dk_electron_repulsion_8(buffer, 3500, 3, 2678, 2039, 2042, 2708, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_9(buffer, 3503, 0, 3, 2681, 3485, 2042, 2081, 2774, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_10(buffer, 3599, 3, 2699, 2081, 2123, 2834, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_3(buffer, 3671, 0, 3, 2708, 3503, 2201, 2252, 2948, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_4(buffer, 3898, 0, 3, 2774, 3599, 2252, 2306, 3050, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_1(buffer, 4027, 0, 3, 3500, 3503, 2948, 3898, 2408, 2486, 3332, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_2(buffer, 4342, 3, 3479, 2702, 2705, 3500, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_3(buffer, 4345, 3, 3482, 2705, 2708, 3503, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_4(buffer, 4348, 3, 3485, 2708, 2774, 3599, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_1(buffer, 4420, 0, 3, 3503, 4348, 2888, 2948, 3898, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 4585, 0, 3, 4342, 4345, 3671, 4420, 3116, 3224, 4027, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 5260, 4585, 675, ncols);
        }
    }

    simdtrf::transform_gl(values, nvalues, buffer, 5260, nmax);
}

}  // namespace simdt2ceri
