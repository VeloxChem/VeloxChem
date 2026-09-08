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


#include "SimdElectronRepulsionRecHL.hpp"

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
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHL.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
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
#include "SimdTransformHL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hl_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(9133, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 59, 0, 7, 8, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 62, 0, 8, 9, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 65, 0, 9, 10, 32, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 68, 0, 10, 11, 35, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 71, 0, 11, 12, 38, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 74, 0, 12, 13, 41, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 77, 0, 13, 14, 44, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 80, 0, 14, 15, 47, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 83, 0, 15, 16, 50, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 86, 0, 16, 17, 53, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 89, 0, 17, 18, 56, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 92, 0, 20, 23, 59, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 98, 0, 23, 26, 62, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 104, 0, 26, 29, 65, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 110, 0, 29, 32, 68, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 116, 0, 32, 35, 71, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 122, 0, 35, 38, 74, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 128, 0, 38, 41, 77, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 134, 0, 41, 44, 80, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 140, 0, 44, 47, 83, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 146, 0, 47, 50, 86, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 152, 0, 50, 53, 89, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 158, 0, 59, 62, 104, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 167, 0, 62, 65, 110, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 176, 0, 65, 68, 116, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 185, 0, 68, 71, 122, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 194, 0, 71, 74, 128, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 203, 0, 74, 77, 134, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 212, 0, 77, 80, 140, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 221, 0, 80, 83, 146, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 230, 0, 83, 86, 152, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 239, 0, 92, 98, 158, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 248, 0, 98, 104, 167, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 257, 0, 104, 110, 176, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 266, 0, 110, 116, 185, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 275, 0, 116, 122, 194, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 284, 0, 122, 128, 203, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 293, 0, 128, 134, 212, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 302, 0, 134, 140, 221, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 311, 0, 140, 146, 230, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 320, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 323, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 326, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 329, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 332, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 335, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 338, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 341, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 344, 3, 18, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 347, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 350, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 353, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 356, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 359, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 362, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 365, 3, 15, 47, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 368, 3, 16, 50, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 371, 3, 17, 53, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 374, 3, 26, 62, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 377, 3, 29, 65, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 380, 3, 32, 68, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 389, 3, 35, 71, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 398, 3, 38, 74, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 407, 3, 41, 77, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 416, 3, 44, 80, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 425, 3, 47, 83, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 434, 3, 50, 86, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 443, 3, 53, 89, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 446, 3, 62, 104, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 449, 3, 65, 110, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 458, 3, 68, 116, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 467, 3, 71, 122, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 476, 3, 74, 128, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 485, 3, 77, 134, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 494, 3, 80, 140, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 503, 3, 83, 146, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 512, 3, 86, 152, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 521, 3, 104, 167, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 533, 3, 110, 176, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 545, 3, 116, 185, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 557, 3, 122, 194, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 569, 3, 128, 203, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 581, 3, 134, 212, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 593, 3, 140, 221, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 605, 3, 146, 230, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 617, 3, 167, 257, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 632, 3, 176, 266, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 647, 3, 185, 275, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 662, 3, 194, 284, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 677, 3, 203, 293, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 692, 3, 212, 302, ncols, p);

            compute_prim_hp_electron_repulsion_3(buffer, 707, 3, 221, 311, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 722, 3, 9, 10, 323, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 725, 3, 10, 11, 326, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 728, 3, 11, 12, 329, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 731, 3, 12, 13, 332, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 734, 3, 13, 14, 335, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 737, 3, 14, 15, 338, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 740, 3, 15, 16, 341, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 743, 3, 16, 17, 344, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 746, 0, 320, 722, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 749, 0, 323, 725, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 752, 0, 326, 728, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 755, 0, 329, 731, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 758, 0, 332, 734, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 761, 0, 335, 737, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 764, 0, 338, 740, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 767, 0, 341, 743, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 770, 3, 347, 59, 62, 377, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 773, 3, 350, 62, 65, 380, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 776, 3, 353, 65, 68, 389, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 785, 3, 356, 68, 71, 398, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 794, 3, 359, 71, 74, 407, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 803, 3, 362, 74, 77, 416, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 812, 3, 365, 77, 80, 425, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 821, 3, 368, 80, 83, 434, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 830, 3, 371, 83, 86, 443, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 833, 0, 3, 374, 770, 92, 98, 446, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 842, 0, 3, 377, 773, 98, 104, 449, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 851, 0, 3, 380, 776, 104, 110, 458, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 866, 0, 3, 389, 785, 110, 116, 467, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 881, 0, 3, 398, 794, 116, 122, 476, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 896, 0, 3, 407, 803, 122, 128, 485, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 911, 0, 3, 416, 812, 128, 134, 494, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 926, 0, 3, 425, 821, 134, 140, 503, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_19(buffer, 941, 0, 3, 434, 830, 140, 146, 512, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 956, 0, 3, 770, 773, 449, 851, 158, 167, 533, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_19(buffer, 980, 0, 3, 773, 776, 458, 866, 167, 176, 545, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1004, 0, 3, 776, 785, 467, 881, 176, 185, 557, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1028, 0, 3, 785, 794, 476, 896, 185, 194, 569, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1052, 0, 3, 794, 803, 485, 911, 194, 203, 581, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1076, 0, 3, 803, 812, 494, 926, 203, 212, 593, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_17(buffer, 1100, 0, 3, 812, 821, 503, 941, 212, 221, 605, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_2(buffer, 1124, 0, 3, 833, 842, 521, 956, 239, 248, 617, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_5(buffer, 1151, 0, 3, 842, 851, 533, 980, 248, 257, 632, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1178, 0, 3, 851, 866, 545, 1004, 257, 266, 647, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1205, 0, 3, 866, 881, 557, 1028, 266, 275, 662, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1232, 0, 3, 881, 896, 569, 1052, 275, 284, 677, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1259, 0, 3, 896, 911, 581, 1076, 284, 293, 692, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_8(buffer, 1286, 0, 3, 911, 926, 593, 1100, 293, 302, 707, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1313, 3, 320, 323, 725, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1316, 3, 323, 326, 728, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1319, 3, 326, 329, 731, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1322, 3, 329, 332, 734, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1325, 3, 332, 335, 737, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1328, 3, 335, 338, 740, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 1331, 3, 338, 341, 743, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1334, 0, 722, 1313, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1337, 0, 725, 1316, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1340, 0, 728, 1319, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1343, 0, 731, 1322, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1346, 0, 734, 1325, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1349, 0, 737, 1328, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 1352, 0, 740, 1331, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 1355, 3, 746, 374, 377, 773, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_14(buffer, 1358, 3, 749, 377, 380, 776, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1361, 3, 752, 380, 389, 785, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1379, 3, 755, 389, 398, 794, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1397, 3, 758, 398, 407, 803, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1415, 3, 761, 407, 416, 812, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 1433, 3, 764, 416, 425, 821, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_21(buffer, 1451, 3, 767, 425, 434, 830, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_21(buffer, 1454, 0, 3, 773, 1358, 446, 449, 851, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1463, 0, 3, 776, 1361, 449, 458, 866, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1490, 0, 3, 785, 1379, 458, 467, 881, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1517, 0, 3, 794, 1397, 467, 476, 896, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1544, 0, 3, 803, 1415, 476, 485, 911, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_14(buffer, 1571, 0, 3, 812, 1433, 485, 494, 926, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_20(buffer, 1598, 0, 3, 821, 1451, 494, 503, 941, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_16(buffer, 1622, 0, 3, 1355, 1358, 851, 1463, 521, 533, 980, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_21(buffer, 1664, 0, 3, 1358, 1361, 866, 1490, 533, 545, 1004, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_22(buffer, 1709, 0, 3, 1361, 1379, 881, 1517, 545, 557, 1028, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_22(buffer, 1754, 0, 3, 1379, 1397, 896, 1544, 557, 569, 1052, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_22(buffer, 1799, 0, 3, 1397, 1415, 911, 1571, 569, 581, 1076, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_19(buffer, 1844, 0, 3, 1415, 1433, 926, 1598, 581, 593, 1100, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_5(buffer, 1886, 0, 3, 1454, 1463, 980, 1664, 617, 632, 1178, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_8(buffer, 1940, 0, 3, 1463, 1490, 1004, 1709, 632, 647, 1205, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_8(buffer, 1994, 0, 3, 1490, 1517, 1028, 1754, 647, 662, 1232, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_8(buffer, 2048, 0, 3, 1517, 1544, 1052, 1799, 662, 677, 1259, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_9(buffer, 2102, 0, 3, 1544, 1571, 1076, 1844, 677, 692, 1286, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2156, 3, 722, 725, 1316, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2159, 3, 725, 728, 1319, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2162, 3, 728, 731, 1322, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2165, 3, 731, 734, 1325, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2168, 3, 734, 737, 1328, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 2171, 3, 737, 740, 1331, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2174, 0, 1313, 2156, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2177, 0, 1316, 2159, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2180, 0, 1319, 2162, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2183, 0, 1322, 2165, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2186, 0, 1325, 2168, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 2189, 0, 1328, 2171, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 2192, 3, 1334, 770, 773, 1358, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_14(buffer, 2195, 3, 1337, 773, 776, 1361, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 2198, 3, 1340, 776, 785, 1379, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 2225, 3, 1343, 785, 794, 1397, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 2252, 3, 1346, 794, 803, 1415, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 2279, 3, 1349, 803, 812, 1433, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_22(buffer, 2306, 3, 1352, 812, 821, 1451, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_16(buffer, 2309, 0, 3, 1355, 2192, 833, 842, 1454, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_22(buffer, 2318, 0, 3, 1358, 2195, 842, 851, 1463, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 2327, 0, 3, 1361, 2198, 851, 866, 1490, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 2366, 0, 3, 1379, 2225, 866, 881, 1517, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 2405, 0, 3, 1397, 2252, 881, 896, 1544, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_15(buffer, 2444, 0, 3, 1415, 2279, 896, 911, 1571, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_25(buffer, 2483, 0, 3, 1433, 2306, 911, 926, 1598, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_20(buffer, 2516, 0, 3, 2192, 2195, 1463, 2327, 956, 980, 1664, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_21(buffer, 2576, 0, 3, 2195, 2198, 1490, 2366, 980, 1004, 1709, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_22(buffer, 2642, 0, 3, 2198, 2225, 1517, 2405, 1004, 1028, 1754, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_22(buffer, 2708, 0, 3, 2225, 2252, 1544, 2444, 1028, 1052, 1799, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_19(buffer, 2774, 0, 3, 2252, 2279, 1571, 2483, 1052, 1076, 1844, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_2(buffer, 2834, 0, 3, 2309, 2318, 1622, 2516, 1124, 1151, 1886, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_5(buffer, 2915, 0, 3, 2318, 2327, 1664, 2576, 1151, 1178, 1940, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_8(buffer, 2996, 0, 3, 2327, 2366, 1709, 2642, 1178, 1205, 1994, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_8(buffer, 3077, 0, 3, 2366, 2405, 1754, 2708, 1205, 1232, 2048, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_9(buffer, 3158, 0, 3, 2405, 2444, 1799, 2774, 1232, 1259, 2102, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3239, 3, 1313, 1316, 2159, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3242, 3, 1316, 1319, 2162, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3245, 3, 1319, 1322, 2165, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3248, 3, 1322, 1325, 2168, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 3251, 3, 1325, 1328, 2171, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3254, 0, 2156, 3239, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3257, 0, 2159, 3242, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3260, 0, 2162, 3245, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3263, 0, 2165, 3248, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 3266, 0, 2168, 3251, ncols, p);

            compute_prim_dh_electron_repulsion_9(buffer, 3269, 3, 2174, 1355, 1358, 2195, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_15(buffer, 3272, 3, 2177, 1358, 1361, 2198, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_27(buffer, 3275, 3, 2180, 1361, 1379, 2225, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_26(buffer, 3314, 3, 2183, 1379, 1397, 2252, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_26(buffer, 3350, 3, 2186, 1397, 1415, 2279, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_23(buffer, 3386, 3, 2189, 1415, 1433, 2306, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_20(buffer, 3389, 0, 3, 2195, 3272, 1454, 1463, 2327, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_24(buffer, 3398, 0, 3, 2198, 3275, 1463, 1490, 2366, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_25(buffer, 3452, 0, 3, 2225, 3314, 1490, 1517, 2405, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_26(buffer, 3512, 0, 3, 2252, 3350, 1517, 1544, 2444, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_23(buffer, 3566, 0, 3, 2279, 3386, 1544, 1571, 2483, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_14(buffer, 3611, 0, 3, 3269, 3272, 2327, 3398, 1622, 1664, 2576, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_15(buffer, 3701, 0, 3, 3272, 3275, 2366, 3452, 1664, 1709, 2642, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_16(buffer, 3794, 0, 3, 3275, 3314, 2405, 3512, 1709, 1754, 2708, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_17(buffer, 3890, 0, 3, 3314, 3350, 2444, 3566, 1754, 1799, 2774, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_5(buffer, 3974, 0, 3, 3389, 3398, 2576, 3701, 1886, 1940, 2996, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_6(buffer, 4091, 0, 3, 3398, 3452, 2642, 3794, 1940, 1994, 3077, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_7(buffer, 4208, 0, 3, 3452, 3512, 2708, 3890, 1994, 2048, 3158, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 4325, 3, 2156, 2159, 3242, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 4328, 3, 2159, 2162, 3245, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 4331, 3, 2162, 2165, 3248, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 4334, 3, 2165, 2168, 3251, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 4337, 0, 3239, 4325, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 4340, 0, 3242, 4328, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 4343, 0, 3245, 4331, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 4346, 0, 3248, 4334, ncols, p);

            compute_prim_di_electron_repulsion_9(buffer, 4349, 3, 3254, 2192, 2195, 3272, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_15(buffer, 4352, 3, 3257, 2195, 2198, 3275, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_23(buffer, 4355, 0, 3, 3260, 4343, 2198, 2225, 3314, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_24(buffer, 4403, 3, 3263, 2225, 2252, 3350, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_22(buffer, 4448, 3, 3266, 2252, 2279, 3386, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_10(buffer, 4451, 0, 3, 3269, 4349, 2309, 2318, 3389, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_16(buffer, 4460, 0, 3, 3272, 4352, 2318, 2327, 3398, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_17(buffer, 4469, 0, 3, 3275, 4355, 2327, 2366, 3452, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_18(buffer, 4571, 0, 3, 3314, 4403, 2366, 2405, 3512, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_19(buffer, 4652, 0, 3, 3350, 4448, 2405, 2444, 3566, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_8(buffer, 4712, 0, 3, 4349, 4352, 3398, 4469, 2516, 2576, 3701, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_9(buffer, 4826, 0, 3, 4352, 4355, 3452, 4571, 2576, 2642, 3794, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_10(buffer, 4985, 0, 3, 4355, 4403, 3512, 4652, 2642, 2708, 3890, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_2(buffer, 5105, 0, 3, 4451, 4460, 3611, 4712, 2834, 2915, 3974, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_3(buffer, 5267, 0, 3, 4460, 4469, 3701, 4826, 2915, 2996, 4091, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_4(buffer, 5429, 0, 3, 4469, 4571, 3794, 4985, 2996, 3077, 4208, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 5615, 3, 3239, 3242, 4328, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 5618, 3, 3242, 3245, 4331, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 5621, 3, 3245, 3248, 4334, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_8(buffer, 5624, 0, 4325, 5615, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 5627, 0, 4328, 5618, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 5630, 0, 4331, 5621, ncols, p);

            compute_prim_dk_electron_repulsion_5(buffer, 5633, 3, 4337, 3269, 3272, 4352, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_13(buffer, 5636, 3, 4340, 3272, 3275, 4355, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_14(buffer, 5639, 3, 4343, 3275, 3314, 4403, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_15(buffer, 5684, 3, 4346, 3314, 3350, 4448, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_8(buffer, 5687, 0, 3, 4352, 5636, 3389, 3398, 4469, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_9(buffer, 5696, 0, 3, 4355, 5639, 3398, 3452, 4571, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_10(buffer, 5827, 0, 3, 4403, 5684, 3452, 3512, 4652, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_3(buffer, 5905, 0, 3, 5633, 5636, 4469, 5696, 3611, 3701, 4826, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_4(buffer, 6245, 0, 3, 5636, 5639, 4571, 5827, 3701, 3794, 4985, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_1(buffer, 6448, 0, 3, 5687, 5696, 4826, 6245, 3974, 4091, 5429, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_2(buffer, 6895, 3, 5624, 4349, 4352, 5636, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_5(buffer, 6898, 3, 5627, 4352, 4355, 5639, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_6(buffer, 6901, 3, 5630, 4355, 4403, 5684, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_2(buffer, 6904, 0, 3, 5633, 6895, 4451, 4460, 5687, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_3(buffer, 6913, 0, 3, 5636, 6898, 4460, 4469, 5696, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_4(buffer, 6922, 0, 3, 5639, 6901, 4469, 4571, 5827, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_1(buffer, 7000, 0, 3, 6895, 6898, 5696, 6922, 4712, 4826, 6245, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 7243, 0, 3, 6904, 6913, 5905, 7000, 5105, 5267, 6448, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 8188, 7243, 945, ncols);
        }
    }

    simdtrf::transform_hl(values, nvalues, buffer, 8188, nmax);
}

}  // namespace simdt2ceri
