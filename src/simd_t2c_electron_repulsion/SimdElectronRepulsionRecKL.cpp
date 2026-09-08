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


#include "SimdElectronRepulsionRecKL.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIL.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKL.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
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
#include "SimdTransformKL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kl_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(18849, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 67, 0, 7, 8, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 70, 0, 8, 9, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 73, 0, 9, 10, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 76, 0, 10, 11, 37, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 79, 0, 11, 12, 40, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 82, 0, 12, 13, 43, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 85, 0, 13, 14, 46, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 88, 0, 14, 15, 49, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 91, 0, 15, 16, 52, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 94, 0, 16, 17, 55, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 97, 0, 17, 18, 58, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 100, 0, 18, 19, 61, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 103, 0, 19, 20, 64, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 106, 0, 22, 25, 67, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 112, 0, 25, 28, 70, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 118, 0, 28, 31, 73, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 124, 0, 31, 34, 76, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 130, 0, 34, 37, 79, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 136, 0, 37, 40, 82, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 142, 0, 40, 43, 85, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 148, 0, 43, 46, 88, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 154, 0, 46, 49, 91, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 160, 0, 49, 52, 94, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 166, 0, 52, 55, 97, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 172, 0, 55, 58, 100, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 178, 0, 58, 61, 103, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 184, 0, 67, 70, 118, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 193, 0, 70, 73, 124, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 202, 0, 73, 76, 130, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 211, 0, 76, 79, 136, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 220, 0, 79, 82, 142, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 229, 0, 82, 85, 148, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 238, 0, 85, 88, 154, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 247, 0, 88, 91, 160, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 256, 0, 91, 94, 166, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 265, 0, 94, 97, 172, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 274, 0, 97, 100, 178, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 283, 0, 106, 112, 184, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 295, 0, 112, 118, 193, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 307, 0, 118, 124, 202, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 319, 0, 124, 130, 211, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 331, 0, 130, 136, 220, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 343, 0, 136, 142, 229, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 355, 0, 142, 148, 238, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 367, 0, 148, 154, 247, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 379, 0, 154, 160, 256, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 391, 0, 160, 166, 265, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 403, 0, 166, 172, 274, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 415, 0, 184, 193, 307, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 430, 0, 193, 202, 319, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 445, 0, 202, 211, 331, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 460, 0, 211, 220, 343, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 475, 0, 220, 229, 355, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 490, 0, 229, 238, 367, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 505, 0, 238, 247, 379, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 520, 0, 247, 256, 391, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 535, 0, 256, 265, 403, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 550, 0, 283, 295, 415, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 565, 0, 295, 307, 430, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 580, 0, 307, 319, 445, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 595, 0, 319, 331, 460, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 610, 0, 331, 343, 475, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 625, 0, 343, 355, 490, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 640, 0, 355, 367, 505, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 655, 0, 367, 379, 520, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 670, 0, 379, 391, 535, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 685, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 688, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 691, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 694, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 697, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 700, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 703, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 706, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 709, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 712, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 715, 3, 20, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 718, 3, 9, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 721, 3, 10, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 724, 3, 11, 37, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 727, 3, 12, 40, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 730, 3, 13, 43, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 733, 3, 14, 46, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 736, 3, 15, 49, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 739, 3, 16, 52, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 742, 3, 17, 55, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 745, 3, 18, 58, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 748, 3, 19, 61, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 751, 3, 28, 70, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 754, 3, 31, 73, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 757, 3, 34, 76, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 760, 3, 37, 79, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 763, 3, 40, 82, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 766, 3, 43, 85, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 769, 3, 46, 88, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 772, 3, 49, 91, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 775, 3, 52, 94, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 778, 3, 55, 97, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 781, 3, 58, 100, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 784, 3, 61, 103, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 787, 3, 70, 118, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 790, 3, 73, 124, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 793, 3, 76, 130, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 796, 3, 79, 136, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 805, 3, 82, 142, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 814, 3, 85, 148, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 823, 3, 88, 154, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 832, 3, 91, 160, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 841, 3, 94, 166, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 850, 3, 97, 172, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 853, 3, 100, 178, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 856, 3, 118, 193, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 859, 3, 124, 202, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 862, 3, 130, 211, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 874, 3, 136, 220, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 886, 3, 142, 229, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 898, 3, 148, 238, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 910, 3, 154, 247, ncols, p);

            compute_prim_gp_electron_repulsion_8(buffer, 922, 3, 160, 256, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 934, 3, 166, 265, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 943, 3, 172, 274, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 946, 3, 193, 307, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 949, 3, 202, 319, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 964, 3, 211, 331, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 979, 3, 220, 343, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 994, 3, 229, 355, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1009, 3, 238, 367, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1024, 3, 247, 379, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1039, 3, 256, 391, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 1054, 3, 265, 403, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1063, 3, 307, 430, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1081, 3, 319, 445, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1099, 3, 331, 460, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1117, 3, 343, 475, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1135, 3, 355, 490, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1153, 3, 367, 505, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1171, 3, 379, 520, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1189, 3, 391, 535, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1207, 3, 430, 580, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1228, 3, 445, 595, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1249, 3, 460, 610, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1270, 3, 475, 625, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1291, 3, 490, 640, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1312, 3, 505, 655, ncols, p);

            compute_prim_kp_electron_repulsion_3(buffer, 1333, 3, 520, 670, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1354, 3, 9, 10, 688, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1357, 3, 10, 11, 691, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1360, 3, 11, 12, 694, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1363, 3, 12, 13, 697, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1366, 3, 13, 14, 700, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1369, 3, 14, 15, 703, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1372, 3, 15, 16, 706, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1375, 3, 16, 17, 709, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1378, 3, 17, 18, 712, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1381, 3, 18, 19, 715, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1384, 0, 685, 1354, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1387, 0, 688, 1357, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1390, 0, 691, 1360, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1393, 0, 694, 1363, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1396, 0, 697, 1366, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1399, 0, 700, 1369, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1402, 0, 703, 1372, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1405, 0, 706, 1375, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1408, 0, 709, 1378, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1411, 0, 712, 1381, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1414, 3, 718, 67, 70, 754, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1417, 3, 721, 70, 73, 757, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1420, 3, 724, 73, 76, 760, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1423, 3, 727, 76, 79, 763, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1426, 3, 730, 79, 82, 766, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1429, 3, 733, 82, 85, 769, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1432, 3, 736, 85, 88, 772, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1435, 3, 739, 88, 91, 775, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1438, 3, 742, 91, 94, 778, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1441, 3, 745, 94, 97, 781, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1444, 3, 748, 97, 100, 784, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1447, 0, 3, 751, 1414, 106, 112, 787, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1456, 0, 3, 754, 1417, 112, 118, 790, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1465, 0, 3, 757, 1420, 118, 124, 793, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_18(buffer, 1474, 0, 3, 760, 1423, 124, 130, 796, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1483, 0, 3, 763, 1426, 130, 136, 805, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1498, 0, 3, 766, 1429, 136, 142, 814, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1513, 0, 3, 769, 1432, 142, 148, 823, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1528, 0, 3, 772, 1435, 148, 154, 832, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_26(buffer, 1543, 0, 3, 775, 1438, 154, 160, 841, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1558, 0, 3, 778, 1441, 160, 166, 850, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1567, 0, 3, 781, 1444, 166, 172, 853, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1576, 0, 3, 1414, 1417, 790, 1465, 184, 193, 859, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_26(buffer, 1591, 0, 3, 1417, 1420, 793, 1474, 193, 202, 862, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1606, 0, 3, 1420, 1423, 796, 1483, 202, 211, 874, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1630, 0, 3, 1423, 1426, 805, 1498, 211, 220, 886, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1654, 0, 3, 1426, 1429, 814, 1513, 220, 229, 898, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1678, 0, 3, 1429, 1432, 823, 1528, 229, 238, 910, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_18(buffer, 1702, 0, 3, 1432, 1435, 832, 1543, 238, 247, 922, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_37(buffer, 1726, 0, 3, 1435, 1438, 841, 1558, 247, 256, 934, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1747, 0, 3, 1438, 1441, 850, 1567, 256, 265, 943, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 1762, 0, 3, 1447, 1456, 856, 1576, 283, 295, 946, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 1783, 0, 3, 1456, 1465, 859, 1591, 295, 307, 949, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_19(buffer, 1804, 0, 3, 1465, 1474, 862, 1606, 307, 319, 964, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_20(buffer, 1837, 0, 3, 1474, 1483, 874, 1630, 319, 331, 979, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1870, 0, 3, 1483, 1498, 886, 1654, 331, 343, 994, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1903, 0, 3, 1498, 1513, 898, 1678, 343, 355, 1009, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_21(buffer, 1936, 0, 3, 1513, 1528, 910, 1702, 355, 367, 1024, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_18(buffer, 1969, 0, 3, 1528, 1543, 922, 1726, 367, 379, 1039, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_34(buffer, 2002, 0, 3, 1543, 1558, 934, 1747, 379, 391, 1054, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_19(buffer, 2029, 0, 3, 1576, 1591, 949, 1804, 415, 430, 1081, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_20(buffer, 2071, 0, 3, 1591, 1606, 964, 1837, 430, 445, 1099, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 2113, 0, 3, 1606, 1630, 979, 1870, 445, 460, 1117, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 2155, 0, 3, 1630, 1654, 994, 1903, 460, 475, 1135, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 2197, 0, 3, 1654, 1678, 1009, 1936, 475, 490, 1153, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_21(buffer, 2239, 0, 3, 1678, 1702, 1024, 1969, 490, 505, 1171, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_18(buffer, 2281, 0, 3, 1702, 1726, 1039, 2002, 505, 520, 1189, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_2(buffer, 2323, 0, 3, 1762, 1783, 1063, 2029, 550, 565, 1207, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_5(buffer, 2368, 0, 3, 1783, 1804, 1081, 2071, 565, 580, 1228, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2413, 0, 3, 1804, 1837, 1099, 2113, 580, 595, 1249, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2458, 0, 3, 1837, 1870, 1117, 2155, 595, 610, 1270, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2503, 0, 3, 1870, 1903, 1135, 2197, 610, 625, 1291, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2548, 0, 3, 1903, 1936, 1153, 2239, 625, 640, 1312, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_8(buffer, 2593, 0, 3, 1936, 1969, 1171, 2281, 640, 655, 1333, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2638, 3, 685, 688, 1357, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2641, 3, 688, 691, 1360, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2644, 3, 691, 694, 1363, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2647, 3, 694, 697, 1366, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2650, 3, 697, 700, 1369, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2653, 3, 700, 703, 1372, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2656, 3, 703, 706, 1375, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2659, 3, 706, 709, 1378, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 2662, 3, 709, 712, 1381, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2665, 0, 1354, 2638, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2668, 0, 1357, 2641, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2671, 0, 1360, 2644, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2674, 0, 1363, 2647, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2677, 0, 1366, 2650, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2680, 0, 1369, 2653, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2683, 0, 1372, 2656, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2686, 0, 1375, 2659, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 2689, 0, 1378, 2662, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 2692, 3, 1384, 751, 754, 1417, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2695, 3, 1387, 754, 757, 1420, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2698, 3, 1390, 757, 760, 1423, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2704, 3, 1393, 760, 763, 1426, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2710, 3, 1396, 763, 766, 1429, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2716, 3, 1399, 766, 769, 1432, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2722, 3, 1402, 769, 772, 1435, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2728, 3, 1405, 772, 775, 1438, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 2734, 3, 1408, 775, 778, 1441, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 2740, 3, 1411, 778, 781, 1444, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 2743, 0, 3, 1417, 2695, 787, 790, 1465, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 2752, 0, 3, 1420, 2698, 790, 793, 1474, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_42(buffer, 2767, 0, 3, 1423, 2704, 793, 796, 1483, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2782, 0, 3, 1426, 2710, 796, 805, 1498, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2803, 0, 3, 1429, 2716, 805, 814, 1513, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2824, 0, 3, 1432, 2722, 814, 823, 1528, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_43(buffer, 2845, 0, 3, 1435, 2728, 823, 832, 1543, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_44(buffer, 2866, 0, 3, 1438, 2734, 832, 841, 1558, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_40(buffer, 2881, 0, 3, 1441, 2740, 841, 850, 1567, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_50(buffer, 2890, 0, 3, 2692, 2695, 1465, 2752, 856, 859, 1591, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_46(buffer, 2911, 0, 3, 2695, 2698, 1474, 2767, 859, 862, 1606, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_52(buffer, 2932, 0, 3, 2698, 2704, 1483, 2782, 862, 874, 1630, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_52(buffer, 2971, 0, 3, 2704, 2710, 1498, 2803, 874, 886, 1654, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_52(buffer, 3010, 0, 3, 2710, 2716, 1513, 2824, 886, 898, 1678, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_52(buffer, 3049, 0, 3, 2716, 2722, 1528, 2845, 898, 910, 1702, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_48(buffer, 3088, 0, 3, 2722, 2728, 1543, 2866, 910, 922, 1726, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_49(buffer, 3124, 0, 3, 2728, 2734, 1558, 2881, 922, 934, 1747, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_30(buffer, 3142, 0, 3, 2743, 2752, 1591, 2911, 946, 949, 1804, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_38(buffer, 3169, 0, 3, 2752, 2767, 1606, 2932, 949, 964, 1837, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_39(buffer, 3232, 0, 3, 2767, 2782, 1630, 2971, 964, 979, 1870, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_40(buffer, 3295, 0, 3, 2782, 2803, 1654, 3010, 979, 994, 1903, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_40(buffer, 3358, 0, 3, 2803, 2824, 1678, 3049, 994, 1009, 1936, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_36(buffer, 3421, 0, 3, 2824, 2845, 1702, 3088, 1009, 1024, 1969, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_37(buffer, 3484, 0, 3, 2845, 2866, 1726, 3124, 1024, 1039, 2002, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_16(buffer, 3526, 0, 3, 2890, 2911, 1804, 3169, 1063, 1081, 2071, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_21(buffer, 3604, 0, 3, 2911, 2932, 1837, 3232, 1081, 1099, 2113, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_22(buffer, 3685, 0, 3, 2932, 2971, 1870, 3295, 1099, 1117, 2155, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_22(buffer, 3766, 0, 3, 2971, 3010, 1903, 3358, 1117, 1135, 2197, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_22(buffer, 3847, 0, 3, 3010, 3049, 1936, 3421, 1135, 1153, 2239, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_19(buffer, 3928, 0, 3, 3049, 3088, 1969, 3484, 1153, 1171, 2281, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_5(buffer, 4006, 0, 3, 3142, 3169, 2071, 3604, 1207, 1228, 2413, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_8(buffer, 4096, 0, 3, 3169, 3232, 2113, 3685, 1228, 1249, 2458, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_8(buffer, 4186, 0, 3, 3232, 3295, 2155, 3766, 1249, 1270, 2503, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_8(buffer, 4276, 0, 3, 3295, 3358, 2197, 3847, 1270, 1291, 2548, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_9(buffer, 4366, 0, 3, 3358, 3421, 2239, 3928, 1291, 1312, 2593, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4456, 3, 1354, 1357, 2641, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4459, 3, 1357, 1360, 2644, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4462, 3, 1360, 1363, 2647, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4465, 3, 1363, 1366, 2650, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4468, 3, 1366, 1369, 2653, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4471, 3, 1369, 1372, 2656, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4474, 3, 1372, 1375, 2659, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 4477, 3, 1375, 1378, 2662, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4480, 0, 2638, 4456, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4483, 0, 2641, 4459, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4486, 0, 2644, 4462, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4489, 0, 2647, 4465, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4492, 0, 2650, 4468, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4495, 0, 2653, 4471, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4498, 0, 2656, 4474, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 4501, 0, 2659, 4477, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 4504, 3, 2665, 1414, 1417, 2695, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_27(buffer, 4507, 3, 2668, 1417, 1420, 2698, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4510, 3, 2671, 1420, 1423, 2704, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4519, 3, 2674, 1423, 1426, 2710, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4528, 3, 2677, 1426, 1429, 2716, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4537, 3, 2680, 1429, 1432, 2722, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4546, 3, 2683, 1432, 1435, 2728, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 4555, 3, 2686, 1435, 1438, 2734, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 4564, 3, 2689, 1438, 1441, 2740, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_16(buffer, 4567, 0, 3, 2692, 4504, 1447, 1456, 2743, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_44(buffer, 4576, 0, 3, 2695, 4507, 1456, 1465, 2752, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 4585, 0, 3, 2698, 4510, 1465, 1474, 2767, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_50(buffer, 4606, 0, 3, 2704, 4519, 1474, 1483, 2782, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_51(buffer, 4627, 0, 3, 2710, 4528, 1483, 1498, 2803, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_51(buffer, 4654, 0, 3, 2716, 4537, 1498, 1513, 2824, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_51(buffer, 4681, 0, 3, 2722, 4546, 1513, 1528, 2845, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_52(buffer, 4708, 0, 3, 2728, 4555, 1528, 1543, 2866, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_48(buffer, 4729, 0, 3, 2734, 4564, 1543, 1558, 2881, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_46(buffer, 4738, 0, 3, 4504, 4507, 2752, 4585, 1576, 1591, 2911, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_52(buffer, 4759, 0, 3, 4507, 4510, 2767, 4606, 1591, 1606, 2932, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_54(buffer, 4786, 0, 3, 4510, 4519, 2782, 4627, 1606, 1630, 2971, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_54(buffer, 4840, 0, 3, 4519, 4528, 2803, 4654, 1630, 1654, 3010, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_54(buffer, 4894, 0, 3, 4528, 4537, 2824, 4681, 1654, 1678, 3049, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_49(buffer, 4948, 0, 3, 4537, 4546, 2845, 4708, 1678, 1702, 3088, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_50(buffer, 4999, 0, 3, 4546, 4555, 2866, 4729, 1702, 1726, 3124, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_29(buffer, 5020, 0, 3, 4567, 4576, 2890, 4738, 1762, 1783, 3142, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_30(buffer, 5054, 0, 3, 4576, 4585, 2911, 4759, 1783, 1804, 3169, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_38(buffer, 5088, 0, 3, 4585, 4606, 2932, 4786, 1804, 1837, 3232, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_39(buffer, 5182, 0, 3, 4606, 4627, 2971, 4840, 1837, 1870, 3295, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_40(buffer, 5276, 0, 3, 4627, 4654, 3010, 4894, 1870, 1903, 3358, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_36(buffer, 5370, 0, 3, 4654, 4681, 3049, 4948, 1903, 1936, 3421, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_37(buffer, 5464, 0, 3, 4681, 4708, 3088, 4999, 1936, 1969, 3484, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_20(buffer, 5522, 0, 3, 4738, 4759, 3169, 5088, 2029, 2071, 3604, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_21(buffer, 5639, 0, 3, 4759, 4786, 3232, 5182, 2071, 2113, 3685, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_22(buffer, 5762, 0, 3, 4786, 4840, 3295, 5276, 2113, 2155, 3766, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_22(buffer, 5885, 0, 3, 4840, 4894, 3358, 5370, 2155, 2197, 3847, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_19(buffer, 6008, 0, 3, 4894, 4948, 3421, 5464, 2197, 2239, 3928, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_2(buffer, 6125, 0, 3, 5020, 5054, 3526, 5522, 2323, 2368, 4006, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_5(buffer, 6263, 0, 3, 5054, 5088, 3604, 5639, 2368, 2413, 4096, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_8(buffer, 6401, 0, 3, 5088, 5182, 3685, 5762, 2413, 2458, 4186, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_8(buffer, 6539, 0, 3, 5182, 5276, 3766, 5885, 2458, 2503, 4276, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_9(buffer, 6677, 0, 3, 5276, 5370, 3847, 6008, 2503, 2548, 4366, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6815, 3, 2638, 2641, 4459, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6818, 3, 2641, 2644, 4462, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6821, 3, 2644, 2647, 4465, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6824, 3, 2647, 2650, 4468, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6827, 3, 2650, 2653, 4471, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6830, 3, 2653, 2656, 4474, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 6833, 3, 2656, 2659, 4477, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6836, 0, 4456, 6815, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6839, 0, 4459, 6818, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6842, 0, 4462, 6821, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6845, 0, 4465, 6824, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6848, 0, 4468, 6827, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6851, 0, 4471, 6830, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 6854, 0, 4474, 6833, ncols, p);

            compute_prim_dh_electron_repulsion_9(buffer, 6857, 3, 4480, 2692, 2695, 4507, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_32(buffer, 6860, 3, 4483, 2695, 2698, 4510, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 6863, 3, 4486, 2698, 2704, 4519, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 6875, 3, 4489, 2704, 2710, 4528, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 6887, 3, 4492, 2710, 2716, 4537, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 6899, 3, 4495, 2716, 2722, 4546, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 6911, 3, 4498, 2722, 2728, 4555, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_40(buffer, 6923, 3, 4501, 2728, 2734, 4564, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_46(buffer, 6926, 0, 3, 4507, 6860, 2743, 2752, 4585, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 6935, 0, 3, 4510, 6863, 2752, 2767, 4606, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_52(buffer, 6962, 0, 3, 4519, 6875, 2767, 2782, 4627, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_53(buffer, 6989, 0, 3, 4528, 6887, 2782, 2803, 4654, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_53(buffer, 7022, 0, 3, 4537, 6899, 2803, 2824, 4681, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_54(buffer, 7055, 0, 3, 4546, 6911, 2824, 2845, 4708, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_50(buffer, 7082, 0, 3, 4555, 6923, 2845, 2866, 4729, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_44(buffer, 7091, 0, 3, 6857, 6860, 4585, 6935, 2890, 2911, 4759, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_45(buffer, 7124, 0, 3, 6860, 6863, 4606, 6962, 2911, 2932, 4786, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_46(buffer, 7157, 0, 3, 6863, 6875, 4627, 6989, 2932, 2971, 4840, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_47(buffer, 7230, 0, 3, 6875, 6887, 4654, 7022, 2971, 3010, 4894, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_42(buffer, 7300, 0, 3, 6887, 6899, 4681, 7055, 3010, 3049, 4948, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_43(buffer, 7366, 0, 3, 6899, 6911, 4708, 7082, 3049, 3088, 4999, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_25(buffer, 7390, 0, 3, 6926, 6935, 4759, 7124, 3142, 3169, 5088, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_29(buffer, 7431, 0, 3, 6935, 6962, 4786, 7157, 3169, 3232, 5182, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_30(buffer, 7561, 0, 3, 6962, 6989, 4840, 7230, 3232, 3295, 5276, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_31(buffer, 7698, 0, 3, 6989, 7022, 4894, 7300, 3295, 3358, 5370, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_32(buffer, 7828, 0, 3, 7022, 7055, 4948, 7366, 3358, 3421, 5464, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_14(buffer, 7905, 0, 3, 7091, 7124, 5088, 7431, 3526, 3604, 5639, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_15(buffer, 8082, 0, 3, 7124, 7157, 5182, 7561, 3604, 3685, 5762, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_16(buffer, 8262, 0, 3, 7157, 7230, 5276, 7698, 3685, 3766, 5885, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_17(buffer, 8445, 0, 3, 7230, 7300, 5370, 7828, 3766, 3847, 6008, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_5(buffer, 8613, 0, 3, 7390, 7431, 5639, 8082, 4006, 4096, 6401, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_6(buffer, 8817, 0, 3, 7431, 7561, 5762, 8262, 4096, 4186, 6539, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_7(buffer, 9021, 0, 3, 7561, 7698, 5885, 8445, 4186, 4276, 6677, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 9225, 3, 4456, 4459, 6818, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 9228, 3, 4459, 4462, 6821, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 9231, 3, 4462, 4465, 6824, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 9234, 3, 4465, 4468, 6827, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 9237, 3, 4468, 4471, 6830, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 9240, 3, 4471, 4474, 6833, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 9243, 0, 6815, 9225, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 9246, 0, 6818, 9228, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 9249, 0, 6821, 9231, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 9252, 0, 6824, 9234, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 9255, 0, 6827, 9237, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 9258, 0, 6830, 9240, ncols, p);

            compute_prim_di_electron_repulsion_9(buffer, 9261, 3, 6836, 4504, 4507, 6860, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_29(buffer, 9264, 3, 6839, 4507, 4510, 6863, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 9267, 3, 6842, 4510, 4519, 6875, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 9282, 3, 6845, 4519, 4528, 6887, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 9297, 3, 6848, 4528, 4537, 6899, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 9312, 3, 6851, 4537, 4546, 6911, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_37(buffer, 9327, 3, 6854, 4546, 4555, 6923, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_10(buffer, 9330, 0, 3, 6857, 9261, 4567, 4576, 6926, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_35(buffer, 9339, 0, 3, 6860, 9264, 4576, 4585, 6935, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_40(buffer, 9348, 0, 3, 6863, 9267, 4585, 4606, 6962, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_41(buffer, 9381, 0, 3, 6875, 9282, 4606, 4627, 6989, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_42(buffer, 9414, 0, 3, 6887, 9297, 4627, 4654, 7022, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_43(buffer, 9453, 0, 3, 6899, 9312, 4654, 4681, 7055, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_39(buffer, 9486, 0, 3, 6911, 9327, 4681, 4708, 7082, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_29(buffer, 9495, 0, 3, 9261, 9264, 6935, 9348, 4738, 4759, 7124, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_30(buffer, 9522, 0, 3, 9264, 9267, 6962, 9381, 4759, 4786, 7157, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_31(buffer, 9561, 0, 3, 9267, 9282, 6989, 9414, 4786, 4840, 7230, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_32(buffer, 9652, 0, 3, 9282, 9297, 7022, 9453, 4840, 4894, 7300, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_33(buffer, 9733, 0, 3, 9297, 9312, 7055, 9486, 4894, 4948, 7366, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_17(buffer, 9760, 0, 3, 9330, 9339, 7091, 9495, 5020, 5054, 7390, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_18(buffer, 9808, 0, 3, 9339, 9348, 7124, 9522, 5054, 5088, 7431, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_19(buffer, 9856, 0, 3, 9348, 9381, 7157, 9561, 5088, 5182, 7561, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_20(buffer, 10072, 0, 3, 9381, 9414, 7230, 9652, 5182, 5276, 7698, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_21(buffer, 10255, 0, 3, 9414, 9453, 7300, 9733, 5276, 5370, 7828, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_8(buffer, 10357, 0, 3, 9495, 9522, 7431, 9856, 5522, 5639, 8082, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_9(buffer, 10588, 0, 3, 9522, 9561, 7561, 10072, 5639, 5762, 8262, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_10(buffer, 10888, 0, 3, 9561, 9652, 7698, 10255, 5762, 5885, 8445, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_2(buffer, 11134, 0, 3, 9760, 9808, 7905, 10357, 6125, 6263, 8613, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_3(buffer, 11422, 0, 3, 9808, 9856, 8082, 10588, 6263, 6401, 8817, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_4(buffer, 11710, 0, 3, 9856, 10072, 8262, 10888, 6401, 6539, 9021, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 12043, 3, 6815, 6818, 9228, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 12046, 3, 6818, 6821, 9231, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 12049, 3, 6821, 6824, 9234, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 12052, 3, 6824, 6827, 9237, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 12055, 3, 6827, 6830, 9240, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_8(buffer, 12058, 0, 9225, 12043, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 12061, 0, 9228, 12046, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 12064, 0, 9231, 12049, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 12067, 0, 9234, 12052, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 12070, 0, 9237, 12055, ncols, p);

            compute_prim_dk_electron_repulsion_5(buffer, 12073, 3, 9243, 6857, 6860, 9264, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_20(buffer, 12076, 3, 9246, 6860, 6863, 9267, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 12079, 3, 9249, 6863, 6875, 9282, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 12094, 3, 9252, 6875, 6887, 9297, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 12109, 3, 9255, 6887, 6899, 9312, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_24(buffer, 12124, 3, 9258, 6899, 6911, 9327, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_22(buffer, 12127, 0, 3, 9264, 12076, 6926, 6935, 9348, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_23(buffer, 12136, 0, 3, 9267, 12079, 6935, 6962, 9381, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_24(buffer, 12172, 0, 3, 9282, 12094, 6962, 6989, 9414, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_25(buffer, 12208, 0, 3, 9297, 12109, 6989, 7022, 9453, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_26(buffer, 12244, 0, 3, 9312, 12124, 7022, 7055, 9486, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_14(buffer, 12253, 0, 3, 12073, 12076, 9348, 12136, 7091, 7124, 9522, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_15(buffer, 12298, 0, 3, 12076, 12079, 9381, 12172, 7124, 7157, 9561, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_16(buffer, 12343, 0, 3, 12079, 12094, 9414, 12208, 7157, 7230, 9652, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_17(buffer, 12430, 0, 3, 12094, 12109, 9453, 12244, 7230, 7300, 9733, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_8(buffer, 12460, 0, 3, 12127, 12136, 9522, 12298, 7390, 7431, 9856, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_9(buffer, 12526, 0, 3, 12136, 12172, 9561, 12343, 7431, 7561, 10072, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_10(buffer, 12778, 0, 3, 12172, 12208, 9652, 12430, 7561, 7698, 10255, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_3(buffer, 12913, 0, 3, 12253, 12298, 9856, 12526, 7905, 8082, 10588, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_4(buffer, 13551, 0, 3, 12298, 12343, 10072, 12778, 8082, 8262, 10888, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_1(buffer, 13947, 0, 3, 12460, 12526, 10588, 13551, 8613, 8817, 11710, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_2(buffer, 14730, 3, 12058, 9261, 9264, 12076, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_7(buffer, 14733, 3, 12061, 9264, 9267, 12079, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_8(buffer, 14736, 3, 12064, 9267, 9282, 12094, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_8(buffer, 14739, 3, 12067, 9282, 9297, 12109, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_11(buffer, 14742, 3, 12070, 9297, 9312, 12124, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_2(buffer, 14745, 0, 3, 12073, 14730, 9330, 9339, 12127, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_8(buffer, 14754, 0, 3, 12076, 14733, 9339, 9348, 12136, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_9(buffer, 14763, 0, 3, 12079, 14736, 9348, 9381, 12172, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_10(buffer, 14772, 0, 3, 12094, 14739, 9381, 9414, 12208, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_11(buffer, 14781, 0, 3, 12109, 14742, 9414, 9453, 12244, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_5(buffer, 14790, 0, 3, 14730, 14733, 12136, 14763, 9495, 9522, 12298, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_6(buffer, 14820, 0, 3, 14733, 14736, 12172, 14772, 9522, 9561, 12343, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_7(buffer, 14850, 0, 3, 14736, 14739, 12208, 14781, 9561, 9652, 12430, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_2(buffer, 14880, 0, 3, 14745, 14754, 12253, 14790, 9760, 9808, 12460, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_3(buffer, 14949, 0, 3, 14754, 14763, 12298, 14820, 9808, 9856, 12526, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_4(buffer, 15018, 0, 3, 14763, 14772, 12343, 14850, 9856, 10072, 12778, ncols, alpha, beta, p);

            compute_prim_il_electron_repulsion_1(buffer, 15156, 0, 3, 14790, 14820, 12526, 15018, 10357, 10588, 13551, ncols, alpha, beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 15609, 0, 3, 14880, 14949, 12913, 15156, 11134, 11422, 13947, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 17229, 15609, 1620, ncols);
        }
    }

    simdtrf::transform_kl(values, nvalues, buffer, 17229, nmax);
}

}  // namespace simdt2ceri
