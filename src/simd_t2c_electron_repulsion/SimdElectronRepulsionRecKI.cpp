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


#include "SimdElectronRepulsionRecKI.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ki_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ki_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 37181, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13}, ncols, fj, mu);

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

            compute_prim_ds_electron_repulsion_0(buffer, 59, 0, 7, 8, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 65, 0, 8, 9, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 71, 0, 9, 10, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 77, 0, 10, 11, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 83, 0, 11, 12, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 89, 0, 12, 13, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 95, 0, 13, 14, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 101, 0, 14, 15, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 107, 0, 15, 16, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 113, 0, 16, 17, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 119, 0, 17, 18, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 20, 23, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 23, 26, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 26, 29, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 29, 32, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 32, 35, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 35, 38, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 38, 41, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 195, 0, 41, 44, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 205, 0, 44, 47, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 215, 0, 47, 50, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 225, 0, 50, 53, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 235, 0, 59, 65, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 250, 0, 65, 71, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 265, 0, 71, 77, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 280, 0, 77, 83, 175, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 295, 0, 83, 89, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 310, 0, 89, 95, 195, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 325, 0, 95, 101, 205, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 340, 0, 101, 107, 215, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 355, 0, 107, 113, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 370, 0, 125, 135, 235, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 391, 0, 135, 145, 250, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 412, 0, 145, 155, 265, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 433, 0, 155, 165, 280, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 454, 0, 165, 175, 295, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 475, 0, 175, 185, 310, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 496, 0, 185, 195, 325, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 517, 0, 195, 205, 340, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 538, 0, 205, 215, 355, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 559, 0, 235, 250, 412, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 587, 0, 250, 265, 433, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 615, 0, 265, 280, 454, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 643, 0, 280, 295, 475, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 671, 0, 295, 310, 496, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 699, 0, 310, 325, 517, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 727, 0, 325, 340, 538, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 755, 0, 370, 391, 559, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 791, 0, 391, 412, 587, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 827, 0, 412, 433, 615, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 863, 0, 433, 454, 643, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 899, 0, 454, 475, 671, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 935, 0, 475, 496, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 971, 0, 496, 517, 727, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1007, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1010, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1013, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1016, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1019, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1022, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1025, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1028, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1031, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1034, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1037, 3, 9, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1046, 3, 10, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1055, 3, 11, 35, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1064, 3, 12, 38, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1073, 3, 13, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1082, 3, 14, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1091, 3, 15, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1100, 3, 16, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1109, 3, 17, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1118, 3, 18, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1127, 0, 3, 26, 1037, 65, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1145, 0, 3, 29, 1046, 71, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1163, 0, 3, 32, 1055, 77, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1181, 0, 3, 35, 1064, 83, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1199, 0, 3, 38, 1073, 89, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1217, 0, 3, 41, 1082, 95, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1235, 0, 3, 44, 1091, 101, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1253, 0, 3, 47, 1100, 107, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1271, 0, 3, 50, 1109, 113, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1289, 0, 3, 53, 1118, 119, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1307, 0, 3, 65, 1145, 145, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1337, 0, 3, 71, 1163, 155, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1367, 0, 3, 77, 1181, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1397, 0, 3, 83, 1199, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1427, 0, 3, 89, 1217, 185, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1457, 0, 3, 95, 1235, 195, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1487, 0, 3, 101, 1253, 205, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1517, 0, 3, 107, 1271, 215, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1547, 0, 3, 113, 1289, 225, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1577, 0, 3, 145, 1337, 250, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1622, 0, 3, 155, 1367, 265, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1667, 0, 3, 165, 1397, 280, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1712, 0, 3, 175, 1427, 295, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1757, 0, 3, 185, 1457, 310, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1802, 0, 3, 195, 1487, 325, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1847, 0, 3, 205, 1517, 340, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1892, 0, 3, 215, 1547, 355, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1937, 0, 3, 250, 1622, 412, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2000, 0, 3, 265, 1667, 433, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2063, 0, 3, 280, 1712, 454, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2126, 0, 3, 295, 1757, 475, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2189, 0, 3, 310, 1802, 496, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2252, 0, 3, 325, 1847, 517, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2315, 0, 3, 340, 1892, 538, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2378, 0, 3, 412, 2000, 587, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2462, 0, 3, 433, 2063, 615, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2546, 0, 3, 454, 2126, 643, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2630, 0, 3, 475, 2189, 671, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2714, 0, 3, 496, 2252, 699, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2798, 0, 3, 517, 2315, 727, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2882, 0, 3, 587, 2462, 827, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 2990, 0, 3, 615, 2546, 863, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3098, 0, 3, 643, 2630, 899, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3206, 0, 3, 671, 2714, 935, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3314, 0, 3, 699, 2798, 971, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3422, 3, 9, 10, 1010, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3428, 3, 10, 11, 1013, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3434, 3, 11, 12, 1016, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3440, 3, 12, 13, 1019, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3446, 3, 13, 14, 1022, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3452, 3, 14, 15, 1025, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3458, 3, 15, 16, 1028, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3464, 3, 16, 17, 1031, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3470, 3, 17, 18, 1034, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3476, 0, 3, 1007, 3422, 1046, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3494, 0, 3, 1010, 3428, 1055, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3512, 0, 3, 1013, 3434, 1064, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3530, 0, 3, 1016, 3440, 1073, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3548, 0, 3, 1019, 3446, 1082, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3566, 0, 3, 1022, 3452, 1091, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3584, 0, 3, 1025, 3458, 1100, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3602, 0, 3, 1028, 3464, 1109, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3620, 0, 3, 1031, 3470, 1118, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 3638, 0, 3, 1037, 3476, 59, 65, 1145,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3674, 0, 3, 1046, 3494, 65, 71, 1163,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3710, 0, 3, 1055, 3512, 71, 77, 1181,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3746, 0, 3, 1064, 3530, 77, 83, 1199,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3782, 0, 3, 1073, 3548, 83, 89, 1217,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3818, 0, 3, 1082, 3566, 89, 95, 1235,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3854, 0, 3, 1091, 3584, 95, 101, 1253,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3890, 0, 3, 1100, 3602, 101, 107, 1271,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3926, 0, 3, 1109, 3620, 107, 113, 1289,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3962, 0, 3, 1127, 3638, 125, 135, 1307,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4022, 0, 3, 1145, 3674, 135, 145, 1337,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4082, 0, 3, 1163, 3710, 145, 155, 1367,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4142, 0, 3, 1181, 3746, 155, 165, 1397,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4202, 0, 3, 1199, 3782, 165, 175, 1427,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4262, 0, 3, 1217, 3818, 175, 185, 1457,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4322, 0, 3, 1235, 3854, 185, 195, 1487,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4382, 0, 3, 1253, 3890, 195, 205, 1517,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4442, 0, 3, 1271, 3926, 205, 215, 1547,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4502, 0, 3, 3638, 3674, 1337, 4082, 235,
                                                 250, 1622, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4592, 0, 3, 3674, 3710, 1367, 4142, 250,
                                                 265, 1667, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4682, 0, 3, 3710, 3746, 1397, 4202, 265,
                                                 280, 1712, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4772, 0, 3, 3746, 3782, 1427, 4262, 280,
                                                 295, 1757, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4862, 0, 3, 3782, 3818, 1457, 4322, 295,
                                                 310, 1802, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4952, 0, 3, 3818, 3854, 1487, 4382, 310,
                                                 325, 1847, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5042, 0, 3, 3854, 3890, 1517, 4442, 325,
                                                 340, 1892, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5132, 0, 3, 3962, 4022, 1577, 4502, 370,
                                                 391, 1937, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5258, 0, 3, 4022, 4082, 1622, 4592, 391,
                                                 412, 2000, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5384, 0, 3, 4082, 4142, 1667, 4682, 412,
                                                 433, 2063, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5510, 0, 3, 4142, 4202, 1712, 4772, 433,
                                                 454, 2126, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5636, 0, 3, 4202, 4262, 1757, 4862, 454,
                                                 475, 2189, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5762, 0, 3, 4262, 4322, 1802, 4952, 475,
                                                 496, 2252, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5888, 0, 3, 4322, 4382, 1847, 5042, 496,
                                                 517, 2315, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6014, 0, 3, 4502, 4592, 2000, 5384, 559,
                                                 587, 2462, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6182, 0, 3, 4592, 4682, 2063, 5510, 587,
                                                 615, 2546, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6350, 0, 3, 4682, 4772, 2126, 5636, 615,
                                                 643, 2630, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6518, 0, 3, 4772, 4862, 2189, 5762, 643,
                                                 671, 2714, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6686, 0, 3, 4862, 4952, 2252, 5888, 671,
                                                 699, 2798, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 6854, 0, 3, 5132, 5258, 2378, 6014, 755,
                                                 791, 2882, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7070, 0, 3, 5258, 5384, 2462, 6182, 791,
                                                 827, 2990, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7286, 0, 3, 5384, 5510, 2546, 6350, 827,
                                                 863, 3098, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7502, 0, 3, 5510, 5636, 2630, 6518, 863,
                                                 899, 3206, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 7718, 0, 3, 5636, 5762, 2714, 6686, 899,
                                                 935, 3314, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7934, 3, 1007, 1010, 3428, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7944, 3, 1010, 1013, 3434, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7954, 3, 1013, 1016, 3440, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7964, 3, 1016, 1019, 3446, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7974, 3, 1019, 1022, 3452, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7984, 3, 1022, 1025, 3458, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7994, 3, 1025, 1028, 3464, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 8004, 3, 1028, 1031, 3470, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 8014, 0, 3, 3422, 7934, 3494, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8044, 0, 3, 3428, 7944, 3512, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8074, 0, 3, 3434, 7954, 3530, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8104, 0, 3, 3440, 7964, 3548, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8134, 0, 3, 3446, 7974, 3566, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8164, 0, 3, 3452, 7984, 3584, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8194, 0, 3, 3458, 7994, 3602, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 8224, 0, 3, 3464, 8004, 3620, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 8254, 0, 3, 3476, 8014, 1127, 1145,
                                                 3674, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8314, 0, 3, 3494, 8044, 1145, 1163,
                                                 3710, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8374, 0, 3, 3512, 8074, 1163, 1181,
                                                 3746, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8434, 0, 3, 3530, 8104, 1181, 1199,
                                                 3782, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8494, 0, 3, 3548, 8134, 1199, 1217,
                                                 3818, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8554, 0, 3, 3566, 8164, 1217, 1235,
                                                 3854, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8614, 0, 3, 3584, 8194, 1235, 1253,
                                                 3890, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 8674, 0, 3, 3602, 8224, 1253, 1271,
                                                 3926, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8734, 0, 3, 3674, 8314, 1307, 1337,
                                                 4082, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8834, 0, 3, 3710, 8374, 1337, 1367,
                                                 4142, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8934, 0, 3, 3746, 8434, 1367, 1397,
                                                 4202, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9034, 0, 3, 3782, 8494, 1397, 1427,
                                                 4262, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9134, 0, 3, 3818, 8554, 1427, 1457,
                                                 4322, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9234, 0, 3, 3854, 8614, 1457, 1487,
                                                 4382, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 9334, 0, 3, 3890, 8674, 1487, 1517,
                                                 4442, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9434, 0, 3, 8254, 8314, 4082, 8834,
                                                 1577, 1622, 4592, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9584, 0, 3, 8314, 8374, 4142, 8934,
                                                 1622, 1667, 4682, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9734, 0, 3, 8374, 8434, 4202, 9034,
                                                 1667, 1712, 4772, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9884, 0, 3, 8434, 8494, 4262, 9134,
                                                 1712, 1757, 4862, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10034, 0, 3, 8494, 8554, 4322, 9234,
                                                 1757, 1802, 4952, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 10184, 0, 3, 8554, 8614, 4382, 9334,
                                                 1802, 1847, 5042, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10334, 0, 3, 8734, 8834, 4592, 9584,
                                                 1937, 2000, 5384, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10544, 0, 3, 8834, 8934, 4682, 9734,
                                                 2000, 2063, 5510, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10754, 0, 3, 8934, 9034, 4772, 9884,
                                                 2063, 2126, 5636, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10964, 0, 3, 9034, 9134, 4862, 10034,
                                                 2126, 2189, 5762, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 11174, 0, 3, 9134, 9234, 4952, 10184,
                                                 2189, 2252, 5888, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11384, 0, 3, 9434, 9584, 5384, 10544,
                                                 2378, 2462, 6182, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11664, 0, 3, 9584, 9734, 5510, 10754,
                                                 2462, 2546, 6350, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11944, 0, 3, 9734, 9884, 5636, 10964,
                                                 2546, 2630, 6518, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 12224, 0, 3, 9884, 10034, 5762, 11174,
                                                 2630, 2714, 6686, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 12504, 0, 3, 10334, 10544, 6182, 11664,
                                                 2882, 2990, 7286, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 12864, 0, 3, 10544, 10754, 6350, 11944,
                                                 2990, 3098, 7502, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 13224, 0, 3, 10754, 10964, 6518, 12224,
                                                 3098, 3206, 7718, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13584, 3, 3422, 3428, 7944, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13599, 3, 3428, 3434, 7954, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13614, 3, 3434, 3440, 7964, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13629, 3, 3440, 3446, 7974, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13644, 3, 3446, 3452, 7984, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13659, 3, 3452, 3458, 7994, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 13674, 3, 3458, 3464, 8004, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 13689, 0, 3, 7934, 13584, 8044, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13734, 0, 3, 7944, 13599, 8074, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13779, 0, 3, 7954, 13614, 8104, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13824, 0, 3, 7964, 13629, 8134, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13869, 0, 3, 7974, 13644, 8164, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13914, 0, 3, 7984, 13659, 8194, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13959, 0, 3, 7994, 13674, 8224, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 14004, 0, 3, 8014, 13689, 3638, 3674,
                                                 8314, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 14094, 0, 3, 8044, 13734, 3674, 3710,
                                                 8374, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 14184, 0, 3, 8074, 13779, 3710, 3746,
                                                 8434, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 14274, 0, 3, 8104, 13824, 3746, 3782,
                                                 8494, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 14364, 0, 3, 8134, 13869, 3782, 3818,
                                                 8554, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 14454, 0, 3, 8164, 13914, 3818, 3854,
                                                 8614, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 14544, 0, 3, 8194, 13959, 3854, 3890,
                                                 8674, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14634, 0, 3, 8254, 14004, 3962, 4022,
                                                 8734, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14784, 0, 3, 8314, 14094, 4022, 4082,
                                                 8834, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14934, 0, 3, 8374, 14184, 4082, 4142,
                                                 8934, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 15084, 0, 3, 8434, 14274, 4142, 4202,
                                                 9034, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 15234, 0, 3, 8494, 14364, 4202, 4262,
                                                 9134, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 15384, 0, 3, 8554, 14454, 4262, 4322,
                                                 9234, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 15534, 0, 3, 8614, 14544, 4322, 4382,
                                                 9334, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15684, 0, 3, 14004, 14094, 8834, 14934,
                                                 4502, 4592, 9584, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15909, 0, 3, 14094, 14184, 8934, 15084,
                                                 4592, 4682, 9734, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 16134, 0, 3, 14184, 14274, 9034, 15234,
                                                 4682, 4772, 9884, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 16359, 0, 3, 14274, 14364, 9134, 15384,
                                                 4772, 4862, 10034, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 16584, 0, 3, 14364, 14454, 9234, 15534,
                                                 4862, 4952, 10184, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 16809, 0, 3, 14634, 14784, 9434, 15684,
                                                 5132, 5258, 10334, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 17124, 0, 3, 14784, 14934, 9584, 15909,
                                                 5258, 5384, 10544, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 17439, 0, 3, 14934, 15084, 9734, 16134,
                                                 5384, 5510, 10754, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 17754, 0, 3, 15084, 15234, 9884, 16359,
                                                 5510, 5636, 10964, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 18069, 0, 3, 15234, 15384, 10034, 16584,
                                                 5636, 5762, 11174, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 18384, 0, 3, 15684, 15909, 10544, 17439,
                                                 6014, 6182, 11664, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 18804, 0, 3, 15909, 16134, 10754, 17754,
                                                 6182, 6350, 11944, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 19224, 0, 3, 16134, 16359, 10964, 18069,
                                                 6350, 6518, 12224, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 19644, 0, 3, 16809, 17124, 11384, 18384,
                                                 6854, 7070, 12504, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 20184, 0, 3, 17124, 17439, 11664, 18804,
                                                 7070, 7286, 12864, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 20724, 0, 3, 17439, 17754, 11944, 19224,
                                                 7286, 7502, 13224, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21264, 3, 7934, 7944, 13599, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21285, 3, 7944, 7954, 13614, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21306, 3, 7954, 7964, 13629, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21327, 3, 7964, 7974, 13644, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21348, 3, 7974, 7984, 13659, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 21369, 3, 7984, 7994, 13674, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 21390, 0, 3, 13584, 21264, 13734, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21453, 0, 3, 13599, 21285, 13779, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21516, 0, 3, 13614, 21306, 13824, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21579, 0, 3, 13629, 21327, 13869, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21642, 0, 3, 13644, 21348, 13914, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 21705, 0, 3, 13659, 21369, 13959, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 21768, 0, 3, 13689, 21390, 8254, 8314,
                                                 14094, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 21894, 0, 3, 13734, 21453, 8314, 8374,
                                                 14184, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22020, 0, 3, 13779, 21516, 8374, 8434,
                                                 14274, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22146, 0, 3, 13824, 21579, 8434, 8494,
                                                 14364, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22272, 0, 3, 13869, 21642, 8494, 8554,
                                                 14454, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 22398, 0, 3, 13914, 21705, 8554, 8614,
                                                 14544, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 22524, 0, 3, 14094, 21894, 8734, 8834,
                                                 14934, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 22734, 0, 3, 14184, 22020, 8834, 8934,
                                                 15084, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 22944, 0, 3, 14274, 22146, 8934, 9034,
                                                 15234, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 23154, 0, 3, 14364, 22272, 9034, 9134,
                                                 15384, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 23364, 0, 3, 14454, 22398, 9134, 9234,
                                                 15534, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 23574, 0, 3, 21768, 21894, 14934, 22734,
                                                 9434, 9584, 15909, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 23889, 0, 3, 21894, 22020, 15084, 22944,
                                                 9584, 9734, 16134, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 24204, 0, 3, 22020, 22146, 15234, 23154,
                                                 9734, 9884, 16359, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 24519, 0, 3, 22146, 22272, 15384, 23364,
                                                 9884, 10034, 16584, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 24834, 0, 3, 22524, 22734, 15909, 23889,
                                                 10334, 10544, 17439, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 25275, 0, 3, 22734, 22944, 16134, 24204,
                                                 10544, 10754, 17754, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 25716, 0, 3, 22944, 23154, 16359, 24519,
                                                 10754, 10964, 18069, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 26157, 0, 3, 23574, 23889, 17439, 25275,
                                                 11384, 11664, 18804, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 26745, 0, 3, 23889, 24204, 17754, 25716,
                                                 11664, 11944, 19224, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 27333, 0, 3, 24834, 25275, 18804, 26745,
                                                 12504, 12864, 20724, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28089, 3, 13584, 13599, 21285, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28117, 3, 13599, 13614, 21306, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28145, 3, 13614, 13629, 21327, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28173, 3, 13629, 13644, 21348, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28201, 3, 13644, 13659, 21369, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 28229, 0, 3, 21264, 28089, 21453, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28313, 0, 3, 21285, 28117, 21516, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28397, 0, 3, 21306, 28145, 21579, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28481, 0, 3, 21327, 28173, 21642, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28565, 0, 3, 21348, 28201, 21705, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 28649, 0, 3, 21390, 28229, 14004, 14094,
                                                 21894, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 28817, 0, 3, 21453, 28313, 14094, 14184,
                                                 22020, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 28985, 0, 3, 21516, 28397, 14184, 14274,
                                                 22146, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 29153, 0, 3, 21579, 28481, 14274, 14364,
                                                 22272, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 29321, 0, 3, 21642, 28565, 14364, 14454,
                                                 22398, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 29489, 0, 3, 21768, 28649, 14634, 14784,
                                                 22524, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 29769, 0, 3, 21894, 28817, 14784, 14934,
                                                 22734, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 30049, 0, 3, 22020, 28985, 14934, 15084,
                                                 22944, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 30329, 0, 3, 22146, 29153, 15084, 15234,
                                                 23154, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 30609, 0, 3, 22272, 29321, 15234, 15384,
                                                 23364, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 30889, 0, 3, 28649, 28817, 22734, 30049,
                                                 15684, 15909, 23889, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 31309, 0, 3, 28817, 28985, 22944, 30329,
                                                 15909, 16134, 24204, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 31729, 0, 3, 28985, 29153, 23154, 30609,
                                                 16134, 16359, 24519, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 32149, 0, 3, 29489, 29769, 23574, 30889,
                                                 16809, 17124, 24834, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 32737, 0, 3, 29769, 30049, 23889, 31309,
                                                 17124, 17439, 25275, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 33325, 0, 3, 30049, 30329, 24204, 31729,
                                                 17439, 17754, 25716, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 33913, 0, 3, 30889, 31309, 25275, 33325,
                                                 18384, 18804, 26745, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 34697, 0, 3, 32149, 32737, 26157, 33913,
                                                 19644, 20184, 27333, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 35705, 34697, 1008, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 36713, 35705, 36, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 36713, 13, nmax);
}

}  // namespace simdt2ceri
