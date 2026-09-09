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


#include "SimdElectronRepulsionRecIL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPL.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSL.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_il_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_il_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 54745, 53009, 1260, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 14, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 61, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 79, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 85, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 91, 0, 12, 13, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 97, 0, 13, 14, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 103, 0, 14, 15, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 109, 0, 15, 16, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 115, 0, 16, 17, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 121, 0, 17, 18, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 127, 0, 18, 19, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 133, 0, 19, 20, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 139, 0, 22, 25, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 25, 28, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 28, 31, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 169, 0, 31, 34, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 179, 0, 34, 37, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 189, 0, 37, 40, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 199, 0, 40, 43, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 209, 0, 43, 46, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 219, 0, 46, 49, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 229, 0, 49, 52, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 239, 0, 52, 55, 133, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 249, 0, 61, 67, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 264, 0, 67, 73, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 279, 0, 73, 79, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 294, 0, 79, 85, 169, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 309, 0, 85, 91, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 324, 0, 91, 97, 189, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 339, 0, 97, 103, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 354, 0, 103, 109, 209, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 369, 0, 109, 115, 219, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 384, 0, 115, 121, 229, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 399, 0, 121, 127, 239, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 414, 0, 139, 149, 279, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 435, 0, 149, 159, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 456, 0, 159, 169, 309, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 477, 0, 169, 179, 324, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 498, 0, 179, 189, 339, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 519, 0, 189, 199, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 540, 0, 199, 209, 369, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 561, 0, 209, 219, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 582, 0, 219, 229, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 603, 0, 249, 264, 414, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 631, 0, 264, 279, 435, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 659, 0, 279, 294, 456, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 687, 0, 294, 309, 477, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 715, 0, 309, 324, 498, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 743, 0, 324, 339, 519, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 771, 0, 339, 354, 540, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 799, 0, 354, 369, 561, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 827, 0, 369, 384, 582, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 855, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 858, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 861, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 864, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 867, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 870, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 873, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 876, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 879, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 882, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 885, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 888, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 891, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 900, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 909, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 918, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 927, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 936, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 945, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 954, 3, 16, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 963, 3, 17, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 972, 3, 18, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 981, 3, 19, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 990, 3, 20, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 999, 0, 3, 25, 900, 73, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1017, 0, 3, 28, 909, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1035, 0, 3, 31, 918, 85, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1053, 0, 3, 34, 927, 91, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1071, 0, 3, 37, 936, 97, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1089, 0, 3, 40, 945, 103, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1107, 0, 3, 43, 954, 109, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1125, 0, 3, 46, 963, 115, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1143, 0, 3, 49, 972, 121, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1161, 0, 3, 52, 981, 127, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1179, 0, 3, 55, 990, 133, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1197, 0, 3, 73, 1017, 149, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1227, 0, 3, 79, 1035, 159, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1257, 0, 3, 85, 1053, 169, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1287, 0, 3, 91, 1071, 179, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1317, 0, 3, 97, 1089, 189, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1347, 0, 3, 103, 1107, 199, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1377, 0, 3, 109, 1125, 209, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1407, 0, 3, 115, 1143, 219, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1437, 0, 3, 121, 1161, 229, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1467, 0, 3, 127, 1179, 239, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1497, 0, 3, 149, 1227, 279, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1542, 0, 3, 159, 1257, 294, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1587, 0, 3, 169, 1287, 309, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1632, 0, 3, 179, 1317, 324, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1677, 0, 3, 189, 1347, 339, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1722, 0, 3, 199, 1377, 354, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1767, 0, 3, 209, 1407, 369, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1812, 0, 3, 219, 1437, 384, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1857, 0, 3, 229, 1467, 399, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1902, 0, 3, 279, 1542, 435, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1965, 0, 3, 294, 1587, 456, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2028, 0, 3, 309, 1632, 477, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2091, 0, 3, 324, 1677, 498, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2154, 0, 3, 339, 1722, 519, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2217, 0, 3, 354, 1767, 540, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2280, 0, 3, 369, 1812, 561, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2343, 0, 3, 384, 1857, 582, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2406, 0, 3, 435, 1965, 659, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2490, 0, 3, 456, 2028, 687, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2574, 0, 3, 477, 2091, 715, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2658, 0, 3, 498, 2154, 743, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2742, 0, 3, 519, 2217, 771, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2826, 0, 3, 540, 2280, 799, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2910, 0, 3, 561, 2343, 827, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2994, 3, 9, 10, 858, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 3000, 3, 10, 11, 861, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3006, 3, 11, 12, 864, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3012, 3, 12, 13, 867, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3018, 3, 13, 14, 870, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3024, 3, 14, 15, 873, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3030, 3, 15, 16, 876, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3036, 3, 16, 17, 879, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3042, 3, 17, 18, 882, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3048, 3, 18, 19, 885, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3054, 3, 19, 20, 888, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3060, 0, 3, 855, 2994, 900, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3078, 0, 3, 858, 3000, 909, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3096, 0, 3, 861, 3006, 918, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3114, 0, 3, 864, 3012, 927, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3132, 0, 3, 867, 3018, 936, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3150, 0, 3, 870, 3024, 945, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3168, 0, 3, 873, 3030, 954, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3186, 0, 3, 876, 3036, 963, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3204, 0, 3, 879, 3042, 972, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3222, 0, 3, 882, 3048, 981, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3240, 0, 3, 885, 3054, 990, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3258, 0, 3, 891, 3060, 61, 67, 999,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3294, 0, 3, 900, 3078, 67, 73, 1017,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3330, 0, 3, 909, 3096, 73, 79, 1035,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3366, 0, 3, 918, 3114, 79, 85, 1053,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3402, 0, 3, 927, 3132, 85, 91, 1071,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3438, 0, 3, 936, 3150, 91, 97, 1089,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3474, 0, 3, 945, 3168, 97, 103, 1107,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3510, 0, 3, 954, 3186, 103, 109, 1125,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3546, 0, 3, 963, 3204, 109, 115, 1143,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3582, 0, 3, 972, 3222, 115, 121, 1161,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3618, 0, 3, 981, 3240, 121, 127, 1179,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3654, 0, 3, 1017, 3330, 139, 149, 1227,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3714, 0, 3, 1035, 3366, 149, 159, 1257,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3774, 0, 3, 1053, 3402, 159, 169, 1287,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3834, 0, 3, 1071, 3438, 169, 179, 1317,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3894, 0, 3, 1089, 3474, 179, 189, 1347,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3954, 0, 3, 1107, 3510, 189, 199, 1377,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4014, 0, 3, 1125, 3546, 199, 209, 1407,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4074, 0, 3, 1143, 3582, 209, 219, 1437,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4134, 0, 3, 1161, 3618, 219, 229, 1467,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4194, 0, 3, 3258, 3294, 1197, 3654, 249,
                                                 264, 1497, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4284, 0, 3, 3294, 3330, 1227, 3714, 264,
                                                 279, 1542, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4374, 0, 3, 3330, 3366, 1257, 3774, 279,
                                                 294, 1587, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4464, 0, 3, 3366, 3402, 1287, 3834, 294,
                                                 309, 1632, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4554, 0, 3, 3402, 3438, 1317, 3894, 309,
                                                 324, 1677, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4644, 0, 3, 3438, 3474, 1347, 3954, 324,
                                                 339, 1722, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4734, 0, 3, 3474, 3510, 1377, 4014, 339,
                                                 354, 1767, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4824, 0, 3, 3510, 3546, 1407, 4074, 354,
                                                 369, 1812, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4914, 0, 3, 3546, 3582, 1437, 4134, 369,
                                                 384, 1857, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5004, 0, 3, 3654, 3714, 1542, 4374, 414,
                                                 435, 1965, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5130, 0, 3, 3714, 3774, 1587, 4464, 435,
                                                 456, 2028, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5256, 0, 3, 3774, 3834, 1632, 4554, 456,
                                                 477, 2091, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5382, 0, 3, 3834, 3894, 1677, 4644, 477,
                                                 498, 2154, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5508, 0, 3, 3894, 3954, 1722, 4734, 498,
                                                 519, 2217, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5634, 0, 3, 3954, 4014, 1767, 4824, 519,
                                                 540, 2280, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 5760, 0, 3, 4014, 4074, 1812, 4914, 540,
                                                 561, 2343, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 5886, 0, 3, 4194, 4284, 1902, 5004, 603,
                                                 631, 2406, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6054, 0, 3, 4284, 4374, 1965, 5130, 631,
                                                 659, 2490, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6222, 0, 3, 4374, 4464, 2028, 5256, 659,
                                                 687, 2574, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6390, 0, 3, 4464, 4554, 2091, 5382, 687,
                                                 715, 2658, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6558, 0, 3, 4554, 4644, 2154, 5508, 715,
                                                 743, 2742, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6726, 0, 3, 4644, 4734, 2217, 5634, 743,
                                                 771, 2826, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 6894, 0, 3, 4734, 4824, 2280, 5760, 771,
                                                 799, 2910, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7062, 3, 855, 858, 3000, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7072, 3, 858, 861, 3006, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7082, 3, 861, 864, 3012, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7092, 3, 864, 867, 3018, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7102, 3, 867, 870, 3024, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7112, 3, 870, 873, 3030, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7122, 3, 873, 876, 3036, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7132, 3, 876, 879, 3042, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7142, 3, 879, 882, 3048, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 7152, 3, 882, 885, 3054, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 7162, 0, 3, 2994, 7062, 3078, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7192, 0, 3, 3000, 7072, 3096, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7222, 0, 3, 3006, 7082, 3114, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7252, 0, 3, 3012, 7092, 3132, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7282, 0, 3, 3018, 7102, 3150, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7312, 0, 3, 3024, 7112, 3168, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7342, 0, 3, 3030, 7122, 3186, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7372, 0, 3, 3036, 7132, 3204, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7402, 0, 3, 3042, 7142, 3222, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 7432, 0, 3, 3048, 7152, 3240, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 7462, 0, 3, 3078, 7192, 999, 1017, 3330,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7522, 0, 3, 3096, 7222, 1017, 1035,
                                                 3366, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7582, 0, 3, 3114, 7252, 1035, 1053,
                                                 3402, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7642, 0, 3, 3132, 7282, 1053, 1071,
                                                 3438, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7702, 0, 3, 3150, 7312, 1071, 1089,
                                                 3474, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7762, 0, 3, 3168, 7342, 1089, 1107,
                                                 3510, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7822, 0, 3, 3186, 7372, 1107, 1125,
                                                 3546, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7882, 0, 3, 3204, 7402, 1125, 1143,
                                                 3582, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 7942, 0, 3, 3222, 7432, 1143, 1161,
                                                 3618, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8002, 0, 3, 3330, 7522, 1197, 1227,
                                                 3714, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8102, 0, 3, 3366, 7582, 1227, 1257,
                                                 3774, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8202, 0, 3, 3402, 7642, 1257, 1287,
                                                 3834, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8302, 0, 3, 3438, 7702, 1287, 1317,
                                                 3894, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8402, 0, 3, 3474, 7762, 1317, 1347,
                                                 3954, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8502, 0, 3, 3510, 7822, 1347, 1377,
                                                 4014, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8602, 0, 3, 3546, 7882, 1377, 1407,
                                                 4074, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8702, 0, 3, 3582, 7942, 1407, 1437,
                                                 4134, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8802, 0, 3, 7462, 7522, 3714, 8102,
                                                 1497, 1542, 4374, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8952, 0, 3, 7522, 7582, 3774, 8202,
                                                 1542, 1587, 4464, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9102, 0, 3, 7582, 7642, 3834, 8302,
                                                 1587, 1632, 4554, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9252, 0, 3, 7642, 7702, 3894, 8402,
                                                 1632, 1677, 4644, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9402, 0, 3, 7702, 7762, 3954, 8502,
                                                 1677, 1722, 4734, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9552, 0, 3, 7762, 7822, 4014, 8602,
                                                 1722, 1767, 4824, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9702, 0, 3, 7822, 7882, 4074, 8702,
                                                 1767, 1812, 4914, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 9852, 0, 3, 8002, 8102, 4374, 8952,
                                                 1902, 1965, 5130, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10062, 0, 3, 8102, 8202, 4464, 9102,
                                                 1965, 2028, 5256, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10272, 0, 3, 8202, 8302, 4554, 9252,
                                                 2028, 2091, 5382, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10482, 0, 3, 8302, 8402, 4644, 9402,
                                                 2091, 2154, 5508, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10692, 0, 3, 8402, 8502, 4734, 9552,
                                                 2154, 2217, 5634, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 10902, 0, 3, 8502, 8602, 4824, 9702,
                                                 2217, 2280, 5760, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11112, 0, 3, 8802, 8952, 5130, 10062,
                                                 2406, 2490, 6222, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11392, 0, 3, 8952, 9102, 5256, 10272,
                                                 2490, 2574, 6390, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11672, 0, 3, 9102, 9252, 5382, 10482,
                                                 2574, 2658, 6558, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 11952, 0, 3, 9252, 9402, 5508, 10692,
                                                 2658, 2742, 6726, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 12232, 0, 3, 9402, 9552, 5634, 10902,
                                                 2742, 2826, 6894, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12512, 3, 2994, 3000, 7072, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12527, 3, 3000, 3006, 7082, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12542, 3, 3006, 3012, 7092, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12557, 3, 3012, 3018, 7102, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12572, 3, 3018, 3024, 7112, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12587, 3, 3024, 3030, 7122, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12602, 3, 3030, 3036, 7132, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12617, 3, 3036, 3042, 7142, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 12632, 3, 3042, 3048, 7152, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 12647, 0, 3, 7062, 12512, 7192, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12692, 0, 3, 7072, 12527, 7222, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12737, 0, 3, 7082, 12542, 7252, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12782, 0, 3, 7092, 12557, 7282, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12827, 0, 3, 7102, 12572, 7312, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12872, 0, 3, 7112, 12587, 7342, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12917, 0, 3, 7122, 12602, 7372, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 12962, 0, 3, 7132, 12617, 7402, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 13007, 0, 3, 7142, 12632, 7432, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 13052, 0, 3, 7162, 12647, 3258, 3294,
                                                 7462, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13142, 0, 3, 7192, 12692, 3294, 3330,
                                                 7522, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13232, 0, 3, 7222, 12737, 3330, 3366,
                                                 7582, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13322, 0, 3, 7252, 12782, 3366, 3402,
                                                 7642, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13412, 0, 3, 7282, 12827, 3402, 3438,
                                                 7702, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13502, 0, 3, 7312, 12872, 3438, 3474,
                                                 7762, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13592, 0, 3, 7342, 12917, 3474, 3510,
                                                 7822, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13682, 0, 3, 7372, 12962, 3510, 3546,
                                                 7882, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 13772, 0, 3, 7402, 13007, 3546, 3582,
                                                 7942, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13862, 0, 3, 7522, 13232, 3654, 3714,
                                                 8102, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14012, 0, 3, 7582, 13322, 3714, 3774,
                                                 8202, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14162, 0, 3, 7642, 13412, 3774, 3834,
                                                 8302, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14312, 0, 3, 7702, 13502, 3834, 3894,
                                                 8402, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14462, 0, 3, 7762, 13592, 3894, 3954,
                                                 8502, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14612, 0, 3, 7822, 13682, 3954, 4014,
                                                 8602, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 14762, 0, 3, 7882, 13772, 4014, 4074,
                                                 8702, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14912, 0, 3, 13052, 13142, 8002, 13862,
                                                 4194, 4284, 8802, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15137, 0, 3, 13142, 13232, 8102, 14012,
                                                 4284, 4374, 8952, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15362, 0, 3, 13232, 13322, 8202, 14162,
                                                 4374, 4464, 9102, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15587, 0, 3, 13322, 13412, 8302, 14312,
                                                 4464, 4554, 9252, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15812, 0, 3, 13412, 13502, 8402, 14462,
                                                 4554, 4644, 9402, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 16037, 0, 3, 13502, 13592, 8502, 14612,
                                                 4644, 4734, 9552, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 16262, 0, 3, 13592, 13682, 8602, 14762,
                                                 4734, 4824, 9702, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 16487, 0, 3, 13862, 14012, 8952, 15362,
                                                 5004, 5130, 10062, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 16802, 0, 3, 14012, 14162, 9102, 15587,
                                                 5130, 5256, 10272, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 17117, 0, 3, 14162, 14312, 9252, 15812,
                                                 5256, 5382, 10482, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 17432, 0, 3, 14312, 14462, 9402, 16037,
                                                 5382, 5508, 10692, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 17747, 0, 3, 14462, 14612, 9552, 16262,
                                                 5508, 5634, 10902, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 18062, 0, 3, 14912, 15137, 9852, 16487,
                                                 5886, 6054, 11112, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 18482, 0, 3, 15137, 15362, 10062, 16802,
                                                 6054, 6222, 11392, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 18902, 0, 3, 15362, 15587, 10272, 17117,
                                                 6222, 6390, 11672, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 19322, 0, 3, 15587, 15812, 10482, 17432,
                                                 6390, 6558, 11952, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 19742, 0, 3, 15812, 16037, 10692, 17747,
                                                 6558, 6726, 12232, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 20162, 3, 7062, 7072, 12527, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 20183, 3, 7072, 7082, 12542, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 20204, 3, 7082, 7092, 12557, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 20225, 3, 7092, 7102, 12572, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 20246, 3, 7102, 7112, 12587, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 20267, 3, 7112, 7122, 12602, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 20288, 3, 7122, 7132, 12617, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 20309, 3, 7132, 7142, 12632, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 20330, 0, 3, 12512, 20162, 12692, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 20393, 0, 3, 12527, 20183, 12737, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 20456, 0, 3, 12542, 20204, 12782, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 20519, 0, 3, 12557, 20225, 12827, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 20582, 0, 3, 12572, 20246, 12872, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 20645, 0, 3, 12587, 20267, 12917, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 20708, 0, 3, 12602, 20288, 12962, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 20771, 0, 3, 12617, 20309, 13007, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 20834, 0, 3, 12692, 20393, 7462, 7522,
                                                 13232, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 20960, 0, 3, 12737, 20456, 7522, 7582,
                                                 13322, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 21086, 0, 3, 12782, 20519, 7582, 7642,
                                                 13412, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 21212, 0, 3, 12827, 20582, 7642, 7702,
                                                 13502, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 21338, 0, 3, 12872, 20645, 7702, 7762,
                                                 13592, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 21464, 0, 3, 12917, 20708, 7762, 7822,
                                                 13682, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 21590, 0, 3, 12962, 20771, 7822, 7882,
                                                 13772, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 21716, 0, 3, 13232, 20960, 8002, 8102,
                                                 14012, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 21926, 0, 3, 13322, 21086, 8102, 8202,
                                                 14162, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 22136, 0, 3, 13412, 21212, 8202, 8302,
                                                 14312, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 22346, 0, 3, 13502, 21338, 8302, 8402,
                                                 14462, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 22556, 0, 3, 13592, 21464, 8402, 8502,
                                                 14612, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 22766, 0, 3, 13682, 21590, 8502, 8602,
                                                 14762, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 22976, 0, 3, 20834, 20960, 14012, 21926,
                                                 8802, 8952, 15362, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 23291, 0, 3, 20960, 21086, 14162, 22136,
                                                 8952, 9102, 15587, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 23606, 0, 3, 21086, 21212, 14312, 22346,
                                                 9102, 9252, 15812, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 23921, 0, 3, 21212, 21338, 14462, 22556,
                                                 9252, 9402, 16037, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 24236, 0, 3, 21338, 21464, 14612, 22766,
                                                 9402, 9552, 16262, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 24551, 0, 3, 21716, 21926, 15362, 23291,
                                                 9852, 10062, 16802, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 24992, 0, 3, 21926, 22136, 15587, 23606,
                                                 10062, 10272, 17117, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 25433, 0, 3, 22136, 22346, 15812, 23921,
                                                 10272, 10482, 17432, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 25874, 0, 3, 22346, 22556, 16037, 24236,
                                                 10482, 10692, 17747, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 26315, 0, 3, 22976, 23291, 16802, 24992,
                                                 11112, 11392, 18902, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 26903, 0, 3, 23291, 23606, 17117, 25433,
                                                 11392, 11672, 19322, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 27491, 0, 3, 23606, 23921, 17432, 25874,
                                                 11672, 11952, 19742, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28079, 3, 12512, 12527, 20183, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28107, 3, 12527, 12542, 20204, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28135, 3, 12542, 12557, 20225, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28163, 3, 12557, 12572, 20246, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28191, 3, 12572, 12587, 20267, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28219, 3, 12587, 12602, 20288, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 28247, 3, 12602, 12617, 20309, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 28275, 0, 3, 20162, 28079, 20393, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28359, 0, 3, 20183, 28107, 20456, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28443, 0, 3, 20204, 28135, 20519, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28527, 0, 3, 20225, 28163, 20582, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28611, 0, 3, 20246, 28191, 20645, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28695, 0, 3, 20267, 28219, 20708, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 28779, 0, 3, 20288, 28247, 20771, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 28863, 0, 3, 20330, 28275, 13052, 13142,
                                                 20834, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 29031, 0, 3, 20393, 28359, 13142, 13232,
                                                 20960, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 29199, 0, 3, 20456, 28443, 13232, 13322,
                                                 21086, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 29367, 0, 3, 20519, 28527, 13322, 13412,
                                                 21212, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 29535, 0, 3, 20582, 28611, 13412, 13502,
                                                 21338, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 29703, 0, 3, 20645, 28695, 13502, 13592,
                                                 21464, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 29871, 0, 3, 20708, 28779, 13592, 13682,
                                                 21590, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 30039, 0, 3, 20960, 29199, 13862, 14012,
                                                 21926, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 30319, 0, 3, 21086, 29367, 14012, 14162,
                                                 22136, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 30599, 0, 3, 21212, 29535, 14162, 14312,
                                                 22346, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 30879, 0, 3, 21338, 29703, 14312, 14462,
                                                 22556, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 31159, 0, 3, 21464, 29871, 14462, 14612,
                                                 22766, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 31439, 0, 3, 28863, 29031, 21716, 30039,
                                                 14912, 15137, 22976, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 31859, 0, 3, 29031, 29199, 21926, 30319,
                                                 15137, 15362, 23291, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 32279, 0, 3, 29199, 29367, 22136, 30599,
                                                 15362, 15587, 23606, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 32699, 0, 3, 29367, 29535, 22346, 30879,
                                                 15587, 15812, 23921, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 33119, 0, 3, 29535, 29703, 22556, 31159,
                                                 15812, 16037, 24236, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 33539, 0, 3, 30039, 30319, 23291, 32279,
                                                 16487, 16802, 24992, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 34127, 0, 3, 30319, 30599, 23606, 32699,
                                                 16802, 17117, 25433, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 34715, 0, 3, 30599, 30879, 23921, 33119,
                                                 17117, 17432, 25874, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 35303, 0, 3, 31439, 31859, 24551, 33539,
                                                 18062, 18482, 26315, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 36087, 0, 3, 31859, 32279, 24992, 34127,
                                                 18482, 18902, 26903, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 36871, 0, 3, 32279, 32699, 25433, 34715,
                                                 18902, 19322, 27491, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 37655, 3, 20162, 20183, 28107, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 37691, 3, 20183, 20204, 28135, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 37727, 3, 20204, 20225, 28163, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 37763, 3, 20225, 20246, 28191, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 37799, 3, 20246, 20267, 28219, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 37835, 3, 20267, 20288, 28247, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 37871, 0, 3, 28079, 37655, 28359, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 37979, 0, 3, 28107, 37691, 28443, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 38087, 0, 3, 28135, 37727, 28527, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 38195, 0, 3, 28163, 37763, 28611, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 38303, 0, 3, 28191, 37799, 28695, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 38411, 0, 3, 28219, 37835, 28779, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 38519, 0, 3, 28359, 37979, 20834, 20960,
                                                 29199, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 38735, 0, 3, 28443, 38087, 20960, 21086,
                                                 29367, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 38951, 0, 3, 28527, 38195, 21086, 21212,
                                                 29535, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 39167, 0, 3, 28611, 38303, 21212, 21338,
                                                 29703, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 39383, 0, 3, 28695, 38411, 21338, 21464,
                                                 29871, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 39599, 0, 3, 29199, 38735, 21716, 21926,
                                                 30319, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 39959, 0, 3, 29367, 38951, 21926, 22136,
                                                 30599, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 40319, 0, 3, 29535, 39167, 22136, 22346,
                                                 30879, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 40679, 0, 3, 29703, 39383, 22346, 22556,
                                                 31159, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 41039, 0, 3, 38519, 38735, 30319, 39959,
                                                 22976, 23291, 32279, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 41579, 0, 3, 38735, 38951, 30599, 40319,
                                                 23291, 23606, 32699, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 42119, 0, 3, 38951, 39167, 30879, 40679,
                                                 23606, 23921, 33119, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 42659, 0, 3, 39599, 39959, 32279, 41579,
                                                 24551, 24992, 34127, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 43415, 0, 3, 39959, 40319, 32699, 42119,
                                                 24992, 25433, 34715, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 44171, 0, 3, 41039, 41579, 34127, 43415,
                                                 26315, 26903, 36871, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 45179, 3, 28079, 28107, 37691, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 45224, 3, 28107, 28135, 37727, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 45269, 3, 28135, 28163, 37763, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 45314, 3, 28163, 28191, 37799, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 45359, 3, 28191, 28219, 37835, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 45404, 0, 3, 37655, 45179, 37979, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 45539, 0, 3, 37691, 45224, 38087, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 45674, 0, 3, 37727, 45269, 38195, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 45809, 0, 3, 37763, 45314, 38303, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 45944, 0, 3, 37799, 45359, 38411, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 46079, 0, 3, 37871, 45404, 28863, 29031,
                                                 38519, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 46349, 0, 3, 37979, 45539, 29031, 29199,
                                                 38735, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 46619, 0, 3, 38087, 45674, 29199, 29367,
                                                 38951, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 46889, 0, 3, 38195, 45809, 29367, 29535,
                                                 39167, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 47159, 0, 3, 38303, 45944, 29535, 29703,
                                                 39383, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 47429, 0, 3, 38735, 46619, 30039, 30319,
                                                 39959, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 47879, 0, 3, 38951, 46889, 30319, 30599,
                                                 40319, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 48329, 0, 3, 39167, 47159, 30599, 30879,
                                                 40679, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 48779, 0, 3, 46079, 46349, 39599, 47429,
                                                 31439, 31859, 41039, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 49454, 0, 3, 46349, 46619, 39959, 47879,
                                                 31859, 32279, 41579, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 50129, 0, 3, 46619, 46889, 40319, 48329,
                                                 32279, 32699, 42119, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 50804, 0, 3, 47429, 47879, 41579, 50129,
                                                 33539, 34127, 43415, ncols, alpha, beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 51749, 0, 3, 48779, 49454, 42659, 50804,
                                                 35303, 36087, 44171, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 53009, 51749, 1260, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 54269, 53009, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 54269, 17, nmax);
}

}  // namespace simdt2ceri
