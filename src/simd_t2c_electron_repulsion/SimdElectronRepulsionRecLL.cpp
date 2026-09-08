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


#include "SimdElectronRepulsionRecLL.hpp"

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
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLK.hpp"
#include "SimdElectronRepulsionVrrRecLL.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
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
#include "SimdTransformLL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ll_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ll_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(25732, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 16, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 23, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 7, 8, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 72, 0, 8, 9, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 75, 0, 9, 10, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 78, 0, 10, 11, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 81, 0, 11, 12, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 84, 0, 12, 13, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 87, 0, 13, 14, 42, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 90, 0, 14, 15, 45, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 93, 0, 15, 16, 48, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 96, 0, 16, 17, 51, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 99, 0, 17, 18, 54, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 102, 0, 18, 19, 57, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 105, 0, 19, 20, 60, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 108, 0, 20, 21, 63, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 111, 0, 21, 22, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 114, 0, 24, 27, 75, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 120, 0, 27, 30, 78, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 126, 0, 30, 33, 81, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 132, 0, 33, 36, 84, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 138, 0, 36, 39, 87, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 144, 0, 39, 42, 90, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 150, 0, 42, 45, 93, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 156, 0, 45, 48, 96, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 162, 0, 48, 51, 99, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 168, 0, 51, 54, 102, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 174, 0, 54, 57, 105, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 180, 0, 57, 60, 108, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 186, 0, 60, 63, 111, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 192, 0, 69, 72, 114, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 201, 0, 72, 75, 120, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 210, 0, 75, 78, 126, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 219, 0, 78, 81, 132, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 228, 0, 81, 84, 138, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 237, 0, 84, 87, 144, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 246, 0, 87, 90, 150, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 255, 0, 90, 93, 156, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 264, 0, 93, 96, 162, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 273, 0, 96, 99, 168, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 282, 0, 99, 102, 174, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 291, 0, 102, 105, 180, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 300, 0, 105, 108, 186, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 309, 0, 114, 120, 210, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 321, 0, 120, 126, 219, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 333, 0, 126, 132, 228, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 345, 0, 132, 138, 237, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 357, 0, 138, 144, 246, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 369, 0, 144, 150, 255, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 381, 0, 150, 156, 264, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 393, 0, 156, 162, 273, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 405, 0, 162, 168, 282, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 417, 0, 168, 174, 291, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 429, 0, 174, 180, 300, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 441, 0, 192, 201, 309, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 456, 0, 201, 210, 321, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 471, 0, 210, 219, 333, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 486, 0, 219, 228, 345, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 501, 0, 228, 237, 357, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 516, 0, 237, 246, 369, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 531, 0, 246, 255, 381, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 546, 0, 255, 264, 393, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 561, 0, 264, 273, 405, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 576, 0, 273, 282, 417, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 591, 0, 282, 291, 429, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 606, 0, 309, 321, 471, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 624, 0, 321, 333, 486, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 642, 0, 333, 345, 501, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 660, 0, 345, 357, 516, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 678, 0, 357, 369, 531, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 696, 0, 369, 381, 546, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 714, 0, 381, 393, 561, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 732, 0, 393, 405, 576, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_15(buffer, 750, 0, 405, 417, 591, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 768, 0, 441, 456, 606, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 786, 0, 456, 471, 624, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 804, 0, 471, 486, 642, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 822, 0, 486, 501, 660, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 840, 0, 501, 516, 678, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 858, 0, 516, 531, 696, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 876, 0, 531, 546, 714, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 894, 0, 546, 561, 732, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 912, 0, 561, 576, 750, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 930, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 933, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 936, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 939, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 942, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 945, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 948, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 951, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 954, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 957, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 960, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 963, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 966, 3, 22, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 969, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 972, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 975, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 978, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 981, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 984, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 987, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 990, 3, 16, 48, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 993, 3, 17, 51, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 996, 3, 18, 54, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 999, 3, 19, 57, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 1002, 3, 20, 60, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 1005, 3, 21, 63, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1008, 3, 27, 75, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1011, 3, 30, 78, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1014, 3, 33, 81, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1017, 3, 36, 84, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1020, 3, 39, 87, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1023, 3, 42, 90, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1026, 3, 45, 93, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1029, 3, 48, 96, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1032, 3, 51, 99, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1035, 3, 54, 102, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1038, 3, 57, 105, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1041, 3, 60, 108, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 1044, 3, 63, 111, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1047, 3, 75, 120, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1050, 3, 78, 126, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1053, 3, 81, 132, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1056, 3, 84, 138, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1059, 3, 87, 144, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1062, 3, 90, 150, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1065, 3, 93, 156, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1068, 3, 96, 162, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1071, 3, 99, 168, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1074, 3, 102, 174, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1077, 3, 105, 180, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 1080, 3, 108, 186, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 1083, 3, 120, 210, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 1086, 3, 126, 219, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 1089, 3, 132, 228, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 1092, 3, 138, 237, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 1101, 3, 144, 246, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 1110, 3, 150, 255, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 1119, 3, 156, 264, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 1128, 3, 162, 273, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 1137, 3, 168, 282, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 1146, 3, 174, 291, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 1149, 3, 180, 300, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 1152, 3, 210, 321, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 1155, 3, 219, 333, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1158, 3, 228, 345, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1173, 3, 237, 357, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1188, 3, 246, 369, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1203, 3, 255, 381, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1218, 3, 264, 393, ncols, p);

            compute_prim_hp_electron_repulsion_8(buffer, 1233, 3, 273, 405, ncols, p);

            compute_prim_hp_electron_repulsion_17(buffer, 1248, 3, 282, 417, ncols, p);

            compute_prim_hp_electron_repulsion_15(buffer, 1257, 3, 291, 429, ncols, p);

            compute_prim_ip_electron_repulsion_15(buffer, 1260, 3, 321, 471, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1263, 3, 333, 486, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1281, 3, 345, 501, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1299, 3, 357, 516, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1317, 3, 369, 531, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1335, 3, 381, 546, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1353, 3, 393, 561, ncols, p);

            compute_prim_ip_electron_repulsion_8(buffer, 1371, 3, 405, 576, ncols, p);

            compute_prim_ip_electron_repulsion_17(buffer, 1389, 3, 417, 591, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1398, 3, 471, 624, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1419, 3, 486, 642, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1440, 3, 501, 660, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1461, 3, 516, 678, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1482, 3, 531, 696, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1503, 3, 546, 714, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1524, 3, 561, 732, ncols, p);

            compute_prim_kp_electron_repulsion_8(buffer, 1545, 3, 576, 750, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1566, 3, 624, 804, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1590, 3, 642, 822, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1614, 3, 660, 840, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1638, 3, 678, 858, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1662, 3, 696, 876, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1686, 3, 714, 894, ncols, p);

            compute_prim_lp_electron_repulsion_3(buffer, 1710, 3, 732, 912, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1734, 3, 9, 10, 933, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1737, 3, 10, 11, 936, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1740, 3, 11, 12, 939, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1743, 3, 12, 13, 942, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1746, 3, 13, 14, 945, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1749, 3, 14, 15, 948, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1752, 3, 15, 16, 951, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1755, 3, 16, 17, 954, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1758, 3, 17, 18, 957, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1761, 3, 18, 19, 960, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1764, 3, 19, 20, 963, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 1767, 3, 20, 21, 966, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1770, 0, 933, 1737, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1773, 0, 936, 1740, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1776, 0, 939, 1743, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1779, 0, 942, 1746, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1782, 0, 945, 1749, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1785, 0, 948, 1752, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1788, 0, 951, 1755, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1791, 0, 954, 1758, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1794, 0, 957, 1761, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1797, 0, 960, 1764, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 1800, 0, 963, 1767, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1803, 3, 969, 69, 72, 1008, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1806, 3, 972, 72, 75, 1011, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1809, 3, 975, 75, 78, 1014, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1812, 3, 978, 78, 81, 1017, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1815, 3, 981, 81, 84, 1020, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1818, 3, 984, 84, 87, 1023, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1821, 3, 987, 87, 90, 1026, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1824, 3, 990, 90, 93, 1029, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1827, 3, 993, 93, 96, 1032, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1830, 3, 996, 96, 99, 1035, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1833, 3, 999, 99, 102, 1038, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1836, 3, 1002, 102, 105, 1041, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 1839, 3, 1005, 105, 108, 1044, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1842, 0, 3, 1011, 1809, 114, 120, 1050, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1851, 0, 3, 1014, 1812, 120, 126, 1053, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1860, 0, 3, 1017, 1815, 126, 132, 1056, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1869, 0, 3, 1020, 1818, 132, 138, 1059, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1878, 0, 3, 1023, 1821, 138, 144, 1062, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1887, 0, 3, 1026, 1824, 144, 150, 1065, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1896, 0, 3, 1029, 1827, 150, 156, 1068, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1905, 0, 3, 1032, 1830, 156, 162, 1071, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1914, 0, 3, 1035, 1833, 162, 168, 1074, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1923, 0, 3, 1038, 1836, 168, 174, 1077, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 1932, 0, 3, 1041, 1839, 174, 180, 1080, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1941, 0, 3, 1803, 1806, 1047, 1842, 192, 201, 1083, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1956, 0, 3, 1806, 1809, 1050, 1851, 201, 210, 1086, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 1971, 0, 3, 1809, 1812, 1053, 1860, 210, 219, 1089, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_39(buffer, 1986, 0, 3, 1812, 1815, 1056, 1869, 219, 228, 1092, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 2001, 0, 3, 1815, 1818, 1059, 1878, 228, 237, 1101, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 2022, 0, 3, 1818, 1821, 1062, 1887, 237, 246, 1110, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 2043, 0, 3, 1821, 1824, 1065, 1896, 246, 255, 1119, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 2064, 0, 3, 1824, 1827, 1068, 1905, 255, 264, 1128, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_40(buffer, 2085, 0, 3, 1827, 1830, 1071, 1914, 264, 273, 1137, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 2106, 0, 3, 1830, 1833, 1074, 1923, 273, 282, 1146, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 2121, 0, 3, 1833, 1836, 1077, 1932, 282, 291, 1149, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 2136, 0, 3, 1842, 1851, 1086, 1971, 309, 321, 1155, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_29(buffer, 2157, 0, 3, 1851, 1860, 1089, 1986, 321, 333, 1158, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 2178, 0, 3, 1860, 1869, 1092, 2001, 333, 345, 1173, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 2211, 0, 3, 1869, 1878, 1101, 2022, 345, 357, 1188, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 2244, 0, 3, 1878, 1887, 1110, 2043, 357, 369, 1203, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 2277, 0, 3, 1887, 1896, 1119, 2064, 369, 381, 1218, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_44(buffer, 2310, 0, 3, 1896, 1905, 1128, 2085, 381, 393, 1233, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_45(buffer, 2343, 0, 3, 1905, 1914, 1137, 2106, 393, 405, 1248, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_28(buffer, 2370, 0, 3, 1914, 1923, 1146, 2121, 405, 417, 1257, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_28(buffer, 2391, 0, 3, 1941, 1956, 1152, 2136, 441, 456, 1260, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_33(buffer, 2418, 0, 3, 1956, 1971, 1155, 2157, 456, 471, 1263, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_19(buffer, 2445, 0, 3, 1971, 1986, 1158, 2178, 471, 486, 1281, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_37(buffer, 2487, 0, 3, 1986, 2001, 1173, 2211, 486, 501, 1299, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_38(buffer, 2529, 0, 3, 2001, 2022, 1188, 2244, 501, 516, 1317, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_38(buffer, 2571, 0, 3, 2022, 2043, 1203, 2277, 516, 531, 1335, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_38(buffer, 2613, 0, 3, 2043, 2064, 1218, 2310, 531, 546, 1353, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_35(buffer, 2655, 0, 3, 2064, 2085, 1233, 2343, 546, 561, 1371, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_36(buffer, 2697, 0, 3, 2085, 2106, 1248, 2370, 561, 576, 1389, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_19(buffer, 2730, 0, 3, 2136, 2157, 1263, 2445, 606, 624, 1419, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_20(buffer, 2781, 0, 3, 2157, 2178, 1281, 2487, 624, 642, 1440, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2832, 0, 3, 2178, 2211, 1299, 2529, 642, 660, 1461, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2883, 0, 3, 2211, 2244, 1317, 2571, 660, 678, 1482, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2934, 0, 3, 2244, 2277, 1335, 2613, 678, 696, 1503, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_21(buffer, 2985, 0, 3, 2277, 2310, 1353, 2655, 696, 714, 1524, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_18(buffer, 3036, 0, 3, 2310, 2343, 1371, 2697, 714, 732, 1545, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_2(buffer, 3087, 0, 3, 2391, 2418, 1398, 2730, 768, 786, 1566, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_5(buffer, 3141, 0, 3, 2418, 2445, 1419, 2781, 786, 804, 1590, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 3195, 0, 3, 2445, 2487, 1440, 2832, 804, 822, 1614, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 3249, 0, 3, 2487, 2529, 1461, 2883, 822, 840, 1638, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 3303, 0, 3, 2529, 2571, 1482, 2934, 840, 858, 1662, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 3357, 0, 3, 2571, 2613, 1503, 2985, 858, 876, 1686, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_8(buffer, 3411, 0, 3, 2613, 2655, 1524, 3036, 876, 894, 1710, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3465, 3, 930, 933, 1737, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3468, 3, 933, 936, 1740, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3471, 3, 936, 939, 1743, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3474, 3, 939, 942, 1746, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3477, 3, 942, 945, 1749, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3480, 3, 945, 948, 1752, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3483, 3, 948, 951, 1755, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3486, 3, 951, 954, 1758, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3489, 3, 954, 957, 1761, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3492, 3, 957, 960, 1764, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 3495, 3, 960, 963, 1767, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3498, 0, 1734, 3465, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3501, 0, 1737, 3468, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3504, 0, 1740, 3471, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3507, 0, 1743, 3474, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3510, 0, 1746, 3477, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3513, 0, 1749, 3480, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3516, 0, 1752, 3483, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3519, 0, 1755, 3486, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3522, 0, 1758, 3489, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3525, 0, 1761, 3492, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 3528, 0, 1764, 3495, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 3531, 3, 1770, 1008, 1011, 1809, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3534, 3, 1773, 1011, 1014, 1812, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3540, 3, 1776, 1014, 1017, 1815, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3546, 3, 1779, 1017, 1020, 1818, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3552, 3, 1782, 1020, 1023, 1821, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3558, 3, 1785, 1023, 1026, 1824, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3564, 3, 1788, 1026, 1029, 1827, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3570, 3, 1791, 1029, 1032, 1830, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3576, 3, 1794, 1032, 1035, 1833, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_24(buffer, 3582, 3, 1797, 1035, 1038, 1836, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 3588, 3, 1800, 1038, 1041, 1839, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_35(buffer, 3591, 0, 3, 1809, 3534, 1047, 1050, 1851, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_35(buffer, 3603, 0, 3, 1812, 3540, 1050, 1053, 1860, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3615, 0, 3, 1815, 3546, 1053, 1056, 1869, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3630, 0, 3, 1818, 3552, 1056, 1059, 1878, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3645, 0, 3, 1821, 3558, 1059, 1062, 1887, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3660, 0, 3, 1824, 3564, 1062, 1065, 1896, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3675, 0, 3, 1827, 3570, 1065, 1068, 1905, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3690, 0, 3, 1830, 3576, 1068, 1071, 1914, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_41(buffer, 3705, 0, 3, 1833, 3582, 1071, 1074, 1923, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_15(buffer, 3720, 0, 3, 1836, 3588, 1074, 1077, 1932, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_57(buffer, 3729, 0, 3, 3531, 3534, 1851, 3603, 1083, 1086, 1971, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_58(buffer, 3747, 0, 3, 3534, 3540, 1860, 3615, 1086, 1089, 1986, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_59(buffer, 3771, 0, 3, 3540, 3546, 1869, 3630, 1089, 1092, 2001, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 3795, 0, 3, 3546, 3552, 1878, 3645, 1092, 1101, 2022, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 3825, 0, 3, 3552, 3558, 1887, 3660, 1101, 1110, 2043, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 3855, 0, 3, 3558, 3564, 1896, 3675, 1110, 1119, 2064, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_60(buffer, 3885, 0, 3, 3564, 3570, 1905, 3690, 1119, 1128, 2085, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_61(buffer, 3915, 0, 3, 3570, 3576, 1914, 3705, 1128, 1137, 2106, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_62(buffer, 3939, 0, 3, 3576, 3582, 1923, 3720, 1137, 1146, 2121, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_52(buffer, 3957, 0, 3, 3591, 3603, 1971, 3747, 1152, 1155, 2157, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_53(buffer, 3987, 0, 3, 3603, 3615, 1986, 3771, 1155, 1158, 2178, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_55(buffer, 4017, 0, 3, 3615, 3630, 2001, 3795, 1158, 1173, 2211, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_55(buffer, 4068, 0, 3, 3630, 3645, 2022, 3825, 1173, 1188, 2244, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_55(buffer, 4119, 0, 3, 3645, 3660, 2043, 3855, 1188, 1203, 2277, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_55(buffer, 4170, 0, 3, 3660, 3675, 2064, 3885, 1203, 1218, 2310, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_56(buffer, 4221, 0, 3, 3675, 3690, 2085, 3915, 1218, 1233, 2343, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_57(buffer, 4266, 0, 3, 3690, 3705, 2106, 3939, 1233, 1248, 2370, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_34(buffer, 4293, 0, 3, 3729, 3747, 2157, 3987, 1260, 1263, 2445, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_39(buffer, 4329, 0, 3, 3747, 3771, 2178, 4017, 1263, 1281, 2487, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_40(buffer, 4410, 0, 3, 3771, 3795, 2211, 4068, 1281, 1299, 2529, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_41(buffer, 4491, 0, 3, 3795, 3825, 2244, 4119, 1299, 1317, 2571, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_41(buffer, 4572, 0, 3, 3825, 3855, 2277, 4170, 1317, 1335, 2613, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_37(buffer, 4653, 0, 3, 3855, 3885, 2310, 4221, 1335, 1353, 2655, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_38(buffer, 4734, 0, 3, 3885, 3915, 2343, 4266, 1353, 1371, 2697, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_16(buffer, 4785, 0, 3, 3957, 3987, 2445, 4329, 1398, 1419, 2781, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_21(buffer, 4881, 0, 3, 3987, 4017, 2487, 4410, 1419, 1440, 2832, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_22(buffer, 4980, 0, 3, 4017, 4068, 2529, 4491, 1440, 1461, 2883, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_22(buffer, 5079, 0, 3, 4068, 4119, 2571, 4572, 1461, 1482, 2934, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_22(buffer, 5178, 0, 3, 4119, 4170, 2613, 4653, 1482, 1503, 2985, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_19(buffer, 5277, 0, 3, 4170, 4221, 2655, 4734, 1503, 1524, 3036, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_5(buffer, 5373, 0, 3, 4293, 4329, 2781, 4881, 1566, 1590, 3195, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_8(buffer, 5481, 0, 3, 4329, 4410, 2832, 4980, 1590, 1614, 3249, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_8(buffer, 5589, 0, 3, 4410, 4491, 2883, 5079, 1614, 1638, 3303, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_8(buffer, 5697, 0, 3, 4491, 4572, 2934, 5178, 1638, 1662, 3357, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_9(buffer, 5805, 0, 3, 4572, 4653, 2985, 5277, 1662, 1686, 3411, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5913, 3, 1734, 1737, 3468, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5916, 3, 1737, 1740, 3471, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5919, 3, 1740, 1743, 3474, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5922, 3, 1743, 1746, 3477, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5925, 3, 1746, 1749, 3480, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5928, 3, 1749, 1752, 3483, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5931, 3, 1752, 1755, 3486, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5934, 3, 1755, 1758, 3489, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5937, 3, 1758, 1761, 3492, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 5940, 3, 1761, 1764, 3495, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5943, 0, 3468, 5916, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5946, 0, 3471, 5919, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5949, 0, 3474, 5922, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5952, 0, 3477, 5925, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5955, 0, 3480, 5928, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5958, 0, 3483, 5931, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5961, 0, 3486, 5934, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5964, 0, 3489, 5937, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 5967, 0, 3492, 5940, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 5970, 3, 3498, 1803, 1806, 3531, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_27(buffer, 5973, 3, 3501, 1806, 1809, 3534, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5976, 3, 3504, 1809, 1812, 3540, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5985, 3, 3507, 1812, 1815, 3546, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 5994, 3, 3510, 1815, 1818, 3552, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 6003, 3, 3513, 1818, 1821, 3558, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 6012, 3, 3516, 1821, 1824, 3564, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 6021, 3, 3519, 1824, 1827, 3570, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 6030, 3, 3522, 1827, 1830, 3576, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_31(buffer, 6039, 3, 3525, 1830, 1833, 3582, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 6048, 3, 3528, 1833, 1836, 3588, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_34(buffer, 6051, 0, 3, 3534, 5976, 1842, 1851, 3603, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_55(buffer, 6060, 0, 3, 3540, 5985, 1851, 1860, 3615, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 6075, 0, 3, 3546, 5994, 1860, 1869, 3630, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 6096, 0, 3, 3552, 6003, 1869, 1878, 3645, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 6117, 0, 3, 3558, 6012, 1878, 1887, 3660, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 6138, 0, 3, 3564, 6021, 1887, 1896, 3675, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 6159, 0, 3, 3570, 6030, 1896, 1905, 3690, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_49(buffer, 6180, 0, 3, 3576, 6039, 1905, 1914, 3705, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_54(buffer, 6201, 0, 3, 3582, 6048, 1914, 1923, 3720, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_29(buffer, 6210, 0, 3, 5970, 5973, 3591, 6051, 1941, 1956, 3729, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_62(buffer, 6231, 0, 3, 5973, 5976, 3603, 6060, 1956, 1971, 3747, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_63(buffer, 6252, 0, 3, 5976, 5985, 3615, 6075, 1971, 1986, 3771, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_64(buffer, 6285, 0, 3, 5985, 5994, 3630, 6096, 1986, 2001, 3795, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_65(buffer, 6318, 0, 3, 5994, 6003, 3645, 6117, 2001, 2022, 3825, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_65(buffer, 6357, 0, 3, 6003, 6012, 3660, 6138, 2022, 2043, 3855, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_65(buffer, 6396, 0, 3, 6012, 6021, 3675, 6159, 2043, 2064, 3885, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_66(buffer, 6435, 0, 3, 6021, 6030, 3690, 6180, 2064, 2085, 3915, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_67(buffer, 6468, 0, 3, 6030, 6039, 3705, 6201, 2085, 2106, 3939, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_48(buffer, 6489, 0, 3, 6051, 6060, 3747, 6252, 2136, 2157, 3987, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_54(buffer, 6523, 0, 3, 6060, 6075, 3771, 6285, 2157, 2178, 4017, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_56(buffer, 6563, 0, 3, 6075, 6096, 3795, 6318, 2178, 2211, 4068, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_56(buffer, 6633, 0, 3, 6096, 6117, 3825, 6357, 2211, 2244, 4119, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_56(buffer, 6703, 0, 3, 6117, 6138, 3855, 6396, 2244, 2277, 4170, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_57(buffer, 6773, 0, 3, 6138, 6159, 3885, 6435, 2277, 2310, 4221, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_58(buffer, 6837, 0, 3, 6159, 6180, 3915, 6468, 2310, 2343, 4266, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_29(buffer, 6871, 0, 3, 6210, 6231, 3957, 6489, 2391, 2418, 4293, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_34(buffer, 6919, 0, 3, 6231, 6252, 3987, 6523, 2418, 2445, 4329, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_39(buffer, 6967, 0, 3, 6252, 6285, 4017, 6563, 2445, 2487, 4410, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_40(buffer, 7090, 0, 3, 6285, 6318, 4068, 6633, 2487, 2529, 4491, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_41(buffer, 7213, 0, 3, 6318, 6357, 4119, 6703, 2529, 2571, 4572, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_37(buffer, 7336, 0, 3, 6357, 6396, 4170, 6773, 2571, 2613, 4653, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_38(buffer, 7459, 0, 3, 6396, 6435, 4221, 6837, 2613, 2655, 4734, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_20(buffer, 7531, 0, 3, 6489, 6523, 4329, 6967, 2730, 2781, 4881, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_21(buffer, 7678, 0, 3, 6523, 6563, 4410, 7090, 2781, 2832, 4980, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_22(buffer, 7831, 0, 3, 6563, 6633, 4491, 7213, 2832, 2883, 5079, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_22(buffer, 7984, 0, 3, 6633, 6703, 4572, 7336, 2883, 2934, 5178, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_19(buffer, 8137, 0, 3, 6703, 6773, 4653, 7459, 2934, 2985, 5277, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_2(buffer, 8284, 0, 3, 6871, 6919, 4785, 7531, 3087, 3141, 5373, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_5(buffer, 8452, 0, 3, 6919, 6967, 4881, 7678, 3141, 3195, 5481, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_8(buffer, 8620, 0, 3, 6967, 7090, 4980, 7831, 3195, 3249, 5589, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_8(buffer, 8788, 0, 3, 7090, 7213, 5079, 7984, 3249, 3303, 5697, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_9(buffer, 8956, 0, 3, 7213, 7336, 5178, 8137, 3303, 3357, 5805, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9124, 3, 3465, 3468, 5916, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9127, 3, 3468, 3471, 5919, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9130, 3, 3471, 3474, 5922, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9133, 3, 3474, 3477, 5925, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9136, 3, 3477, 3480, 5928, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9139, 3, 3480, 3483, 5931, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9142, 3, 3483, 3486, 5934, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9145, 3, 3486, 3489, 5937, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_12(buffer, 9148, 3, 3489, 3492, 5940, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9151, 0, 5913, 9124, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9154, 0, 5916, 9127, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9157, 0, 5919, 9130, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9160, 0, 5922, 9133, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9163, 0, 5925, 9136, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9166, 0, 5928, 9139, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9169, 0, 5931, 9142, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9172, 0, 5934, 9145, ncols, p);

            compute_prim_ph_electron_repulsion_17(buffer, 9175, 0, 5937, 9148, ncols, p);

            compute_prim_dh_electron_repulsion_32(buffer, 9178, 3, 5943, 3531, 3534, 5976, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 9181, 3, 5946, 3534, 3540, 5985, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 9193, 3, 5949, 3540, 3546, 5994, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 9205, 3, 5952, 3546, 3552, 6003, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 9217, 3, 5955, 3552, 3558, 6012, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 9229, 3, 5958, 3558, 3564, 6021, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 9241, 3, 5961, 3564, 3570, 6030, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_36(buffer, 9253, 3, 5964, 3570, 3576, 6039, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_40(buffer, 9265, 3, 5967, 3576, 3582, 6048, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_38(buffer, 9268, 0, 3, 5976, 9181, 3591, 3603, 6060, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_59(buffer, 9286, 0, 3, 5985, 9193, 3603, 3615, 6075, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 9304, 0, 3, 5994, 9205, 3615, 3630, 6096, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 9331, 0, 3, 6003, 9217, 3630, 3645, 6117, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 9358, 0, 3, 6012, 9229, 3645, 3660, 6138, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 9385, 0, 3, 6021, 9241, 3660, 3675, 6159, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_51(buffer, 9412, 0, 3, 6030, 9253, 3675, 3690, 6180, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_58(buffer, 9439, 0, 3, 6039, 9265, 3690, 3705, 6201, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_57(buffer, 9448, 0, 3, 9178, 9181, 6060, 9286, 3729, 3747, 6252, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_58(buffer, 9472, 0, 3, 9181, 9193, 6075, 9304, 3747, 3771, 6285, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_59(buffer, 9514, 0, 3, 9193, 9205, 6096, 9331, 3771, 3795, 6318, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_60(buffer, 9556, 0, 3, 9205, 9217, 6117, 9358, 3795, 3825, 6357, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_60(buffer, 9604, 0, 3, 9217, 9229, 6138, 9385, 3825, 3855, 6396, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_61(buffer, 9652, 0, 3, 9229, 9241, 6159, 9412, 3855, 3885, 6435, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_62(buffer, 9694, 0, 3, 9241, 9253, 6180, 9439, 3885, 3915, 6468, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_45(buffer, 9718, 0, 3, 9268, 9286, 6252, 9472, 3957, 3987, 6523, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_46(buffer, 9768, 0, 3, 9286, 9304, 6285, 9514, 3987, 4017, 6563, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_47(buffer, 9818, 0, 3, 9304, 9331, 6318, 9556, 4017, 4068, 6633, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_48(buffer, 9913, 0, 3, 9331, 9358, 6357, 9604, 4068, 4119, 6703, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_49(buffer, 10004, 0, 3, 9358, 9385, 6396, 9652, 4119, 4170, 6773, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_50(buffer, 10087, 0, 3, 9385, 9412, 6435, 9694, 4170, 4221, 6837, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_29(buffer, 10128, 0, 3, 9448, 9472, 6523, 9768, 4293, 4329, 6967, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_30(buffer, 10188, 0, 3, 9472, 9514, 6563, 9818, 4329, 4410, 7090, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_31(buffer, 10359, 0, 3, 9514, 9556, 6633, 9913, 4410, 4491, 7213, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_32(buffer, 10539, 0, 3, 9556, 9604, 6703, 10004, 4491, 4572, 7336, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_33(buffer, 10710, 0, 3, 9604, 9652, 6773, 10087, 4572, 4653, 7459, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_14(buffer, 10806, 0, 3, 9718, 9768, 6967, 10188, 4785, 4881, 7678, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_15(buffer, 11031, 0, 3, 9768, 9818, 7090, 10359, 4881, 4980, 7831, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_16(buffer, 11259, 0, 3, 9818, 9913, 7213, 10539, 4980, 5079, 7984, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_17(buffer, 11490, 0, 3, 9913, 10004, 7336, 10710, 5079, 5178, 8137, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_5(buffer, 11703, 0, 3, 10128, 10188, 7678, 11031, 5373, 5481, 8620, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_6(buffer, 11955, 0, 3, 10188, 10359, 7831, 11259, 5481, 5589, 8788, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_7(buffer, 12207, 0, 3, 10359, 10539, 7984, 11490, 5589, 5697, 8956, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 12459, 3, 5913, 5916, 9127, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 12462, 3, 5916, 5919, 9130, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 12465, 3, 5919, 5922, 9133, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 12468, 3, 5922, 5925, 9136, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 12471, 3, 5925, 5928, 9139, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 12474, 3, 5928, 5931, 9142, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 12477, 3, 5931, 5934, 9145, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_11(buffer, 12480, 3, 5934, 5937, 9148, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_14(buffer, 12483, 0, 9127, 12462, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 12486, 0, 9130, 12465, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 12489, 0, 9133, 12468, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 12492, 0, 9136, 12471, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 12495, 0, 9139, 12474, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 12498, 0, 9142, 12477, ncols, p);

            compute_prim_pi_electron_repulsion_14(buffer, 12501, 0, 9145, 12480, ncols, p);

            compute_prim_di_electron_repulsion_9(buffer, 12504, 3, 9151, 5970, 5973, 9178, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_29(buffer, 12507, 3, 9154, 5973, 5976, 9181, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 12510, 3, 9157, 5976, 5985, 9193, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 12525, 3, 9160, 5985, 5994, 9205, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 12540, 3, 9163, 5994, 6003, 9217, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 12555, 3, 9166, 6003, 6012, 9229, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 12570, 3, 9169, 6012, 6021, 9241, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_33(buffer, 12585, 3, 9172, 6021, 6030, 9253, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_37(buffer, 12600, 3, 9175, 6030, 6039, 9265, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_27(buffer, 12603, 0, 3, 9181, 12510, 6051, 6060, 9286, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_48(buffer, 12612, 0, 3, 9193, 12525, 6060, 6075, 9304, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_40(buffer, 12633, 0, 3, 9205, 12540, 6075, 6096, 9331, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_40(buffer, 12666, 0, 3, 9217, 12555, 6096, 6117, 9358, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_40(buffer, 12699, 0, 3, 9229, 12570, 6117, 6138, 9385, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_40(buffer, 12732, 0, 3, 9241, 12585, 6138, 6159, 9412, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_47(buffer, 12765, 0, 3, 9253, 12600, 6159, 6180, 9439, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_17(buffer, 12774, 0, 3, 12504, 12507, 9268, 12603, 6210, 6231, 9448, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_43(buffer, 12801, 0, 3, 12507, 12510, 9286, 12612, 6231, 6252, 9472, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_44(buffer, 12828, 0, 3, 12510, 12525, 9304, 12633, 6252, 6285, 9514, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_45(buffer, 12879, 0, 3, 12525, 12540, 9331, 12666, 6285, 6318, 9556, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_46(buffer, 12930, 0, 3, 12540, 12555, 9358, 12699, 6318, 6357, 9604, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_47(buffer, 12987, 0, 3, 12555, 12570, 9385, 12732, 6357, 6396, 9652, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_48(buffer, 13038, 0, 3, 12570, 12585, 9412, 12765, 6396, 6435, 9694, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_29(buffer, 13065, 0, 3, 12603, 12612, 9472, 12828, 6489, 6523, 9768, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_30(buffer, 13113, 0, 3, 12612, 12633, 9514, 12879, 6523, 6563, 9818, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_31(buffer, 13176, 0, 3, 12633, 12666, 9556, 12930, 6563, 6633, 9913, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_32(buffer, 13295, 0, 3, 12666, 12699, 9604, 12987, 6633, 6703, 10004, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_33(buffer, 13400, 0, 3, 12699, 12732, 9652, 13038, 6703, 6773, 10087, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_17(buffer, 13451, 0, 3, 12774, 12801, 9718, 13065, 6871, 6919, 10128, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_18(buffer, 13523, 0, 3, 12801, 12828, 9768, 13113, 6919, 6967, 10188, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_19(buffer, 13595, 0, 3, 12828, 12879, 9818, 13176, 6967, 7090, 10359, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_20(buffer, 13883, 0, 3, 12879, 12930, 9913, 13295, 7090, 7213, 10539, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_21(buffer, 14126, 0, 3, 12930, 12987, 10004, 13400, 7213, 7336, 10710, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_8(buffer, 14258, 0, 3, 13065, 13113, 10188, 13595, 7531, 7678, 11031, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_9(buffer, 14552, 0, 3, 13113, 13176, 10359, 13883, 7678, 7831, 11259, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_10(buffer, 14936, 0, 3, 13176, 13295, 10539, 14126, 7831, 7984, 11490, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_2(buffer, 15254, 0, 3, 13451, 13523, 10806, 14258, 8284, 8452, 11703, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_3(buffer, 15614, 0, 3, 13523, 13595, 11031, 14552, 8452, 8620, 11955, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_4(buffer, 15974, 0, 3, 13595, 13883, 11259, 14936, 8620, 8788, 12207, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 16394, 3, 9124, 9127, 12462, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 16397, 3, 9127, 9130, 12465, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 16400, 3, 9130, 9133, 12468, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 16403, 3, 9133, 9136, 12471, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 16406, 3, 9136, 9139, 12474, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 16409, 3, 9139, 9142, 12477, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_9(buffer, 16412, 3, 9142, 9145, 12480, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_8(buffer, 16415, 0, 12459, 16394, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 16418, 0, 12462, 16397, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 16421, 0, 12465, 16400, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 16424, 0, 12468, 16403, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 16427, 0, 12471, 16406, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 16430, 0, 12474, 16409, ncols, p);

            compute_prim_pk_electron_repulsion_8(buffer, 16433, 0, 12477, 16412, ncols, p);

            compute_prim_dk_electron_repulsion_20(buffer, 16436, 3, 12483, 9178, 9181, 12510, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 16439, 3, 12486, 9181, 9193, 12525, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 16454, 3, 12489, 9193, 9205, 12540, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 16469, 3, 12492, 9205, 9217, 12555, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 16484, 3, 12495, 9217, 9229, 12570, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_21(buffer, 16499, 3, 12498, 9229, 9241, 12585, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_24(buffer, 16514, 3, 12501, 9241, 9253, 12600, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_14(buffer, 16517, 0, 3, 12510, 16439, 9268, 9286, 12612, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_29(buffer, 16541, 0, 3, 12525, 16454, 9286, 9304, 12633, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_23(buffer, 16565, 0, 3, 12540, 16469, 9304, 9331, 12666, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_23(buffer, 16601, 0, 3, 12555, 16484, 9331, 9358, 12699, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_23(buffer, 16637, 0, 3, 12570, 16499, 9358, 9385, 12732, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_30(buffer, 16673, 0, 3, 12585, 16514, 9385, 9412, 12765, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_22(buffer, 16682, 0, 3, 16436, 16439, 12612, 16541, 9448, 9472, 12828, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_23(buffer, 16712, 0, 3, 16439, 16454, 12633, 16565, 9472, 9514, 12879, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_24(buffer, 16769, 0, 3, 16454, 16469, 12666, 16601, 9514, 9556, 12930, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_25(buffer, 16826, 0, 3, 16469, 16484, 12699, 16637, 9556, 9604, 12987, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_26(buffer, 16883, 0, 3, 16484, 16499, 12732, 16673, 9604, 9652, 13038, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_14(buffer, 16913, 0, 3, 16517, 16541, 12828, 16712, 9718, 9768, 13113, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_15(buffer, 16994, 0, 3, 16541, 16565, 12879, 16769, 9768, 9818, 13176, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_16(buffer, 17075, 0, 3, 16565, 16601, 12930, 16826, 9818, 9913, 13295, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_17(buffer, 17198, 0, 3, 16601, 16637, 12987, 16883, 9913, 10004, 13400, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_8(buffer, 17264, 0, 3, 16682, 16712, 13113, 16994, 10128, 10188, 13595, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_9(buffer, 17381, 0, 3, 16712, 16769, 13176, 17075, 10188, 10359, 13883, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_10(buffer, 17716, 0, 3, 16769, 16826, 13295, 17198, 10359, 10539, 14126, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_3(buffer, 17902, 0, 3, 16913, 16994, 13595, 17381, 10806, 11031, 14552, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_4(buffer, 18725, 0, 3, 16994, 17075, 13883, 17716, 11031, 11259, 14936, ncols, alpha, beta, p);

            compute_prim_lk_electron_repulsion_1(buffer, 19240, 0, 3, 17264, 17381, 14552, 18725, 11703, 11955, 15974, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_2(buffer, 20227, 3, 16415, 12504, 12507, 16436, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_7(buffer, 20230, 3, 16418, 12507, 12510, 16439, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_8(buffer, 20233, 3, 16421, 12510, 12525, 16454, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_8(buffer, 20236, 3, 16424, 12525, 12540, 16469, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_8(buffer, 20239, 3, 16427, 12540, 12555, 16484, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_8(buffer, 20242, 3, 16430, 12555, 12570, 16499, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_11(buffer, 20245, 3, 16433, 12570, 12585, 16514, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_5(buffer, 20248, 0, 3, 16439, 20233, 12603, 12612, 16541, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_12(buffer, 20257, 0, 3, 16454, 20236, 12612, 12633, 16565, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_9(buffer, 20266, 0, 3, 16469, 20239, 12633, 12666, 16601, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_9(buffer, 20275, 0, 3, 16484, 20242, 12666, 12699, 16637, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_13(buffer, 20284, 0, 3, 16499, 20245, 12699, 12732, 16673, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_2(buffer, 20293, 0, 3, 20227, 20230, 16517, 20248, 12774, 12801, 16682, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_8(buffer, 20323, 0, 3, 20230, 20233, 16541, 20257, 12801, 12828, 16712, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_9(buffer, 20353, 0, 3, 20233, 20236, 16565, 20266, 12828, 12879, 16769, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_10(buffer, 20383, 0, 3, 20236, 20239, 16601, 20275, 12879, 12930, 16826, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_11(buffer, 20413, 0, 3, 20239, 20242, 16637, 20284, 12930, 12987, 16883, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_5(buffer, 20443, 0, 3, 20248, 20257, 16712, 20353, 13065, 13113, 16994, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_6(buffer, 20512, 0, 3, 20257, 20266, 16769, 20383, 13113, 13176, 17075, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_7(buffer, 20581, 0, 3, 20266, 20275, 16826, 20413, 13176, 13295, 17198, ncols, alpha, beta, p);

            compute_prim_il_electron_repulsion_2(buffer, 20650, 0, 3, 20293, 20323, 16913, 20443, 13451, 13523, 17264, ncols, alpha, beta, p);

            compute_prim_il_electron_repulsion_3(buffer, 20776, 0, 3, 20323, 20353, 16994, 20512, 13523, 13595, 17381, ncols, alpha, beta, p);

            compute_prim_il_electron_repulsion_4(buffer, 20902, 0, 3, 20353, 20383, 17075, 20581, 13595, 13883, 17716, ncols, alpha, beta, p);

            compute_prim_kl_electron_repulsion_1(buffer, 21097, 0, 3, 20443, 20512, 17381, 20902, 14258, 14552, 18725, ncols, alpha, beta, p);

            compute_prim_ll_electron_repulsion_0(buffer, 21682, 0, 3, 20650, 20776, 17902, 21097, 15254, 15614, 19240, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 23707, 21682, 2025, ncols);
        }
    }

    simdtrf::transform_ll_tri(values, nvalues, buffer, 23707, nmax);
}

}  // namespace simdt2ceri
