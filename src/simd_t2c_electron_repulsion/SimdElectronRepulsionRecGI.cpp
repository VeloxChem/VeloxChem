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


#include "SimdElectronRepulsionRecGI.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformGI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gi_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(2920, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 45, 0, 7, 8, 18, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 48, 0, 8, 9, 21, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 9, 10, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 54, 0, 10, 11, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 57, 0, 11, 12, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 60, 0, 12, 13, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 63, 0, 13, 14, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 66, 0, 14, 15, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 15, 16, 42, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 72, 0, 18, 21, 51, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 78, 0, 21, 24, 54, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 84, 0, 24, 27, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 90, 0, 27, 30, 60, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 96, 0, 30, 33, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 102, 0, 33, 36, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 108, 0, 36, 39, 69, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 114, 0, 45, 48, 72, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 120, 0, 48, 51, 78, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 126, 0, 51, 54, 84, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 132, 0, 54, 57, 90, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 138, 0, 57, 60, 96, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 144, 0, 60, 63, 102, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 150, 0, 63, 66, 108, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 156, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 159, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 162, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 165, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 168, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 171, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 174, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 177, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 180, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 183, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 186, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 189, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 192, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 195, 3, 15, 39, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 198, 3, 21, 51, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 201, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 210, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 219, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 228, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 237, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 246, 3, 39, 69, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 255, 3, 51, 78, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 264, 3, 54, 84, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 273, 3, 57, 90, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 282, 3, 60, 96, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 291, 3, 63, 102, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 300, 3, 66, 108, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 309, 3, 78, 126, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 321, 3, 84, 132, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 333, 3, 90, 138, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 345, 3, 96, 144, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 357, 3, 102, 150, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 369, 3, 9, 10, 159, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 372, 3, 10, 11, 162, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 375, 3, 11, 12, 165, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 378, 3, 12, 13, 168, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 381, 3, 13, 14, 171, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 384, 3, 14, 15, 174, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 387, 0, 159, 372, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 390, 0, 162, 375, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 393, 0, 165, 378, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 396, 0, 168, 381, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 399, 0, 171, 384, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 402, 3, 177, 45, 48, 198, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 405, 3, 180, 48, 51, 201, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 408, 3, 183, 51, 54, 210, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 417, 3, 186, 54, 57, 219, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 426, 3, 189, 57, 60, 228, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 435, 3, 192, 60, 63, 237, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 444, 3, 195, 63, 66, 246, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 453, 0, 3, 201, 408, 72, 78, 264, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 468, 0, 3, 210, 417, 78, 84, 273, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 483, 0, 3, 219, 426, 84, 90, 282, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 498, 0, 3, 228, 435, 90, 96, 291, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 513, 0, 3, 237, 444, 96, 102, 300, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_2(buffer, 528, 0, 3, 402, 405, 255, 453, 114, 120, 309, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_5(buffer, 546, 0, 3, 405, 408, 264, 468, 120, 126, 321, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 564, 0, 3, 408, 417, 273, 483, 126, 132, 333, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 582, 0, 3, 417, 426, 282, 498, 132, 138, 345, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 600, 0, 3, 426, 435, 291, 513, 138, 144, 357, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 618, 3, 156, 159, 372, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 621, 3, 159, 162, 375, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_3(buffer, 624, 3, 162, 165, 378, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 633, 3, 165, 168, 381, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 639, 3, 168, 171, 384, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 645, 0, 369, 618, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 648, 0, 372, 621, ncols, p);

            compute_prim_pf_electron_repulsion_9(buffer, 651, 0, 375, 624, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 654, 0, 378, 633, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 657, 0, 381, 639, ncols, p);

            compute_prim_df_electron_repulsion_14(buffer, 660, 3, 387, 198, 201, 408, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 663, 3, 390, 201, 210, 417, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_17(buffer, 681, 3, 393, 210, 219, 426, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 702, 3, 396, 219, 228, 435, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 720, 3, 399, 228, 237, 444, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_10(buffer, 738, 0, 3, 408, 663, 255, 264, 468, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_11(buffer, 762, 0, 3, 417, 681, 264, 273, 483, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_12(buffer, 789, 0, 3, 426, 702, 273, 282, 498, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_13(buffer, 819, 0, 3, 435, 720, 282, 291, 513, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_5(buffer, 843, 0, 3, 660, 663, 468, 762, 309, 321, 564, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_6(buffer, 879, 0, 3, 663, 681, 483, 789, 321, 333, 582, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_7(buffer, 915, 0, 3, 681, 702, 498, 819, 333, 345, 600, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 951, 3, 369, 372, 621, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_8(buffer, 954, 3, 372, 375, 624, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_9(buffer, 957, 3, 375, 378, 633, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 966, 3, 378, 381, 639, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_14(buffer, 975, 0, 621, 954, ncols, p);

            compute_prim_pg_electron_repulsion_15(buffer, 978, 0, 3, 624, 957, 654, ncols, p);

            compute_prim_pg_electron_repulsion_13(buffer, 990, 0, 633, 966, ncols, p);

            compute_prim_dg_electron_repulsion_9(buffer, 993, 3, 645, 402, 405, 660, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_14(buffer, 996, 3, 648, 405, 408, 663, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_15(buffer, 999, 0, 3, 651, 978, 408, 417, 681, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_16(buffer, 1032, 0, 3, 654, 990, 417, 426, 702, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 1065, 3, 657, 426, 435, 720, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_7(buffer, 1092, 0, 3, 663, 999, 453, 468, 762, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_8(buffer, 1125, 0, 3, 681, 1032, 468, 483, 789, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_9(buffer, 1184, 0, 3, 702, 1065, 483, 498, 819, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_2(buffer, 1223, 0, 3, 993, 996, 738, 1092, 528, 546, 843, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_3(buffer, 1277, 0, 3, 996, 999, 762, 1125, 546, 564, 879, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_4(buffer, 1331, 0, 3, 999, 1032, 789, 1184, 564, 582, 915, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_8(buffer, 1403, 3, 618, 621, 954, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_9(buffer, 1406, 3, 621, 624, 957, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_10(buffer, 1409, 3, 624, 633, 966, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_8(buffer, 1418, 0, 951, 1403, ncols, p);

            compute_prim_ph_electron_repulsion_8(buffer, 1421, 0, 954, 1406, ncols, p);

            compute_prim_ph_electron_repulsion_5(buffer, 1424, 0, 957, 1409, ncols, p);

            compute_prim_dh_electron_repulsion_8(buffer, 1433, 3, 975, 660, 663, 999, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_9(buffer, 1436, 0, 3, 978, 1424, 663, 681, 1032, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_10(buffer, 1493, 3, 990, 681, 702, 1065, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_3(buffer, 1532, 0, 3, 999, 1436, 738, 762, 1125, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_4(buffer, 1663, 0, 3, 1032, 1493, 762, 789, 1184, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_1(buffer, 1741, 0, 3, 1433, 1436, 1125, 1663, 843, 879, 1331, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_2(buffer, 1936, 3, 1418, 993, 996, 1433, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_3(buffer, 1939, 3, 1421, 996, 999, 1436, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_4(buffer, 1942, 3, 1424, 999, 1032, 1493, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_1(buffer, 1981, 0, 3, 1436, 1942, 1092, 1125, 1663, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 2080, 0, 3, 1936, 1939, 1532, 1981, 1223, 1277, 1741, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2500, 2080, 420, ncols);
        }
    }

    simdtrf::transform_gi(values, nvalues, buffer, 2500, nmax);
}

}  // namespace simdt2ceri
