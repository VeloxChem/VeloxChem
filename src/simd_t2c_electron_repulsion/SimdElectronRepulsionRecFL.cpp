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


#include "SimdElectronRepulsionRecFL.hpp"

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
#include "SimdTransformFL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fl_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(3434, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 7, 8, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 54, 0, 8, 9, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 57, 0, 9, 10, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 60, 0, 10, 11, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 63, 0, 11, 12, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 66, 0, 12, 13, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 13, 14, 42, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 72, 0, 14, 15, 45, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 75, 0, 15, 16, 48, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 78, 0, 18, 21, 51, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 81, 0, 21, 24, 54, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 84, 0, 24, 27, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 87, 0, 27, 30, 60, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 90, 0, 30, 33, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 93, 0, 33, 36, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 96, 0, 36, 39, 69, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 99, 0, 39, 42, 72, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 102, 0, 42, 45, 75, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 105, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 108, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 111, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 114, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 117, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 120, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 123, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 126, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 129, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 132, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 135, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 138, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 141, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 144, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 147, 3, 15, 45, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 150, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 159, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 168, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 177, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 186, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 195, 3, 39, 69, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 204, 3, 42, 72, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 213, 3, 45, 75, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 222, 3, 54, 84, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 231, 3, 57, 87, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 240, 3, 60, 90, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 249, 3, 63, 93, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 258, 3, 66, 96, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 267, 3, 69, 99, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 276, 3, 72, 102, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 285, 3, 9, 10, 108, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 288, 3, 10, 11, 111, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 291, 3, 11, 12, 114, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 294, 3, 12, 13, 117, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 297, 3, 13, 14, 120, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 300, 3, 14, 15, 123, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 303, 3, 15, 16, 126, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 306, 0, 105, 285, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 309, 0, 108, 288, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 312, 0, 111, 291, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 315, 0, 114, 294, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 318, 0, 117, 297, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 321, 0, 120, 300, ncols, p);

            compute_prim_dd_electron_repulsion_2(buffer, 324, 3, 129, 51, 54, 159, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 333, 3, 132, 54, 57, 168, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 342, 3, 135, 57, 60, 177, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 351, 3, 138, 60, 63, 186, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 360, 3, 141, 63, 66, 195, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 369, 3, 144, 66, 69, 204, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 378, 3, 147, 69, 72, 213, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 387, 3, 150, 78, 81, 222, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 396, 3, 159, 81, 84, 231, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 405, 3, 168, 84, 87, 240, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 414, 3, 177, 87, 90, 249, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 423, 3, 186, 90, 93, 258, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 432, 3, 195, 93, 96, 267, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 441, 3, 204, 96, 99, 276, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 450, 3, 105, 108, 288, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 456, 3, 108, 111, 291, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 462, 3, 111, 114, 294, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 468, 3, 114, 117, 297, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 474, 3, 117, 120, 300, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 480, 3, 120, 123, 303, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_8(buffer, 486, 0, 285, 450, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 489, 0, 288, 456, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 492, 0, 291, 462, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 495, 0, 294, 468, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 498, 0, 297, 474, ncols, p);

            compute_prim_df_electron_repulsion_2(buffer, 501, 3, 306, 150, 159, 333, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 519, 3, 309, 159, 168, 342, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 537, 3, 312, 168, 177, 351, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 555, 3, 315, 177, 186, 360, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 573, 3, 318, 186, 195, 369, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 591, 3, 321, 195, 204, 378, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 609, 3, 333, 222, 231, 405, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 627, 3, 342, 231, 240, 414, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 645, 3, 351, 240, 249, 423, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 663, 3, 360, 249, 258, 432, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 681, 3, 369, 258, 267, 441, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 699, 3, 285, 288, 456, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 708, 3, 288, 291, 462, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 719, 3, 291, 294, 468, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 728, 3, 294, 297, 474, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 737, 3, 297, 300, 480, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_8(buffer, 746, 0, 450, 699, ncols, p);

            compute_prim_pg_electron_repulsion_16(buffer, 749, 0, 456, 708, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 752, 0, 462, 719, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 755, 0, 468, 728, ncols, p);

            compute_prim_dg_electron_repulsion_2(buffer, 758, 3, 486, 324, 333, 519, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 785, 3, 489, 333, 342, 537, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 812, 3, 492, 342, 351, 555, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 839, 3, 495, 351, 360, 573, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 866, 3, 498, 360, 369, 591, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 893, 3, 501, 387, 396, 609, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 920, 3, 519, 396, 405, 627, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 947, 3, 537, 405, 414, 645, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 974, 3, 555, 414, 423, 663, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 1001, 3, 573, 423, 432, 681, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_11(buffer, 1028, 3, 450, 456, 708, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_4(buffer, 1045, 3, 456, 462, 719, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 1061, 3, 462, 468, 728, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 1074, 3, 468, 474, 737, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_15(buffer, 1087, 0, 699, 1028, ncols, p);

            compute_prim_ph_electron_repulsion_16(buffer, 1090, 0, 3, 708, 1045, 752, ncols, p);

            compute_prim_ph_electron_repulsion_8(buffer, 1098, 0, 719, 1061, ncols, p);

            compute_prim_dh_electron_repulsion_2(buffer, 1101, 3, 746, 501, 519, 785, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 1140, 3, 749, 519, 537, 812, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 1179, 3, 752, 537, 555, 839, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 1218, 3, 755, 555, 573, 866, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_2(buffer, 1257, 3, 785, 609, 627, 947, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_2(buffer, 1296, 3, 812, 627, 645, 974, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_2(buffer, 1335, 3, 839, 645, 663, 1001, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_9(buffer, 1374, 3, 699, 708, 1045, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_10(buffer, 1397, 3, 708, 719, 1061, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 1418, 3, 719, 728, 1074, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_11(buffer, 1436, 0, 3, 1028, 1374, 1090, ncols, p);

            compute_prim_pi_electron_repulsion_12(buffer, 1457, 0, 3, 1045, 1397, 1098, ncols, p);

            compute_prim_pi_electron_repulsion_13(buffer, 1475, 0, 1061, 1418, ncols, p);

            compute_prim_di_electron_repulsion_2(buffer, 1478, 3, 1087, 758, 785, 1140, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_7(buffer, 1532, 0, 3, 1090, 1457, 785, 812, 1179, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_8(buffer, 1598, 0, 3, 1098, 1475, 812, 839, 1218, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_2(buffer, 1655, 3, 1101, 893, 920, 1257, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_2(buffer, 1709, 3, 1140, 920, 947, 1296, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_3(buffer, 1763, 0, 3, 1179, 1598, 947, 974, 1335, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_8(buffer, 1832, 3, 1028, 1045, 1397, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_6(buffer, 1856, 3, 1045, 1061, 1418, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_6(buffer, 1880, 0, 3, 1374, 1832, 1457, ncols, p);

            compute_prim_pk_electron_repulsion_7(buffer, 1925, 0, 1397, 1856, ncols, p);

            compute_prim_dk_electron_repulsion_3(buffer, 1940, 0, 3, 1436, 1880, 1101, 1140, 1532, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_4(buffer, 2081, 0, 3, 1457, 1925, 1140, 1179, 1598, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_1(buffer, 2177, 0, 3, 1532, 2081, 1257, 1296, 1763, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_2(buffer, 2378, 3, 1374, 1397, 1856, ncols, alpha, beta, p);

            compute_prim_pl_electron_repulsion_2(buffer, 2402, 0, 1832, 2378, ncols, p);

            compute_prim_dl_electron_repulsion_1(buffer, 2426, 0, 3, 1880, 2402, 1478, 1532, 2081, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 2534, 0, 3, 1940, 2426, 1655, 1709, 2177, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2984, 2534, 450, ncols);
        }
    }

    simdtrf::transform_fl(values, nvalues, buffer, 2984, nmax);
}

}  // namespace simdt2ceri
