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


#include "SimdElectronRepulsionRecFK.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
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
#include "SimdTransformFK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fk_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(2465, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 16, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 44, 0, 7, 8, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 47, 0, 8, 9, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 50, 0, 9, 10, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 53, 0, 10, 11, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 56, 0, 11, 12, 32, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 59, 0, 12, 13, 35, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 62, 0, 13, 14, 38, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 65, 0, 14, 15, 41, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 68, 0, 17, 20, 47, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 71, 0, 20, 23, 50, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 74, 0, 23, 26, 53, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 77, 0, 26, 29, 56, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 80, 0, 29, 32, 59, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 83, 0, 32, 35, 62, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 86, 0, 35, 38, 65, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 89, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 92, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 95, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 98, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 101, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 104, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 107, 3, 16, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 110, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 113, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 116, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 119, 3, 12, 32, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 122, 3, 13, 35, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 125, 3, 14, 38, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 128, 3, 20, 47, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 137, 3, 23, 50, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 146, 3, 26, 53, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 155, 3, 29, 56, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 164, 3, 32, 59, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 173, 3, 35, 62, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 182, 3, 38, 65, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 191, 3, 44, 68, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 200, 3, 47, 71, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 209, 3, 50, 74, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 218, 3, 53, 77, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 227, 3, 56, 80, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 236, 3, 59, 83, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 245, 3, 62, 86, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 254, 3, 9, 10, 92, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 257, 3, 10, 11, 95, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 260, 3, 11, 12, 98, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 263, 3, 12, 13, 101, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 266, 3, 13, 14, 104, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 269, 3, 14, 15, 107, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 272, 0, 89, 254, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 275, 0, 92, 257, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 278, 0, 95, 260, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 281, 0, 98, 263, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 284, 0, 101, 266, ncols, p);

            compute_prim_dd_electron_repulsion_2(buffer, 287, 3, 110, 44, 47, 137, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 296, 3, 113, 47, 50, 146, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 305, 3, 116, 50, 53, 155, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 314, 3, 119, 53, 56, 164, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 323, 3, 122, 56, 59, 173, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 332, 3, 125, 59, 62, 182, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 341, 3, 137, 68, 71, 209, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 350, 3, 146, 71, 74, 218, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 359, 3, 155, 74, 77, 227, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 368, 3, 164, 77, 80, 236, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 377, 3, 173, 80, 83, 245, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 386, 3, 89, 92, 257, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 392, 3, 92, 95, 260, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 400, 3, 95, 98, 263, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 406, 3, 98, 101, 266, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 412, 3, 101, 104, 269, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_8(buffer, 418, 0, 254, 386, ncols, p);

            compute_prim_pf_electron_repulsion_14(buffer, 421, 0, 257, 392, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 424, 0, 260, 400, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 427, 0, 263, 406, ncols, p);

            compute_prim_df_electron_repulsion_2(buffer, 430, 3, 272, 128, 137, 296, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 448, 3, 275, 137, 146, 305, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 466, 3, 278, 146, 155, 314, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 484, 3, 281, 155, 164, 323, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 502, 3, 284, 164, 173, 332, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 520, 3, 287, 191, 200, 341, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 538, 3, 296, 200, 209, 350, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 556, 3, 305, 209, 218, 359, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 574, 3, 314, 218, 227, 368, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 592, 3, 323, 227, 236, 377, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_6(buffer, 610, 3, 254, 257, 392, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_4(buffer, 622, 3, 257, 260, 400, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 634, 3, 260, 263, 406, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 643, 3, 263, 266, 412, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_9(buffer, 652, 0, 386, 610, ncols, p);

            compute_prim_pg_electron_repulsion_15(buffer, 655, 0, 3, 392, 622, 424, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 661, 0, 400, 634, ncols, p);

            compute_prim_dg_electron_repulsion_2(buffer, 664, 3, 418, 287, 296, 448, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 691, 3, 421, 296, 305, 466, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 718, 3, 424, 305, 314, 484, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 745, 3, 427, 314, 323, 502, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 772, 3, 448, 341, 350, 556, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 799, 3, 466, 350, 359, 574, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 826, 3, 484, 359, 368, 592, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_9(buffer, 853, 3, 386, 392, 622, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_10(buffer, 870, 3, 392, 400, 634, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 886, 3, 400, 406, 643, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_13(buffer, 899, 0, 3, 610, 853, 655, ncols, p);

            compute_prim_ph_electron_repulsion_14(buffer, 915, 0, 3, 622, 870, 661, ncols, p);

            compute_prim_ph_electron_repulsion_8(buffer, 930, 0, 634, 886, ncols, p);

            compute_prim_dh_electron_repulsion_2(buffer, 933, 3, 652, 430, 448, 691, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_7(buffer, 972, 0, 3, 655, 915, 448, 466, 718, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_8(buffer, 1022, 0, 3, 661, 930, 466, 484, 745, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_2(buffer, 1064, 3, 664, 520, 538, 772, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_2(buffer, 1103, 3, 691, 538, 556, 799, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_3(buffer, 1142, 0, 3, 718, 1022, 556, 574, 826, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_8(buffer, 1196, 3, 610, 622, 870, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_7(buffer, 1214, 3, 622, 634, 886, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_9(buffer, 1232, 0, 3, 853, 1196, 915, ncols, p);

            compute_prim_pi_electron_repulsion_10(buffer, 1268, 0, 870, 1214, ncols, p);

            compute_prim_di_electron_repulsion_5(buffer, 1280, 0, 3, 899, 1232, 664, 691, 972, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_6(buffer, 1388, 0, 3, 915, 1268, 691, 718, 1022, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_1(buffer, 1463, 0, 3, 972, 1388, 772, 799, 1142, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_7(buffer, 1625, 3, 853, 870, 1214, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_5(buffer, 1643, 0, 1196, 1625, ncols, p);

            compute_prim_dk_electron_repulsion_2(buffer, 1661, 0, 3, 1232, 1643, 933, 972, 1388, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 1745, 0, 3, 1280, 1661, 1064, 1103, 1463, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2105, 1745, 360, ncols);
        }
    }

    simdtrf::transform_fk(values, nvalues, buffer, 2105, nmax);
}

}  // namespace simdt2ceri
