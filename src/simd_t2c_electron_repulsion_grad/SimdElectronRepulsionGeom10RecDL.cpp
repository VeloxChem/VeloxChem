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


#include "SimdElectronRepulsionGeom10RecDL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionGeom10VrrRecDL.hpp"
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
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_dl_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_dl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 11463, 10551, 810, nvalues);

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
                                            10, 11}, ncols, fj, mu);

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

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 7, 8, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 8, 9, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 9, 10, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 10, 11, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 11, 12, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 12, 13, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 13, 14, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 14, 15, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 99, 0, 15, 16, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 18, 21, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 115, 0, 21, 24, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 24, 27, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 27, 30, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 30, 33, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 33, 36, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 36, 39, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 39, 42, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 42, 45, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 195, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 198, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 201, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 204, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 207, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 210, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 213, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 216, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 219, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 222, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 225, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 234, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 243, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 252, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 261, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 270, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 279, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 288, 3, 16, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 297, 0, 3, 24, 225, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 315, 0, 3, 27, 234, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 333, 0, 3, 30, 243, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 351, 0, 3, 33, 252, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 369, 0, 3, 36, 261, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 387, 0, 3, 39, 270, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 405, 0, 3, 42, 279, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 423, 0, 3, 45, 288, 99, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 441, 0, 3, 57, 315, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 471, 0, 3, 63, 333, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 501, 0, 3, 69, 351, 145, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 531, 0, 3, 75, 369, 155, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 561, 0, 3, 81, 387, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 591, 0, 3, 87, 405, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 621, 0, 3, 93, 423, 185, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 651, 3, 7, 8, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 657, 3, 8, 9, 201, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 663, 3, 9, 10, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 669, 3, 10, 11, 207, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 675, 3, 11, 12, 210, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 681, 3, 12, 13, 213, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 687, 3, 13, 14, 216, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 693, 3, 14, 15, 219, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 699, 3, 15, 16, 222, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 705, 0, 3, 201, 663, 234, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 723, 0, 3, 204, 669, 243, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 741, 0, 3, 207, 675, 252, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 759, 0, 3, 210, 681, 261, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 777, 0, 3, 213, 687, 270, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 795, 0, 3, 216, 693, 279, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 813, 0, 3, 219, 699, 288, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 831, 0, 3, 225, 705, 51, 57, 315, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 867, 0, 3, 234, 723, 57, 63, 333, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 903, 0, 3, 243, 741, 63, 69, 351, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 939, 0, 3, 252, 759, 69, 75, 369, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 975, 0, 3, 261, 777, 75, 81, 387, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1011, 0, 3, 270, 795, 81, 87, 405,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1047, 0, 3, 279, 813, 87, 93, 423,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1083, 0, 3, 297, 831, 105, 115, 441,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1143, 0, 3, 315, 867, 115, 125, 471,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1203, 0, 3, 333, 903, 125, 135, 501,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1263, 0, 3, 351, 939, 135, 145, 531,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1323, 0, 3, 369, 975, 145, 155, 561,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1383, 0, 3, 387, 1011, 155, 165, 591,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1443, 0, 3, 405, 1047, 165, 175, 621,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1503, 3, 195, 198, 657, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1513, 3, 198, 201, 663, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1523, 3, 201, 204, 669, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1533, 3, 204, 207, 675, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1543, 3, 207, 210, 681, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1553, 3, 210, 213, 687, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1563, 3, 213, 216, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1573, 3, 216, 219, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1583, 0, 3, 663, 1523, 723, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1613, 0, 3, 669, 1533, 741, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1643, 0, 3, 675, 1543, 759, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1673, 0, 3, 681, 1553, 777, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1703, 0, 3, 687, 1563, 795, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1733, 0, 3, 693, 1573, 813, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1763, 0, 3, 705, 1583, 297, 315, 867,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1823, 0, 3, 723, 1613, 315, 333, 903,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1883, 0, 3, 741, 1643, 333, 351, 939,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1943, 0, 3, 759, 1673, 351, 369, 975,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2003, 0, 3, 777, 1703, 369, 387, 1011,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2063, 0, 3, 795, 1733, 387, 405, 1047,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2123, 0, 3, 867, 1823, 441, 471, 1203,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2223, 0, 3, 903, 1883, 471, 501, 1263,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2323, 0, 3, 939, 1943, 501, 531, 1323,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2423, 0, 3, 975, 2003, 531, 561, 1383,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2523, 0, 3, 1011, 2063, 561, 591, 1443,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2623, 3, 651, 657, 1513, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2638, 3, 657, 663, 1523, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2653, 3, 663, 669, 1533, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2668, 3, 669, 675, 1543, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2683, 3, 675, 681, 1553, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2698, 3, 681, 687, 1563, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2713, 3, 687, 693, 1573, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2728, 0, 3, 1523, 2653, 1613, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2773, 0, 3, 1533, 2668, 1643, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2818, 0, 3, 1543, 2683, 1673, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2863, 0, 3, 1553, 2698, 1703, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2908, 0, 3, 1563, 2713, 1733, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2953, 0, 3, 1583, 2728, 831, 867, 1823,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3043, 0, 3, 1613, 2773, 867, 903, 1883,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3133, 0, 3, 1643, 2818, 903, 939, 1943,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3223, 0, 3, 1673, 2863, 939, 975, 2003,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3313, 0, 3, 1703, 2908, 975, 1011, 2063,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3403, 0, 3, 1763, 2953, 1083, 1143,
                                                 2123, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3553, 0, 3, 1823, 3043, 1143, 1203,
                                                 2223, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3703, 0, 3, 1883, 3133, 1203, 1263,
                                                 2323, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3853, 0, 3, 1943, 3223, 1263, 1323,
                                                 2423, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 4003, 0, 3, 2003, 3313, 1323, 1383,
                                                 2523, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4153, 3, 1503, 1513, 2638, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4174, 3, 1513, 1523, 2653, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4195, 3, 1523, 1533, 2668, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4216, 3, 1533, 1543, 2683, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4237, 3, 1543, 1553, 2698, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4258, 3, 1553, 1563, 2713, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 4279, 0, 3, 2653, 4195, 2773, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4342, 0, 3, 2668, 4216, 2818, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4405, 0, 3, 2683, 4237, 2863, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4468, 0, 3, 2698, 4258, 2908, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 4531, 0, 3, 2728, 4279, 1763, 1823,
                                                 3043, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4657, 0, 3, 2773, 4342, 1823, 1883,
                                                 3133, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4783, 0, 3, 2818, 4405, 1883, 1943,
                                                 3223, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4909, 0, 3, 2863, 4468, 1943, 2003,
                                                 3313, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 5035, 0, 3, 3043, 4657, 2123, 2223,
                                                 3703, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 5245, 0, 3, 3133, 4783, 2223, 2323,
                                                 3853, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 5455, 0, 3, 3223, 4909, 2323, 2423,
                                                 4003, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5665, 3, 2623, 2638, 4174, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5693, 3, 2638, 2653, 4195, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5721, 3, 2653, 2668, 4216, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5749, 3, 2668, 2683, 4237, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5777, 3, 2683, 2698, 4258, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 5805, 0, 3, 4195, 5721, 4342, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5889, 0, 3, 4216, 5749, 4405, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5973, 0, 3, 4237, 5777, 4468, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 6057, 0, 3, 4279, 5805, 2953, 3043,
                                                 4657, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6225, 0, 3, 4342, 5889, 3043, 3133,
                                                 4783, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6393, 0, 3, 4405, 5973, 3133, 3223,
                                                 4909, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 6561, 0, 3, 4531, 6057, 3403, 3553,
                                                 5035, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 6841, 0, 3, 4657, 6225, 3553, 3703,
                                                 5245, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 7121, 0, 3, 4783, 6393, 3703, 3853,
                                                 5455, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7401, 3, 4153, 4174, 5693, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7437, 3, 4174, 4195, 5721, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7473, 3, 4195, 4216, 5749, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7509, 3, 4216, 4237, 5777, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 7545, 0, 3, 5693, 7437, 5805, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7653, 0, 3, 5721, 7473, 5889, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7761, 0, 3, 5749, 7509, 5973, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 7869, 0, 3, 5805, 7653, 4531, 4657,
                                                 6225, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 8085, 0, 3, 5889, 7761, 4657, 4783,
                                                 6393, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 8301, 0, 3, 6225, 8085, 5035, 5245,
                                                 7121, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8661, 3, 5665, 5693, 7437, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8706, 3, 5721, 5749, 7509, ncols, alpha,
                                                 beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 8751, 0, 3, 7401, 8661, 7545, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 8886, 0, 3, 7473, 8706, 7761, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 9021, 0, 3, 7653, 8886, 6057, 6225,
                                                 8085, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 9291, 0, 3, 7869, 9021, 6561, 6841,
                                                 8301, ncols, alpha, beta, p);

            compute_prim_geom_10_dl_electron_repulsion_0(buffer, 9741, 8751, 9291, ncols,
                                                         alpha);

            compute_prim_geom_10_dl_electron_repulsion_1(buffer, 10011, 8751, 9291, ncols,
                                                         alpha);

            compute_prim_geom_10_dl_electron_repulsion_2(buffer, 10281, 8751, 9291, ncols,
                                                         alpha);

            simdfunc::contract_primitives(buffer, 10551, 9741, 810, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 11361, 10551, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 11361, 17, nmax);

    simdtrf::transform_l_inner(buffer, 11361, 10821, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 85 * nvalues, nvalues, buffer, 11361, 17, nmax);

    simdtrf::transform_l_inner(buffer, 11361, 11091, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 170 * nvalues, nvalues, buffer, 11361, 17, nmax);
}

}  // namespace simdt2ceri
