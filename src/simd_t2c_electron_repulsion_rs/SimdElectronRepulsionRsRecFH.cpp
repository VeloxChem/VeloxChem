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


#include "SimdElectronRepulsionRsRecFH.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_fh_electron_repulsion(double               *values,
                                 const size_t          nvalues,
                                 const CBasisFunction &bra,
                                 const CBasisFunction &ket,
                                 const CSimdMatrix    &coordinates,
                                 CSimdMatrix          &buffer,
                                 const double          omega) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_fh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 5148, 4618, 420, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 15, {1, 2, 3, 4, 5, 6, 7, 8},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 23, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 7, 8, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 8, 9, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 9, 10, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 10, 11, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 11, 12, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 12, 13, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 16, 17, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 17, 18, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 18, 19, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 19, 20, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 20, 21, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 21, 22, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 138, 0, 24, 27, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 148, 0, 27, 30, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 30, 33, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 33, 36, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 36, 39, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 45, 48, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 48, 51, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 51, 54, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 54, 57, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 57, 60, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 238, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 241, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 244, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 247, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 250, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 253, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 256, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 259, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 262, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 265, 3, 23, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 268, 3, 9, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 277, 3, 10, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 286, 3, 11, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 295, 3, 12, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 304, 3, 13, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 313, 3, 18, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 322, 3, 19, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 331, 3, 20, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 340, 3, 21, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 349, 3, 22, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 358, 0, 3, 27, 268, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 376, 0, 3, 30, 277, 78, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 394, 0, 3, 33, 286, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 412, 0, 3, 36, 295, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 430, 0, 3, 39, 304, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 448, 0, 3, 48, 313, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 466, 0, 3, 51, 322, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 484, 0, 3, 54, 331, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 502, 0, 3, 57, 340, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 520, 0, 3, 60, 349, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 538, 0, 3, 66, 358, 138, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 568, 0, 3, 72, 376, 148, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 598, 0, 3, 78, 394, 158, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 628, 0, 3, 84, 412, 168, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 658, 0, 3, 90, 430, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 688, 0, 3, 102, 448, 188, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 718, 0, 3, 108, 466, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 748, 0, 3, 114, 484, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 778, 0, 3, 120, 502, 218, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 808, 0, 3, 126, 520, 228, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 838, 3, 9, 10, 241, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 844, 3, 10, 11, 244, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 850, 3, 11, 12, 247, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 856, 3, 12, 13, 250, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 862, 3, 18, 19, 256, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 868, 3, 19, 20, 259, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 874, 3, 20, 21, 262, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 880, 3, 21, 22, 265, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 886, 0, 3, 238, 838, 277, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 904, 0, 3, 241, 844, 286, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 922, 0, 3, 244, 850, 295, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 940, 0, 3, 247, 856, 304, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 958, 0, 3, 253, 862, 322, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 976, 0, 3, 256, 868, 331, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 994, 0, 3, 259, 874, 340, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1012, 0, 3, 262, 880, 349, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1030, 0, 3, 268, 886, 66, 72, 376,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1066, 0, 3, 277, 904, 72, 78, 394,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1102, 0, 3, 286, 922, 78, 84, 412,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1138, 0, 3, 295, 940, 84, 90, 430,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1174, 0, 3, 313, 958, 102, 108, 466,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1210, 0, 3, 322, 976, 108, 114, 484,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1246, 0, 3, 331, 994, 114, 120, 502,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1282, 0, 3, 340, 1012, 120, 126, 520,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1318, 0, 3, 376, 1066, 138, 148, 598,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1378, 0, 3, 394, 1102, 148, 158, 628,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1438, 0, 3, 412, 1138, 158, 168, 658,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1498, 0, 3, 466, 1210, 188, 198, 748,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1558, 0, 3, 484, 1246, 198, 208, 778,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1618, 0, 3, 502, 1282, 208, 218, 808,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1678, 3, 238, 241, 844, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1688, 3, 241, 244, 850, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1698, 3, 244, 247, 856, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1708, 3, 253, 256, 868, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1718, 3, 256, 259, 874, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1728, 3, 259, 262, 880, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1738, 0, 3, 838, 1678, 904, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1768, 0, 3, 844, 1688, 922, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1798, 0, 3, 850, 1698, 940, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1828, 0, 3, 862, 1708, 976, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1858, 0, 3, 868, 1718, 994, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1888, 0, 3, 874, 1728, 1012, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1918, 0, 3, 886, 1738, 358, 376, 1066,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1978, 0, 3, 904, 1768, 376, 394, 1102,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2038, 0, 3, 922, 1798, 394, 412, 1138,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2098, 0, 3, 958, 1828, 448, 466, 1210,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2158, 0, 3, 976, 1858, 466, 484, 1246,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2218, 0, 3, 994, 1888, 484, 502, 1282,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2278, 0, 3, 1030, 1918, 538, 568, 1318,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2378, 0, 3, 1066, 1978, 568, 598, 1378,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2478, 0, 3, 1102, 2038, 598, 628, 1438,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2578, 0, 3, 1174, 2098, 688, 718, 1498,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2678, 0, 3, 1210, 2158, 718, 748, 1558,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2778, 0, 3, 1246, 2218, 748, 778, 1618,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2878, 3, 838, 844, 1688, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2893, 3, 844, 850, 1698, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2908, 3, 862, 868, 1718, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2923, 3, 868, 874, 1728, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2938, 0, 3, 1678, 2878, 1768, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2983, 0, 3, 1688, 2893, 1798, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3028, 0, 3, 1708, 2908, 1858, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3073, 0, 3, 1718, 2923, 1888, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 3118, 0, 3, 1738, 2938, 1030, 1066,
                                                 1978, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3208, 0, 3, 1768, 2983, 1066, 1102,
                                                 2038, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3298, 0, 3, 1828, 3028, 1174, 1210,
                                                 2158, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3388, 0, 3, 1858, 3073, 1210, 1246,
                                                 2218, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3478, 0, 3, 1978, 3208, 1318, 1378,
                                                 2478, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3628, 0, 3, 2158, 3388, 1498, 1558,
                                                 2778, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3778, 3, 1678, 1688, 2893, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3799, 3, 1708, 1718, 2923, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 3820, 0, 3, 2878, 3778, 2983, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 3883, 0, 3, 2908, 3799, 3073, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 3946, 0, 3, 2938, 3820, 1918, 1978,
                                                 3208, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4072, 0, 3, 3028, 3883, 2098, 2158,
                                                 3388, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4198, 0, 3, 3118, 3946, 2278, 2378,
                                                 3478, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4408, 0, 3, 3298, 4072, 2578, 2678,
                                                 3628, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4618, 4198, 420, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 5038, 4828, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 5038, 11, nmax);

    simdtrf::transform_h_inner(buffer, 5038, 4618, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 77 * nvalues, nvalues, buffer, 5038, 11, nmax);
}

}  // namespace simdt2ceri
