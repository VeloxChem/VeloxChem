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


#include "SimdElectronRepulsionRecGK.hpp"

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
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gk_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 13673, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 48, 0, 7, 8, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 54, 0, 8, 9, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 9, 10, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 10, 11, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 11, 12, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 12, 13, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 13, 14, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 14, 15, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 15, 16, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 18, 21, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 21, 24, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 24, 27, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 27, 30, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 30, 33, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 33, 36, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 36, 39, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 39, 42, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 182, 0, 48, 54, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 197, 0, 54, 60, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 212, 0, 60, 66, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 227, 0, 66, 72, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 242, 0, 72, 78, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 257, 0, 78, 84, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 272, 0, 84, 90, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 287, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 290, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 293, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 296, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 299, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 302, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 305, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 308, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 311, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 314, 3, 8, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 323, 3, 9, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 332, 3, 10, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 341, 3, 11, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 350, 3, 12, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 359, 3, 13, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 368, 3, 14, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 377, 3, 15, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 386, 3, 16, 45, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 395, 0, 3, 18, 314, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 413, 0, 3, 21, 323, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 431, 0, 3, 24, 332, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 449, 0, 3, 27, 341, 66, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 467, 0, 3, 30, 350, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 485, 0, 3, 33, 359, 78, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 503, 0, 3, 36, 368, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 521, 0, 3, 39, 377, 90, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 539, 0, 3, 42, 386, 96, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 557, 0, 3, 54, 431, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 587, 0, 3, 60, 449, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 617, 0, 3, 66, 467, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 647, 0, 3, 72, 485, 142, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 677, 0, 3, 78, 503, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 707, 0, 3, 84, 521, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 737, 0, 3, 90, 539, 172, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 767, 0, 3, 102, 557, 182, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 812, 0, 3, 112, 587, 197, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 857, 0, 3, 122, 617, 212, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 902, 0, 3, 132, 647, 227, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 947, 0, 3, 142, 677, 242, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 992, 0, 3, 152, 707, 257, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1037, 0, 3, 162, 737, 272, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1082, 3, 8, 9, 290, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1088, 3, 9, 10, 293, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1094, 3, 10, 11, 296, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1100, 3, 11, 12, 299, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1106, 3, 12, 13, 302, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1112, 3, 13, 14, 305, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1118, 3, 14, 15, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1124, 3, 15, 16, 311, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1130, 0, 3, 287, 1082, 323, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1148, 0, 3, 290, 1088, 332, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1166, 0, 3, 293, 1094, 341, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1184, 0, 3, 296, 1100, 350, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1202, 0, 3, 299, 1106, 359, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1220, 0, 3, 302, 1112, 368, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1238, 0, 3, 305, 1118, 377, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1256, 0, 3, 308, 1124, 386, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1274, 0, 3, 323, 1148, 48, 54, 431,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1310, 0, 3, 332, 1166, 54, 60, 449,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1346, 0, 3, 341, 1184, 60, 66, 467,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1382, 0, 3, 350, 1202, 66, 72, 485,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1418, 0, 3, 359, 1220, 72, 78, 503,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1454, 0, 3, 368, 1238, 78, 84, 521,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1490, 0, 3, 377, 1256, 84, 90, 539,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1526, 0, 3, 431, 1310, 102, 112, 587,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1586, 0, 3, 449, 1346, 112, 122, 617,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1646, 0, 3, 467, 1382, 122, 132, 647,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1706, 0, 3, 485, 1418, 132, 142, 677,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1766, 0, 3, 503, 1454, 142, 152, 707,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1826, 0, 3, 521, 1490, 152, 162, 737,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1886, 0, 3, 1274, 1310, 587, 1586, 182,
                                                 197, 857, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1976, 0, 3, 1310, 1346, 617, 1646, 197,
                                                 212, 902, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2066, 0, 3, 1346, 1382, 647, 1706, 212,
                                                 227, 947, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2156, 0, 3, 1382, 1418, 677, 1766, 227,
                                                 242, 992, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2246, 0, 3, 1418, 1454, 707, 1826, 242,
                                                 257, 1037, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2336, 3, 287, 290, 1088, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2346, 3, 290, 293, 1094, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2356, 3, 293, 296, 1100, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2366, 3, 296, 299, 1106, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2376, 3, 299, 302, 1112, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2386, 3, 302, 305, 1118, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 2396, 3, 305, 308, 1124, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 2406, 0, 3, 1082, 2336, 1148, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2436, 0, 3, 1088, 2346, 1166, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2466, 0, 3, 1094, 2356, 1184, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2496, 0, 3, 1100, 2366, 1202, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2526, 0, 3, 1106, 2376, 1220, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2556, 0, 3, 1112, 2386, 1238, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2586, 0, 3, 1118, 2396, 1256, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2616, 0, 3, 1130, 2406, 395, 413, 1274,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2676, 0, 3, 1148, 2436, 413, 431, 1310,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2736, 0, 3, 1166, 2466, 431, 449, 1346,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2796, 0, 3, 1184, 2496, 449, 467, 1382,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2856, 0, 3, 1202, 2526, 467, 485, 1418,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2916, 0, 3, 1220, 2556, 485, 503, 1454,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2976, 0, 3, 1238, 2586, 503, 521, 1490,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3036, 0, 3, 1310, 2736, 557, 587, 1586,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3136, 0, 3, 1346, 2796, 587, 617, 1646,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3236, 0, 3, 1382, 2856, 617, 647, 1706,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3336, 0, 3, 1418, 2916, 647, 677, 1766,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 3436, 0, 3, 1454, 2976, 677, 707, 1826,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3536, 0, 3, 2616, 2676, 1526, 3036, 767,
                                                 812, 1886, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3686, 0, 3, 2676, 2736, 1586, 3136, 812,
                                                 857, 1976, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3836, 0, 3, 2736, 2796, 1646, 3236, 857,
                                                 902, 2066, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 3986, 0, 3, 2796, 2856, 1706, 3336, 902,
                                                 947, 2156, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4136, 0, 3, 2856, 2916, 1766, 3436, 947,
                                                 992, 2246, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4286, 3, 1082, 1088, 2346, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4301, 3, 1088, 1094, 2356, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4316, 3, 1094, 1100, 2366, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4331, 3, 1100, 1106, 2376, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4346, 3, 1106, 1112, 2386, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 4361, 3, 1112, 1118, 2396, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 4376, 0, 3, 2336, 4286, 2436, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4421, 0, 3, 2346, 4301, 2466, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4466, 0, 3, 2356, 4316, 2496, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4511, 0, 3, 2366, 4331, 2526, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4556, 0, 3, 2376, 4346, 2556, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 4601, 0, 3, 2386, 4361, 2586, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 4646, 0, 3, 2436, 4421, 1274, 1310,
                                                 2736, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4736, 0, 3, 2466, 4466, 1310, 1346,
                                                 2796, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4826, 0, 3, 2496, 4511, 1346, 1382,
                                                 2856, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4916, 0, 3, 2526, 4556, 1382, 1418,
                                                 2916, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 5006, 0, 3, 2556, 4601, 1418, 1454,
                                                 2976, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5096, 0, 3, 2736, 4736, 1526, 1586,
                                                 3136, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5246, 0, 3, 2796, 4826, 1586, 1646,
                                                 3236, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5396, 0, 3, 2856, 4916, 1646, 1706,
                                                 3336, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 5546, 0, 3, 2916, 5006, 1706, 1766,
                                                 3436, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 5696, 0, 3, 4646, 4736, 3136, 5246,
                                                 1886, 1976, 3836, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 5921, 0, 3, 4736, 4826, 3236, 5396,
                                                 1976, 2066, 3986, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 6146, 0, 3, 4826, 4916, 3336, 5546,
                                                 2066, 2156, 4136, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6371, 3, 2336, 2346, 4301, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6392, 3, 2346, 2356, 4316, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6413, 3, 2356, 2366, 4331, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6434, 3, 2366, 2376, 4346, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 6455, 3, 2376, 2386, 4361, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 6476, 0, 3, 4286, 6371, 4421, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6539, 0, 3, 4301, 6392, 4466, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6602, 0, 3, 4316, 6413, 4511, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6665, 0, 3, 4331, 6434, 4556, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 6728, 0, 3, 4346, 6455, 4601, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 6791, 0, 3, 4376, 6476, 2616, 2676,
                                                 4646, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 6917, 0, 3, 4421, 6539, 2676, 2736,
                                                 4736, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7043, 0, 3, 4466, 6602, 2736, 2796,
                                                 4826, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7169, 0, 3, 4511, 6665, 2796, 2856,
                                                 4916, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 7295, 0, 3, 4556, 6728, 2856, 2916,
                                                 5006, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 7421, 0, 3, 4736, 7043, 3036, 3136,
                                                 5246, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 7631, 0, 3, 4826, 7169, 3136, 3236,
                                                 5396, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 7841, 0, 3, 4916, 7295, 3236, 3336,
                                                 5546, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 8051, 0, 3, 6791, 6917, 5096, 7421,
                                                 3536, 3686, 5696, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 8366, 0, 3, 6917, 7043, 5246, 7631,
                                                 3686, 3836, 5921, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 8681, 0, 3, 7043, 7169, 5396, 7841,
                                                 3836, 3986, 6146, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 8996, 3, 4286, 4301, 6392, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9024, 3, 4301, 4316, 6413, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9052, 3, 4316, 4331, 6434, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 9080, 3, 4331, 4346, 6455, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 9108, 0, 3, 6371, 8996, 6539, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9192, 0, 3, 6392, 9024, 6602, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9276, 0, 3, 6413, 9052, 6665, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 9360, 0, 3, 6434, 9080, 6728, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 9444, 0, 3, 6539, 9192, 4646, 4736,
                                                 7043, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 9612, 0, 3, 6602, 9276, 4736, 4826,
                                                 7169, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 9780, 0, 3, 6665, 9360, 4826, 4916,
                                                 7295, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 9948, 0, 3, 7043, 9612, 5096, 5246,
                                                 7631, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 10228, 0, 3, 7169, 9780, 5246, 5396,
                                                 7841, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 10508, 0, 3, 9444, 9612, 7631, 10228,
                                                 5696, 5921, 8681, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 10928, 3, 6371, 6392, 9024, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 10964, 3, 6392, 6413, 9052, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 11000, 3, 6413, 6434, 9080, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 11036, 0, 3, 8996, 10928, 9192, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 11144, 0, 3, 9024, 10964, 9276, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 11252, 0, 3, 9052, 11000, 9360, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 11360, 0, 3, 9108, 11036, 6791, 6917,
                                                 9444, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 11576, 0, 3, 9192, 11144, 6917, 7043,
                                                 9612, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 11792, 0, 3, 9276, 11252, 7043, 7169,
                                                 9780, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 12008, 0, 3, 9612, 11792, 7421, 7631,
                                                 10228, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 12368, 0, 3, 11360, 11576, 9948, 12008,
                                                 8051, 8366, 10508, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 12908, 12368, 540, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 13448, 12908, 15, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 13448, 15, nmax);
}

}  // namespace simdt2ceri
