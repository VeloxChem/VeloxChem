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


#include "SimdElectronRepulsionRecGH.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gh_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(5556, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 40, 0, 7, 8, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 46, 0, 8, 9, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 9, 10, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 10, 11, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 11, 12, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 12, 13, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 13, 14, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 82, 0, 16, 19, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 92, 0, 19, 22, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 22, 25, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 25, 28, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 28, 31, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 31, 34, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 142, 0, 40, 46, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 157, 0, 46, 52, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 172, 0, 52, 58, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 187, 0, 58, 64, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 202, 0, 64, 70, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 217, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 220, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 223, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 226, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 229, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 232, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 235, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 238, 3, 8, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 247, 3, 9, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 256, 3, 10, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 265, 3, 11, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 274, 3, 12, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 283, 3, 13, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 292, 3, 14, 37, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 301, 0, 3, 16, 238, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 319, 0, 3, 19, 247, 46, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 337, 0, 3, 22, 256, 52, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 355, 0, 3, 25, 265, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 373, 0, 3, 28, 274, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 391, 0, 3, 31, 283, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 409, 0, 3, 34, 292, 76, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 427, 0, 3, 46, 337, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 457, 0, 3, 52, 355, 102, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 487, 0, 3, 58, 373, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 517, 0, 3, 64, 391, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 547, 0, 3, 70, 409, 132, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 577, 0, 3, 82, 427, 142, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 622, 0, 3, 92, 457, 157, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 667, 0, 3, 102, 487, 172, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 712, 0, 3, 112, 517, 187, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 757, 0, 3, 122, 547, 202, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 802, 3, 8, 9, 220, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 808, 3, 9, 10, 223, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 814, 3, 10, 11, 226, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 820, 3, 11, 12, 229, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 826, 3, 12, 13, 232, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 832, 3, 13, 14, 235, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 838, 0, 3, 217, 802, 247, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 856, 0, 3, 220, 808, 256, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 874, 0, 3, 223, 814, 265, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 892, 0, 3, 226, 820, 274, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 910, 0, 3, 229, 826, 283, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 928, 0, 3, 232, 832, 292, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 946, 0, 3, 247, 856, 40, 46, 337, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 982, 0, 3, 256, 874, 46, 52, 355, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1018, 0, 3, 265, 892, 52, 58, 373,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1054, 0, 3, 274, 910, 58, 64, 391,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1090, 0, 3, 283, 928, 64, 70, 409,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1126, 0, 3, 337, 982, 82, 92, 457,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1186, 0, 3, 355, 1018, 92, 102, 487,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1246, 0, 3, 373, 1054, 102, 112, 517,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1306, 0, 3, 391, 1090, 112, 122, 547,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1366, 0, 3, 946, 982, 457, 1186, 142,
                                                 157, 667, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1456, 0, 3, 982, 1018, 487, 1246, 157,
                                                 172, 712, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1546, 0, 3, 1018, 1054, 517, 1306, 172,
                                                 187, 757, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1636, 3, 217, 220, 808, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1646, 3, 220, 223, 814, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1656, 3, 223, 226, 820, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1666, 3, 226, 229, 826, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1676, 3, 229, 232, 832, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1686, 0, 3, 802, 1636, 856, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1716, 0, 3, 808, 1646, 874, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1746, 0, 3, 814, 1656, 892, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1776, 0, 3, 820, 1666, 910, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1806, 0, 3, 826, 1676, 928, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1836, 0, 3, 838, 1686, 301, 319, 946,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1896, 0, 3, 856, 1716, 319, 337, 982,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1956, 0, 3, 874, 1746, 337, 355, 1018,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2016, 0, 3, 892, 1776, 355, 373, 1054,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2076, 0, 3, 910, 1806, 373, 391, 1090,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2136, 0, 3, 982, 1956, 427, 457, 1186,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2236, 0, 3, 1018, 2016, 457, 487, 1246,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2336, 0, 3, 1054, 2076, 487, 517, 1306,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 2436, 0, 3, 1836, 1896, 1126, 2136, 577,
                                                 622, 1366, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 2586, 0, 3, 1896, 1956, 1186, 2236, 622,
                                                 667, 1456, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 2736, 0, 3, 1956, 2016, 1246, 2336, 667,
                                                 712, 1546, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2886, 3, 802, 808, 1646, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2901, 3, 808, 814, 1656, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2916, 3, 814, 820, 1666, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2931, 3, 820, 826, 1676, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2946, 0, 3, 1636, 2886, 1716, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2991, 0, 3, 1646, 2901, 1746, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3036, 0, 3, 1656, 2916, 1776, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3081, 0, 3, 1666, 2931, 1806, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 3126, 0, 3, 1716, 2991, 946, 982, 1956,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3216, 0, 3, 1746, 3036, 982, 1018, 2016,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3306, 0, 3, 1776, 3081, 1018, 1054,
                                                 2076, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3396, 0, 3, 1956, 3216, 1126, 1186,
                                                 2236, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3546, 0, 3, 2016, 3306, 1186, 1246,
                                                 2336, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 3696, 0, 3, 3126, 3216, 2236, 3546,
                                                 1366, 1456, 2736, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3921, 3, 1636, 1646, 2901, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3942, 3, 1646, 1656, 2916, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 3963, 3, 1656, 1666, 2931, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 3984, 0, 3, 2886, 3921, 2991, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4047, 0, 3, 2901, 3942, 3036, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4110, 0, 3, 2916, 3963, 3081, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 4173, 0, 3, 2946, 3984, 1836, 1896,
                                                 3126, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4299, 0, 3, 2991, 4047, 1896, 1956,
                                                 3216, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4425, 0, 3, 3036, 4110, 1956, 2016,
                                                 3306, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4551, 0, 3, 3216, 4425, 2136, 2236,
                                                 3546, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 4761, 0, 3, 4173, 4299, 3396, 4551,
                                                 2436, 2586, 3696, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 5076, 4761, 315, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 5391, 5076, 15, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 5391, 11, nmax);
}

}  // namespace simdt2ceri
