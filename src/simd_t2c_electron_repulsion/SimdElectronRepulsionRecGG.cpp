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


#include "SimdElectronRepulsionRecGG.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gg_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 3246, 2886, 225, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 8, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 37, 0, 7, 8, 16, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 43, 0, 8, 9, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 49, 0, 9, 10, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 55, 0, 10, 11, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 61, 0, 11, 12, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 12, 13, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 13, 14, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 79, 0, 16, 19, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 89, 0, 19, 22, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 99, 0, 22, 25, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 109, 0, 25, 28, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 119, 0, 28, 31, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 129, 0, 37, 43, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 144, 0, 43, 49, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 159, 0, 49, 55, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 174, 0, 55, 61, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 189, 0, 61, 67, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 204, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 207, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 210, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 213, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 216, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 219, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 222, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 231, 3, 10, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 240, 3, 11, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 249, 3, 12, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 258, 3, 13, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 267, 3, 14, 34, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 276, 0, 3, 19, 231, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 294, 0, 3, 22, 240, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 312, 0, 3, 25, 249, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 330, 0, 3, 28, 258, 67, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 348, 0, 3, 31, 267, 73, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 366, 0, 3, 49, 294, 89, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 396, 0, 3, 55, 312, 99, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 426, 0, 3, 61, 330, 109, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 456, 0, 3, 67, 348, 119, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 486, 0, 3, 89, 396, 159, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 531, 0, 3, 99, 426, 174, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 576, 0, 3, 109, 456, 189, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 621, 3, 9, 10, 207, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 627, 3, 10, 11, 210, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 633, 3, 11, 12, 213, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 639, 3, 12, 13, 216, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 645, 3, 13, 14, 219, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 651, 0, 3, 204, 621, 231, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 669, 0, 3, 207, 627, 240, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 687, 0, 3, 210, 633, 249, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 705, 0, 3, 213, 639, 258, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 723, 0, 3, 216, 645, 267, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 741, 0, 3, 222, 651, 37, 43, 276, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 777, 0, 3, 231, 669, 43, 49, 294, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 813, 0, 3, 240, 687, 49, 55, 312, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 849, 0, 3, 249, 705, 55, 61, 330, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 885, 0, 3, 258, 723, 61, 67, 348, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 921, 0, 3, 294, 813, 79, 89, 396, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 981, 0, 3, 312, 849, 89, 99, 426, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1041, 0, 3, 330, 885, 99, 109, 456,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1101, 0, 3, 741, 777, 366, 921, 129,
                                                 144, 486, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1191, 0, 3, 777, 813, 396, 981, 144,
                                                 159, 531, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1281, 0, 3, 813, 849, 426, 1041, 159,
                                                 174, 576, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1371, 3, 204, 207, 627, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1381, 3, 207, 210, 633, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1391, 3, 210, 213, 639, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1401, 3, 213, 216, 645, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1411, 0, 3, 621, 1371, 669, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1441, 0, 3, 627, 1381, 687, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1471, 0, 3, 633, 1391, 705, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1501, 0, 3, 639, 1401, 723, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1531, 0, 3, 669, 1441, 276, 294, 813,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1591, 0, 3, 687, 1471, 294, 312, 849,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1651, 0, 3, 705, 1501, 312, 330, 885,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1711, 0, 3, 813, 1591, 366, 396, 981,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1811, 0, 3, 849, 1651, 396, 426, 1041,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 1911, 0, 3, 1531, 1591, 981, 1811, 486,
                                                 531, 1281, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2061, 3, 621, 627, 1381, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2076, 3, 627, 633, 1391, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2091, 3, 633, 639, 1401, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2106, 0, 3, 1371, 2061, 1441, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2151, 0, 3, 1381, 2076, 1471, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2196, 0, 3, 1391, 2091, 1501, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2241, 0, 3, 1411, 2106, 741, 777, 1531,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2331, 0, 3, 1441, 2151, 777, 813, 1591,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2421, 0, 3, 1471, 2196, 813, 849, 1651,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2511, 0, 3, 1591, 2421, 921, 981, 1811,
                                                 ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 2661, 0, 3, 2241, 2331, 1711, 2511,
                                                 1101, 1191, 1911, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2886, 2661, 225, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 3111, 2886, 15, nmax);

    simdtrf::transform_g_outer_tri(values, nvalues, buffer, 3111, nmax);
}

}  // namespace simdt2ceri
