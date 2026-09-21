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


#include "SimdElectronRepulsionGeom10RsRecGD.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_gd_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_gd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 3473, 2858, 540, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 14, {1, 2, 3, 4, 5, 6, 7},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 7, 8, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 8, 9, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 9, 10, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 10, 11, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 11, 12, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 15, 16, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 16, 17, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 17, 18, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 18, 19, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 19, 20, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 124, 0, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 134, 0, 25, 28, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 144, 0, 28, 31, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 154, 0, 31, 34, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 164, 0, 34, 37, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 174, 0, 43, 46, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 184, 0, 46, 49, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 194, 0, 49, 52, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 52, 55, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 55, 58, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 224, 0, 64, 70, 144, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 239, 0, 70, 76, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 254, 0, 76, 82, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 269, 0, 94, 100, 194, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 284, 0, 100, 106, 204, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 299, 0, 106, 112, 214, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 314, 0, 124, 134, 224, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 335, 0, 134, 144, 239, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 356, 0, 144, 154, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 377, 0, 174, 184, 269, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 398, 0, 184, 194, 284, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 419, 0, 194, 204, 299, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 440, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 443, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 446, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 449, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 452, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 455, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 458, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 461, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 464, 3, 9, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 473, 3, 10, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 482, 3, 11, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 491, 3, 12, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 500, 3, 17, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 509, 3, 18, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 518, 3, 19, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 527, 3, 20, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 536, 0, 3, 28, 464, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 554, 0, 3, 31, 473, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 572, 0, 3, 34, 482, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 590, 0, 3, 37, 491, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 608, 0, 3, 49, 500, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 626, 0, 3, 52, 509, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 644, 0, 3, 55, 518, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 662, 0, 3, 58, 527, 118, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 680, 0, 3, 70, 554, 144, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 710, 0, 3, 76, 572, 154, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 740, 0, 3, 82, 590, 164, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 770, 0, 3, 100, 626, 194, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 800, 0, 3, 106, 644, 204, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 830, 0, 3, 112, 662, 214, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 860, 0, 3, 144, 710, 239, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 905, 0, 3, 154, 740, 254, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 950, 0, 3, 194, 800, 284, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 995, 0, 3, 204, 830, 299, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1040, 0, 3, 239, 905, 356, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1103, 0, 3, 284, 995, 419, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1166, 3, 9, 10, 443, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1172, 3, 10, 11, 446, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1178, 3, 11, 12, 449, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1184, 3, 17, 18, 455, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1190, 3, 18, 19, 458, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1196, 3, 19, 20, 461, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1202, 0, 3, 440, 1166, 473, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1220, 0, 3, 443, 1172, 482, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1238, 0, 3, 446, 1178, 491, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1256, 0, 3, 452, 1184, 509, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1274, 0, 3, 455, 1190, 518, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1292, 0, 3, 458, 1196, 527, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1310, 0, 3, 464, 1202, 64, 70, 554,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1346, 0, 3, 473, 1220, 70, 76, 572,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1382, 0, 3, 482, 1238, 76, 82, 590,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1418, 0, 3, 500, 1256, 94, 100, 626,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1454, 0, 3, 509, 1274, 100, 106, 644,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1490, 0, 3, 518, 1292, 106, 112, 662,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1526, 0, 3, 536, 1310, 124, 134, 680,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1586, 0, 3, 554, 1346, 134, 144, 710,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1646, 0, 3, 572, 1382, 144, 154, 740,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1706, 0, 3, 608, 1418, 174, 184, 770,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1766, 0, 3, 626, 1454, 184, 194, 800,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1826, 0, 3, 644, 1490, 194, 204, 830,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1886, 0, 3, 1310, 1346, 710, 1646, 224,
                                                 239, 905, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1976, 0, 3, 1418, 1454, 800, 1826, 269,
                                                 284, 995, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2066, 0, 3, 1526, 1586, 860, 1886, 314,
                                                 335, 1040, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2192, 0, 3, 1706, 1766, 950, 1976, 377,
                                                 398, 1103, ncols, alpha, beta, p);

            simdgeo::geom_g_x(buffer, 2318, 1706, 2192, 1, 6, ncols, alpha);

            simdgeo::geom_g_y(buffer, 2408, 1706, 2192, 1, 6, ncols, alpha);

            simdgeo::geom_g_z(buffer, 2498, 1706, 2192, 1, 6, ncols, alpha);

            simdgeo::geom_g_x(buffer, 2588, 1526, 2066, 1, 6, ncols, alpha);

            simdgeo::geom_g_y(buffer, 2678, 1526, 2066, 1, 6, ncols, alpha);

            simdgeo::geom_g_z(buffer, 2768, 1526, 2066, 1, 6, ncols, alpha);

            simdfunc::contract_primitives(buffer, 2858, 2588, 270, ncols);

            simdfunc::contract_primitives(buffer, 3128, 2318, 270, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 3398, 3128, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 3398, 5, nmax);

    simdtrf::transform_d_inner(buffer, 3398, 3218, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 45 * nvalues, nvalues, buffer, 3398, 5, nmax);

    simdtrf::transform_d_inner(buffer, 3398, 3308, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 90 * nvalues, nvalues, buffer, 3398, 5, nmax);

    simdtrf::transform_d_inner(buffer, 3398, 2858, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 135 * nvalues, nvalues, buffer, 3398, 5, nmax);

    simdtrf::transform_d_inner(buffer, 3398, 2948, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 180 * nvalues, nvalues, buffer, 3398, 5, nmax);

    simdtrf::transform_d_inner(buffer, 3398, 3038, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 225 * nvalues, nvalues, buffer, 3398, 5, nmax);
}

}  // namespace simdt2ceri
