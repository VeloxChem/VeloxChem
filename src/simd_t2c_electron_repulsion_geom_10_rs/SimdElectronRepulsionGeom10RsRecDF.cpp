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


#include "SimdElectronRepulsionGeom10RsRecDF.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_df_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_df_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 2192, 1790, 360, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 13, {1, 2, 3, 4, 5, 6}, ncols,
                                            fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 35, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 19, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 50, 0, 7, 8, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 56, 0, 8, 9, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 62, 0, 9, 10, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 68, 0, 10, 11, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 74, 0, 14, 15, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 80, 0, 15, 16, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 86, 0, 16, 17, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 92, 0, 17, 18, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 98, 0, 20, 23, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 108, 0, 23, 26, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 118, 0, 26, 29, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 128, 0, 35, 38, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 138, 0, 38, 41, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 148, 0, 41, 44, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 170, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 173, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 176, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 179, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 182, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 185, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 188, 3, 9, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 197, 3, 10, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 206, 3, 11, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 215, 3, 16, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 224, 3, 17, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 233, 3, 18, 47, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 242, 0, 3, 23, 188, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 260, 0, 3, 26, 197, 62, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 278, 0, 3, 29, 206, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 296, 0, 3, 38, 215, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 314, 0, 3, 41, 224, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 332, 0, 3, 44, 233, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 350, 0, 3, 50, 242, 98, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 380, 0, 3, 56, 260, 108, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 410, 0, 3, 62, 278, 118, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 440, 0, 3, 74, 296, 128, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 470, 0, 3, 80, 314, 138, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 500, 0, 3, 86, 332, 148, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 530, 3, 7, 8, 161, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 536, 3, 8, 9, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 542, 3, 9, 10, 167, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 548, 3, 10, 11, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 554, 3, 14, 15, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 560, 3, 15, 16, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 566, 3, 16, 17, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 572, 3, 17, 18, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 578, 0, 3, 161, 536, 188, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 596, 0, 3, 164, 542, 197, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 614, 0, 3, 167, 548, 206, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 632, 0, 3, 176, 560, 215, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 650, 0, 3, 179, 566, 224, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 668, 0, 3, 182, 572, 233, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 686, 0, 3, 188, 596, 50, 56, 260, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 722, 0, 3, 197, 614, 56, 62, 278, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 758, 0, 3, 215, 650, 74, 80, 314, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 794, 0, 3, 224, 668, 80, 86, 332, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 830, 0, 3, 260, 722, 98, 108, 410,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 890, 0, 3, 314, 794, 128, 138, 500,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 950, 3, 158, 161, 536, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 960, 3, 164, 167, 548, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 970, 3, 173, 176, 560, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 980, 3, 179, 182, 572, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 990, 0, 3, 530, 950, 578, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1020, 0, 3, 542, 960, 614, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1050, 0, 3, 554, 970, 632, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1080, 0, 3, 566, 980, 668, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1110, 0, 3, 596, 1020, 242, 260, 722,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1170, 0, 3, 650, 1080, 296, 314, 794,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1230, 0, 3, 686, 1110, 350, 380, 830,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1330, 0, 3, 758, 1170, 440, 470, 890,
                                                 ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 1430, 1050, 1330, 1, 10, ncols, alpha);

            simdgeo::geom_d_y(buffer, 1490, 1050, 1330, 1, 10, ncols, alpha);

            simdgeo::geom_d_z(buffer, 1550, 1050, 1330, 1, 10, ncols, alpha);

            simdgeo::geom_d_x(buffer, 1610, 990, 1230, 1, 10, ncols, alpha);

            simdgeo::geom_d_y(buffer, 1670, 990, 1230, 1, 10, ncols, alpha);

            simdgeo::geom_d_z(buffer, 1730, 990, 1230, 1, 10, ncols, alpha);

            simdfunc::contract_primitives(buffer, 1790, 1610, 180, ncols);

            simdfunc::contract_primitives(buffer, 1970, 1430, 180, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 2150, 1970, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 2150, 7, nmax);

    simdtrf::transform_f_inner(buffer, 2150, 2030, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 35 * nvalues, nvalues, buffer, 2150, 7, nmax);

    simdtrf::transform_f_inner(buffer, 2150, 2090, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 70 * nvalues, nvalues, buffer, 2150, 7, nmax);

    simdtrf::transform_f_inner(buffer, 2150, 1790, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 105 * nvalues, nvalues, buffer, 2150, 7, nmax);

    simdtrf::transform_f_inner(buffer, 2150, 1850, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 140 * nvalues, nvalues, buffer, 2150, 7, nmax);

    simdtrf::transform_f_inner(buffer, 2150, 1910, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 175 * nvalues, nvalues, buffer, 2150, 7, nmax);
}

}  // namespace simdt2ceri
