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


#include "SimdElectronRepulsionGeom10RsRecPH.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ph_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ph_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 3191, 2780, 378, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 7, 8, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 8, 9, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 9, 10, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 10, 11, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 11, 12, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 15, 16, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 16, 17, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 17, 18, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 18, 19, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 19, 20, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 118, 3, 7, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 121, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 124, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 127, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 130, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 133, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 136, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 139, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 142, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 145, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 148, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 151, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 154, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 157, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 160, 3, 8, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 169, 3, 9, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 178, 3, 10, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 187, 3, 11, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 196, 3, 12, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 205, 3, 16, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 214, 3, 17, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 223, 3, 18, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 232, 3, 19, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 241, 3, 20, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 250, 0, 3, 22, 160, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 268, 0, 3, 25, 169, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 286, 0, 3, 28, 178, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 304, 0, 3, 31, 187, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 322, 0, 3, 34, 196, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 340, 0, 3, 40, 205, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 358, 0, 3, 43, 214, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 376, 0, 3, 46, 223, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 394, 0, 3, 49, 232, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 412, 0, 3, 52, 241, 112, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 430, 3, 7, 8, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 436, 3, 8, 9, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 442, 3, 9, 10, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 448, 3, 10, 11, 133, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 454, 3, 11, 12, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 460, 3, 15, 16, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 466, 3, 16, 17, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 472, 3, 17, 18, 151, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 478, 3, 18, 19, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 484, 3, 19, 20, 157, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 490, 0, 3, 124, 436, 169, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 508, 0, 3, 127, 442, 178, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 526, 0, 3, 130, 448, 187, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 544, 0, 3, 133, 454, 196, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 562, 0, 3, 145, 466, 214, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 580, 0, 3, 148, 472, 223, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 598, 0, 3, 151, 478, 232, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 616, 0, 3, 154, 484, 241, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 634, 0, 3, 169, 508, 58, 64, 286, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 670, 0, 3, 178, 526, 64, 70, 304, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 706, 0, 3, 187, 544, 70, 76, 322, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 742, 0, 3, 214, 580, 88, 94, 376, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 778, 0, 3, 223, 598, 94, 100, 394,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 814, 0, 3, 232, 616, 100, 106, 412,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 850, 3, 118, 121, 430, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 860, 3, 121, 124, 436, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 870, 3, 124, 127, 442, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 880, 3, 127, 130, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 890, 3, 130, 133, 454, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 900, 3, 139, 142, 460, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 910, 3, 142, 145, 466, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 920, 3, 145, 148, 472, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 930, 3, 148, 151, 478, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 940, 3, 151, 154, 484, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 950, 0, 3, 436, 870, 508, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 980, 0, 3, 442, 880, 526, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1010, 0, 3, 448, 890, 544, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1040, 0, 3, 466, 920, 580, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1070, 0, 3, 472, 930, 598, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1100, 0, 3, 478, 940, 616, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1130, 0, 3, 490, 950, 250, 268, 634,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1190, 0, 3, 508, 980, 268, 286, 670,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1250, 0, 3, 526, 1010, 286, 304, 706,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1310, 0, 3, 562, 1040, 340, 358, 742,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1370, 0, 3, 580, 1070, 358, 376, 778,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1430, 0, 3, 598, 1100, 376, 394, 814,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1490, 3, 430, 436, 870, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1505, 3, 436, 442, 880, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1520, 3, 442, 448, 890, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1535, 3, 460, 466, 920, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1550, 3, 466, 472, 930, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1565, 3, 472, 478, 940, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1580, 0, 3, 870, 1505, 980, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1625, 0, 3, 880, 1520, 1010, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1670, 0, 3, 920, 1550, 1070, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1715, 0, 3, 930, 1565, 1100, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1760, 0, 3, 980, 1625, 634, 670, 1250,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1850, 0, 3, 1070, 1715, 742, 778, 1430,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1940, 3, 850, 860, 1490, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1961, 3, 870, 880, 1520, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1982, 3, 900, 910, 1535, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2003, 3, 920, 930, 1565, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 2024, 0, 3, 1505, 1961, 1625, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2087, 0, 3, 1550, 2003, 1715, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 2150, 0, 3, 1580, 2024, 1130, 1190,
                                                 1760, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2276, 0, 3, 1670, 2087, 1310, 1370,
                                                 1850, ncols, alpha, beta, p);

            simdgeo::geom_p_x(buffer, 2402, 1982, 2276, 1, 21, ncols, alpha);

            simdgeo::geom_p_y(buffer, 2465, 1982, 2276, 1, 21, ncols, alpha);

            simdgeo::geom_p_z(buffer, 2528, 1982, 2276, 1, 21, ncols, alpha);

            simdgeo::geom_p_x(buffer, 2591, 1940, 2150, 1, 21, ncols, alpha);

            simdgeo::geom_p_y(buffer, 2654, 1940, 2150, 1, 21, ncols, alpha);

            simdgeo::geom_p_z(buffer, 2717, 1940, 2150, 1, 21, ncols, alpha);

            simdfunc::contract_primitives(buffer, 2780, 2591, 189, ncols);

            simdfunc::contract_primitives(buffer, 2969, 2402, 189, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 3158, 2969, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 3158, 11, nmax);

    simdtrf::transform_h_inner(buffer, 3158, 3032, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 33 * nvalues, nvalues, buffer, 3158, 11, nmax);

    simdtrf::transform_h_inner(buffer, 3158, 3095, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 66 * nvalues, nvalues, buffer, 3158, 11, nmax);

    simdtrf::transform_h_inner(buffer, 3158, 2780, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 99 * nvalues, nvalues, buffer, 3158, 11, nmax);

    simdtrf::transform_h_inner(buffer, 3158, 2843, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 132 * nvalues, nvalues, buffer, 3158, 11, nmax);

    simdtrf::transform_h_inner(buffer, 3158, 2906, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 165 * nvalues, nvalues, buffer, 3158, 11, nmax);
}

}  // namespace simdt2ceri
