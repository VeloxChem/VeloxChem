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


#include "SimdElectronRepulsionRsRecDH.hpp"

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
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_dh_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_dh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 2584, 2266, 252, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 118, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 121, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 124, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 127, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 130, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 133, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 136, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 139, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 142, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 145, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 148, 3, 8, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 157, 3, 9, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 166, 3, 10, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 175, 3, 11, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 184, 3, 12, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 193, 3, 16, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 202, 3, 17, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 211, 3, 18, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 220, 3, 19, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 229, 3, 20, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 238, 0, 3, 22, 148, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 256, 0, 3, 25, 157, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 274, 0, 3, 28, 166, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 292, 0, 3, 31, 175, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 310, 0, 3, 34, 184, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 328, 0, 3, 40, 193, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 346, 0, 3, 43, 202, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 364, 0, 3, 46, 211, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 382, 0, 3, 49, 220, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 400, 0, 3, 52, 229, 112, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 418, 3, 8, 9, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 424, 3, 9, 10, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 430, 3, 10, 11, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 436, 3, 11, 12, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 442, 3, 16, 17, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 448, 3, 17, 18, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 454, 3, 18, 19, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 460, 3, 19, 20, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 466, 0, 3, 118, 418, 157, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 484, 0, 3, 121, 424, 166, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 502, 0, 3, 124, 430, 175, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 520, 0, 3, 127, 436, 184, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 538, 0, 3, 133, 442, 202, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 556, 0, 3, 136, 448, 211, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 574, 0, 3, 139, 454, 220, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 592, 0, 3, 142, 460, 229, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 610, 0, 3, 157, 484, 58, 64, 274, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 646, 0, 3, 166, 502, 64, 70, 292, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 682, 0, 3, 175, 520, 70, 76, 310, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 718, 0, 3, 202, 556, 88, 94, 364, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 754, 0, 3, 211, 574, 94, 100, 382,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 790, 0, 3, 220, 592, 100, 106, 400,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 826, 3, 118, 121, 424, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 836, 3, 121, 124, 430, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 846, 3, 124, 127, 436, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 856, 3, 133, 136, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 866, 3, 136, 139, 454, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 876, 3, 139, 142, 460, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 886, 0, 3, 418, 826, 484, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 916, 0, 3, 424, 836, 502, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 946, 0, 3, 430, 846, 520, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 976, 0, 3, 442, 856, 556, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1006, 0, 3, 448, 866, 574, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1036, 0, 3, 454, 876, 592, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1066, 0, 3, 466, 886, 238, 256, 610,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1126, 0, 3, 484, 916, 256, 274, 646,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1186, 0, 3, 502, 946, 274, 292, 682,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1246, 0, 3, 538, 976, 328, 346, 718,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1306, 0, 3, 556, 1006, 346, 364, 754,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1366, 0, 3, 574, 1036, 364, 382, 790,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1426, 3, 418, 424, 836, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1441, 3, 424, 430, 846, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1456, 3, 442, 448, 866, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1471, 3, 448, 454, 876, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1486, 0, 3, 826, 1426, 916, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1531, 0, 3, 836, 1441, 946, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1576, 0, 3, 856, 1456, 1006, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1621, 0, 3, 866, 1471, 1036, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1666, 0, 3, 916, 1531, 610, 646, 1186,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1756, 0, 3, 1006, 1621, 718, 754, 1366,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1846, 3, 826, 836, 1441, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1867, 3, 856, 866, 1471, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 1888, 0, 3, 1426, 1846, 1531, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 1951, 0, 3, 1456, 1867, 1621, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 2014, 0, 3, 1486, 1888, 1066, 1126,
                                                 1666, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2140, 0, 3, 1576, 1951, 1246, 1306,
                                                 1756, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2266, 2014, 252, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 2518, 2392, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 2518, 11, nmax);

    simdtrf::transform_h_inner(buffer, 2518, 2266, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 55 * nvalues, nvalues, buffer, 2518, 11, nmax);
}

}  // namespace simdt2ceri
