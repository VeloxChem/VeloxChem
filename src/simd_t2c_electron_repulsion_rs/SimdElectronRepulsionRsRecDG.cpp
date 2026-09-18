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


#include "SimdElectronRepulsionRsRecDG.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_dg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_dg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1490, 1256, 180, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 6, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 14, 6, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 15, 16, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 16, 17, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 17, 18, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 18, 19, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 19, 20, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 112, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 115, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 118, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 121, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 124, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 127, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 130, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 133, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 136, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 145, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 154, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 163, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 172, 3, 17, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 181, 3, 18, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 190, 3, 19, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 199, 3, 20, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 208, 0, 3, 25, 145, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 226, 0, 3, 28, 154, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 244, 0, 3, 31, 163, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 262, 0, 3, 40, 181, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 280, 0, 3, 43, 190, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 298, 0, 3, 46, 199, 106, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 316, 3, 9, 10, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 322, 3, 10, 11, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 328, 3, 11, 12, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 334, 3, 17, 18, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 340, 3, 18, 19, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 346, 3, 19, 20, 133, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 352, 0, 3, 112, 316, 145, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 370, 0, 3, 115, 322, 154, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 388, 0, 3, 118, 328, 163, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 406, 0, 3, 124, 334, 181, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 424, 0, 3, 127, 340, 190, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 442, 0, 3, 130, 346, 199, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 460, 0, 3, 136, 352, 52, 58, 208, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 496, 0, 3, 145, 370, 58, 64, 226, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 532, 0, 3, 154, 388, 64, 70, 244, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 568, 0, 3, 172, 406, 82, 88, 262, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 604, 0, 3, 181, 424, 88, 94, 280, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 640, 0, 3, 190, 442, 94, 100, 298,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 676, 3, 112, 115, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 686, 3, 115, 118, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 696, 3, 124, 127, 340, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 706, 3, 127, 130, 346, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 716, 0, 3, 316, 676, 370, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 746, 0, 3, 322, 686, 388, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 776, 0, 3, 334, 696, 424, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 806, 0, 3, 340, 706, 442, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 836, 0, 3, 370, 746, 208, 226, 532,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 896, 0, 3, 424, 806, 262, 280, 640,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 956, 3, 316, 322, 686, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 971, 3, 334, 340, 706, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 986, 0, 3, 676, 956, 746, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1031, 0, 3, 696, 971, 806, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1076, 0, 3, 716, 986, 460, 496, 836,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1166, 0, 3, 776, 1031, 568, 604, 896,
                                                 ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1256, 1076, 180, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 1436, 1346, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 1436, 9, nmax);

    simdtrf::transform_g_inner(buffer, 1436, 1256, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 45 * nvalues, nvalues, buffer, 1436, 9, nmax);
}

}  // namespace simdt2ceri
