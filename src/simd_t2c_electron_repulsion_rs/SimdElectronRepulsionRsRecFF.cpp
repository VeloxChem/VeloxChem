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


#include "SimdElectronRepulsionRsRecFF.hpp"

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
#include "SimdTransformF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_ff_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ff_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1548, 1278, 200, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 170, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 173, 3, 19, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 176, 3, 9, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 185, 3, 10, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 194, 3, 11, 32, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 203, 3, 16, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 212, 3, 17, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 221, 3, 18, 47, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 230, 0, 3, 23, 176, 56, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 248, 0, 3, 26, 185, 62, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 266, 0, 3, 29, 194, 68, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 284, 0, 3, 38, 203, 80, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 302, 0, 3, 41, 212, 86, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 320, 0, 3, 44, 221, 92, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 338, 0, 3, 50, 230, 98, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 368, 0, 3, 56, 248, 108, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 398, 0, 3, 62, 266, 118, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 428, 0, 3, 74, 284, 128, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 458, 0, 3, 80, 302, 138, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 488, 0, 3, 86, 320, 148, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 518, 3, 9, 10, 161, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 524, 3, 10, 11, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 530, 3, 16, 17, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 536, 3, 17, 18, 173, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 542, 0, 3, 158, 518, 185, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 560, 0, 3, 161, 524, 194, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 578, 0, 3, 167, 530, 212, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 596, 0, 3, 170, 536, 221, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 614, 0, 3, 176, 542, 50, 56, 248, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 650, 0, 3, 185, 560, 56, 62, 266, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 686, 0, 3, 203, 578, 74, 80, 302, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 722, 0, 3, 212, 596, 80, 86, 320, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 758, 0, 3, 248, 650, 98, 108, 398,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 818, 0, 3, 302, 722, 128, 138, 488,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 878, 3, 158, 161, 524, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 888, 3, 167, 170, 536, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 898, 0, 3, 518, 878, 560, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 928, 0, 3, 530, 888, 596, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 958, 0, 3, 542, 898, 230, 248, 650,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1018, 0, 3, 578, 928, 284, 302, 722,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1078, 0, 3, 614, 958, 338, 368, 758,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1178, 0, 3, 686, 1018, 428, 458, 818,
                                                 ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1278, 1078, 200, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 1478, 1378, 10, 1, nmax);

    simdtrf::transform_f_outer_tri(values, nvalues, buffer, 1478, nmax);

    simdtrf::transform_f_inner(buffer, 1478, 1278, 10, 1, nmax);

    simdtrf::transform_f_outer_tri(values + 49 * nvalues, nvalues, buffer, 1478, nmax);
}

}  // namespace simdt2ceri
