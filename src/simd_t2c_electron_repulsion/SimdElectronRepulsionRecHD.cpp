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


#include "SimdElectronRepulsionRecHD.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

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
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hd_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1393, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7}, ncols,
                                            fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 35, 0, 7, 8, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 41, 0, 8, 9, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 47, 0, 9, 10, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 53, 0, 10, 11, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 59, 0, 11, 12, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 65, 0, 14, 17, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 75, 0, 17, 20, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 85, 0, 20, 23, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 95, 0, 23, 26, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 26, 29, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 115, 0, 35, 41, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 130, 0, 41, 47, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 145, 0, 47, 53, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 160, 0, 65, 75, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 181, 0, 75, 85, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 202, 0, 85, 95, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 223, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 226, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 229, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 232, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 235, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 244, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 253, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 262, 3, 12, 32, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 271, 0, 3, 20, 235, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 289, 0, 3, 23, 244, 47, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 307, 0, 3, 26, 253, 53, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 325, 0, 3, 29, 262, 59, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 343, 0, 3, 41, 289, 85, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 373, 0, 3, 47, 307, 95, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 403, 0, 3, 53, 325, 105, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 433, 0, 3, 85, 373, 130, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 478, 0, 3, 95, 403, 145, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 523, 0, 3, 130, 478, 202, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 586, 3, 9, 10, 226, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 592, 3, 10, 11, 229, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 598, 3, 11, 12, 232, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 604, 0, 3, 223, 586, 244, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 622, 0, 3, 226, 592, 253, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 640, 0, 3, 229, 598, 262, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 658, 0, 3, 235, 604, 35, 41, 289, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 694, 0, 3, 244, 622, 41, 47, 307, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 730, 0, 3, 253, 640, 47, 53, 325, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 766, 0, 3, 271, 658, 65, 75, 343, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 826, 0, 3, 289, 694, 75, 85, 373, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 886, 0, 3, 307, 730, 85, 95, 403, ncols,
                                                 alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 946, 0, 3, 658, 694, 373, 886, 115, 130,
                                                 478, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 1036, 0, 3, 766, 826, 433, 946, 160,
                                                 181, 523, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1162, 1036, 126, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 1288, 1162, 21, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 1288, 5, nmax);
}

}  // namespace simdt2ceri
