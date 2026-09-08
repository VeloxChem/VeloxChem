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


#include "SimdElectronRepulsionRecGD.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gd_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 836, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 6, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 29, 0, 7, 8, 14, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 35, 0, 8, 9, 17, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 41, 0, 9, 10, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 47, 0, 10, 11, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 53, 0, 11, 12, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 59, 0, 14, 17, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 69, 0, 17, 20, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 79, 0, 20, 23, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 89, 0, 29, 35, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 104, 0, 35, 41, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 119, 0, 41, 47, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 134, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 137, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 140, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 143, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 146, 3, 9, 17, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 155, 3, 10, 20, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 164, 3, 11, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 173, 3, 12, 26, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 182, 0, 3, 17, 155, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 200, 0, 3, 20, 164, 47, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 218, 0, 3, 23, 173, 53, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 236, 0, 3, 41, 200, 69, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 266, 0, 3, 47, 218, 79, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 296, 0, 3, 69, 266, 119, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 341, 3, 9, 10, 137, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 347, 3, 10, 11, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 353, 3, 11, 12, 143, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 359, 0, 3, 134, 341, 155, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 377, 0, 3, 137, 347, 164, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 395, 0, 3, 140, 353, 173, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 413, 0, 3, 146, 359, 29, 35, 182, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 449, 0, 3, 155, 377, 35, 41, 200, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 485, 0, 3, 164, 395, 41, 47, 218, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 521, 0, 3, 200, 485, 59, 69, 266, ncols,
                                                 alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 581, 0, 3, 413, 449, 236, 521, 89, 104,
                                                 296, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 671, 581, 90, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 761, 671, 15, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 761, 5, nmax);
}

}  // namespace simdt2ceri
