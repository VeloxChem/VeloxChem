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


#include "SimdElectronRepulsionRecDH.hpp"

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
compute_dh_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_dh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1328, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 32, 0, 7, 8, 17, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 38, 0, 8, 9, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 44, 0, 9, 10, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 50, 0, 10, 11, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 56, 0, 11, 12, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 62, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 65, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 68, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 71, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 74, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 77, 3, 8, 17, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 86, 3, 9, 20, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 95, 3, 10, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 104, 3, 11, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 113, 3, 12, 29, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 122, 0, 3, 14, 77, 32, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 140, 0, 3, 17, 86, 38, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 158, 0, 3, 20, 95, 44, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 176, 0, 3, 23, 104, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 194, 0, 3, 26, 113, 56, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 212, 3, 8, 9, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 218, 3, 9, 10, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 224, 3, 10, 11, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 230, 3, 11, 12, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 236, 0, 3, 62, 212, 86, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 254, 0, 3, 65, 218, 95, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 272, 0, 3, 68, 224, 104, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 290, 0, 3, 71, 230, 113, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 308, 0, 3, 86, 254, 32, 38, 158, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 344, 0, 3, 95, 272, 38, 44, 176, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 380, 0, 3, 104, 290, 44, 50, 194, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 416, 3, 62, 65, 218, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 426, 3, 65, 68, 224, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 436, 3, 68, 71, 230, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 446, 0, 3, 212, 416, 254, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 476, 0, 3, 218, 426, 272, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 506, 0, 3, 224, 436, 290, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 536, 0, 3, 236, 446, 122, 140, 308,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 596, 0, 3, 254, 476, 140, 158, 344,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 656, 0, 3, 272, 506, 158, 176, 380,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 716, 3, 212, 218, 426, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 731, 3, 218, 224, 436, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 746, 0, 3, 416, 716, 476, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 791, 0, 3, 426, 731, 506, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 836, 0, 3, 476, 791, 308, 344, 656,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 926, 3, 416, 426, 731, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 947, 0, 3, 716, 926, 791, ncols, p);

            compute_prim_dh_electron_repulsion_0(buffer, 1010, 0, 3, 746, 947, 536, 596, 836,
                                                 ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1136, 1010, 126, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 1262, 1136, 6, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 1262, 11, nmax);
}

}  // namespace simdt2ceri
