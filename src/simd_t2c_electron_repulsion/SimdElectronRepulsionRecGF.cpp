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


#include "SimdElectronRepulsionRecGF.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gf_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gf_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1741, nvalues);

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

            compute_prim_fs_electron_repulsion_0(buffer, 62, 0, 14, 17, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 72, 0, 17, 20, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 82, 0, 20, 23, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 92, 0, 23, 26, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 102, 0, 32, 38, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 117, 0, 38, 44, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 132, 0, 44, 50, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 147, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 150, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 153, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 156, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 159, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 162, 3, 8, 17, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 171, 3, 9, 20, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 180, 3, 10, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 189, 3, 11, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 198, 3, 12, 29, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 207, 0, 3, 14, 162, 32, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 225, 0, 3, 17, 171, 38, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 243, 0, 3, 20, 180, 44, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 261, 0, 3, 23, 189, 50, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 279, 0, 3, 26, 198, 56, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 297, 0, 3, 38, 243, 72, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 327, 0, 3, 44, 261, 82, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 357, 0, 3, 50, 279, 92, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 387, 0, 3, 62, 297, 102, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 432, 0, 3, 72, 327, 117, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 477, 0, 3, 82, 357, 132, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 522, 3, 8, 9, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 528, 3, 9, 10, 153, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 534, 3, 10, 11, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 540, 3, 11, 12, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 546, 0, 3, 147, 522, 171, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 564, 0, 3, 150, 528, 180, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 582, 0, 3, 153, 534, 189, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 600, 0, 3, 156, 540, 198, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 618, 0, 3, 171, 564, 32, 38, 243, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 654, 0, 3, 180, 582, 38, 44, 261, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 690, 0, 3, 189, 600, 44, 50, 279, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 726, 0, 3, 243, 654, 62, 72, 327, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 786, 0, 3, 261, 690, 72, 82, 357, ncols,
                                                 alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 846, 0, 3, 618, 654, 327, 786, 102, 117,
                                                 477, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 936, 3, 147, 150, 528, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 946, 3, 150, 153, 534, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 956, 3, 153, 156, 540, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 966, 0, 3, 522, 936, 564, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 996, 0, 3, 528, 946, 582, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1026, 0, 3, 534, 956, 600, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1056, 0, 3, 546, 966, 207, 225, 618,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1116, 0, 3, 564, 996, 225, 243, 654,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1176, 0, 3, 582, 1026, 243, 261, 690,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1236, 0, 3, 654, 1176, 297, 327, 786,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 1336, 0, 3, 1056, 1116, 726, 1236, 387,
                                                 432, 846, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1486, 1336, 150, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 1636, 1486, 15, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 1636, 7, nmax);
}

}  // namespace simdt2ceri
