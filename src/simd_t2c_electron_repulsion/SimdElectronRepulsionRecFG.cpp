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


#include "SimdElectronRepulsionRecFG.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fg_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1525, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 115, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 118, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 121, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 124, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 127, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 136, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 145, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 154, 3, 12, 32, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 163, 0, 3, 20, 127, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 181, 0, 3, 23, 136, 47, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 199, 0, 3, 26, 145, 53, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 217, 0, 3, 29, 154, 59, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 235, 0, 3, 41, 181, 85, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 265, 0, 3, 47, 199, 95, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 295, 0, 3, 53, 217, 105, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 325, 3, 9, 10, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 331, 3, 10, 11, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 337, 3, 11, 12, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 343, 0, 3, 115, 325, 136, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 361, 0, 3, 118, 331, 145, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 379, 0, 3, 121, 337, 154, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 397, 0, 3, 127, 343, 35, 41, 181, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 433, 0, 3, 136, 361, 41, 47, 199, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 469, 0, 3, 145, 379, 47, 53, 217, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 505, 0, 3, 163, 397, 65, 75, 235, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 565, 0, 3, 181, 433, 75, 85, 265, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 625, 0, 3, 199, 469, 85, 95, 295, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 685, 3, 115, 118, 331, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 695, 3, 118, 121, 337, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 705, 0, 3, 325, 685, 361, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 735, 0, 3, 331, 695, 379, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 765, 0, 3, 343, 705, 163, 181, 433,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 825, 0, 3, 361, 735, 181, 199, 469,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 885, 0, 3, 433, 825, 235, 265, 625,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 985, 3, 325, 331, 695, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1000, 0, 3, 685, 985, 735, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1045, 0, 3, 705, 1000, 397, 433, 825,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 1135, 0, 3, 765, 1045, 505, 565, 885,
                                                 ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1285, 1135, 150, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 1435, 1285, 10, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 1435, 9, nmax);
}

}  // namespace simdt2ceri
