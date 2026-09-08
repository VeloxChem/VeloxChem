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
#include "SimdTransformFG.hpp"

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

    auto buffer = CSimdMatrix(758, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 35, 0, 7, 8, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 38, 0, 8, 9, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 41, 0, 9, 10, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 44, 0, 10, 11, 29, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 47, 0, 11, 12, 32, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 50, 0, 14, 17, 35, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 53, 0, 17, 20, 38, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 56, 0, 20, 23, 41, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 59, 0, 23, 26, 44, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 62, 0, 26, 29, 47, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 65, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 68, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 71, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 74, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 77, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 80, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 83, 3, 11, 29, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 86, 3, 20, 38, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 95, 3, 23, 41, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 104, 3, 26, 44, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 113, 3, 29, 47, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 122, 3, 38, 56, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 131, 3, 41, 59, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 140, 3, 44, 62, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 149, 3, 9, 10, 68, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 152, 3, 10, 11, 71, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 155, 3, 11, 12, 74, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_3(buffer, 158, 0, 65, 149, ncols, p);

            compute_prim_pd_electron_repulsion_3(buffer, 164, 0, 68, 152, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 170, 0, 71, 155, ncols, p);

            compute_prim_dd_electron_repulsion_2(buffer, 173, 3, 77, 35, 38, 95, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_3(buffer, 182, 0, 3, 80, 164, 38, 41, 104, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_7(buffer, 197, 0, 3, 83, 170, 41, 44, 113, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 209, 3, 86, 50, 53, 122, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 218, 3, 95, 53, 56, 131, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_3(buffer, 227, 0, 3, 104, 197, 56, 59, 140, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_4(buffer, 251, 3, 65, 68, 152, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 257, 3, 68, 71, 155, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_11(buffer, 263, 0, 3, 149, 251, 164, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 278, 0, 152, 257, ncols, p);

            compute_prim_df_electron_repulsion_5(buffer, 281, 0, 3, 158, 263, 86, 95, 182, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_6(buffer, 317, 0, 3, 164, 278, 95, 104, 197, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_1(buffer, 347, 0, 3, 182, 317, 122, 131, 227, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_2(buffer, 410, 3, 149, 152, 257, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_8(buffer, 416, 0, 251, 410, ncols, p);

            compute_prim_dg_electron_repulsion_4(buffer, 422, 0, 3, 263, 416, 173, 182, 317, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 458, 0, 3, 281, 422, 209, 218, 347, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 608, 458, 150, ncols);
        }
    }

    simdtrf::transform_fg(values, nvalues, buffer, 608, nmax);
}

}  // namespace simdt2ceri
