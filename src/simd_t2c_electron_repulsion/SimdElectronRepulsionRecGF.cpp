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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformGF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gf_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
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

    auto buffer = CSimdMatrix(742, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 32, 0, 7, 8, 17, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 35, 0, 8, 9, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 38, 0, 9, 10, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 41, 0, 10, 11, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 44, 0, 11, 12, 29, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 47, 0, 14, 17, 35, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 50, 0, 17, 20, 38, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 56, 0, 20, 23, 41, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_3(buffer, 65, 0, 23, 26, 44, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 73, 0, 32, 35, 50, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_2(buffer, 79, 0, 35, 38, 56, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_3(buffer, 85, 0, 38, 41, 65, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 94, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 97, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 100, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 103, 3, 12, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 106, 3, 9, 20, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 109, 3, 10, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 112, 3, 11, 26, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 115, 3, 14, 32, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 118, 3, 17, 35, ncols, p);

            compute_prim_dp_electron_repulsion_8(buffer, 121, 0, 3, 20, 109, 38, ncols, p);

            compute_prim_dp_electron_repulsion_8(buffer, 131, 0, 3, 23, 112, 41, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 141, 3, 26, 44, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 150, 3, 35, 50, ncols, p);

            compute_prim_fp_electron_repulsion_8(buffer, 159, 0, 3, 38, 131, 56, ncols, p);

            compute_prim_fp_electron_repulsion_9(buffer, 174, 0, 3, 41, 141, 65, ncols, p);

            compute_prim_gp_electron_repulsion_2(buffer, 186, 3, 47, 73, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 198, 3, 50, 79, ncols, p);

            compute_prim_gp_electron_repulsion_4(buffer, 210, 0, 3, 56, 174, 85, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 232, 3, 8, 9, 97, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 235, 3, 9, 10, 100, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 238, 3, 10, 11, 103, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 241, 0, 94, 232, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 244, 0, 97, 235, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 247, 0, 100, 238, ncols, p);

            compute_prim_dd_electron_repulsion_10(buffer, 250, 3, 106, 32, 35, 121, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_6(buffer, 253, 0, 3, 109, 247, 35, 38, 131, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 269, 3, 112, 38, 41, 141, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_5(buffer, 278, 0, 3, 121, 253, 47, 50, 159, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_6(buffer, 308, 0, 3, 131, 269, 50, 56, 174, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_1(buffer, 331, 0, 3, 250, 253, 159, 308, 73, 79, 210, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_9(buffer, 388, 3, 241, 115, 118, 250, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_10(buffer, 391, 3, 244, 118, 121, 253, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_11(buffer, 394, 3, 247, 121, 131, 269, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_4(buffer, 403, 0, 3, 253, 394, 150, 159, 308, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 442, 0, 3, 388, 391, 278, 403, 186, 198, 331, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 592, 442, 150, ncols);
        }
    }

    simdtrf::transform_gf(values, nvalues, buffer, 592, nmax);
}

}  // namespace simdt2ceri
