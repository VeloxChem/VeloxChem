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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdTransformHD.hpp"

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

    auto buffer = CSimdMatrix(632, nvalues);

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

            compute_prim_fs_electron_repulsion_4(buffer, 50, 0, 14, 17, 35, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 56, 0, 17, 20, 38, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 62, 0, 20, 23, 41, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 71, 0, 23, 26, 44, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 80, 0, 26, 29, 47, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_6(buffer, 86, 0, 35, 38, 62, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 95, 0, 38, 41, 71, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_5(buffer, 107, 0, 41, 44, 80, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_1(buffer, 119, 0, 50, 56, 86, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_2(buffer, 128, 0, 56, 62, 95, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_3(buffer, 137, 0, 62, 71, 107, ncols, alpha, beta, p);

            compute_prim_pp_electron_repulsion_2(buffer, 154, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 157, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 160, 3, 11, 29, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 163, 3, 20, 38, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 166, 3, 23, 41, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 169, 3, 26, 44, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 172, 3, 29, 47, ncols, p);

            compute_prim_fp_electron_repulsion_11(buffer, 175, 3, 38, 62, ncols, p);

            compute_prim_fp_electron_repulsion_12(buffer, 178, 0, 3, 41, 169, 71, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 194, 3, 44, 80, ncols, p);

            compute_prim_gp_electron_repulsion_6(buffer, 203, 0, 3, 62, 178, 95, ncols, p);

            compute_prim_gp_electron_repulsion_7(buffer, 229, 0, 3, 71, 194, 107, ncols, p);

            compute_prim_hp_electron_repulsion_1(buffer, 252, 0, 3, 95, 229, 137, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 299, 3, 154, 35, 38, 166, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 302, 3, 157, 38, 41, 169, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 305, 3, 160, 41, 44, 172, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 308, 0, 3, 163, 299, 50, 56, 175, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_14(buffer, 317, 0, 3, 166, 302, 56, 62, 178, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_15(buffer, 326, 0, 3, 169, 305, 62, 71, 194, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_8(buffer, 341, 0, 3, 299, 302, 178, 326, 86, 95, 229, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 380, 0, 3, 308, 317, 203, 341, 119, 128, 252, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 506, 380, 126, ncols);
        }
    }

    simdtrf::transform_hd(values, nvalues, buffer, 506, nmax);
}

}  // namespace simdt2ceri
