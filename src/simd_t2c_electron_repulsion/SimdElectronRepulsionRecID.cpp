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


#include "SimdElectronRepulsionRecID.hpp"

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
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdTransformID.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_id_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_id_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(922, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 8, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 37, 0, 7, 8, 16, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 40, 0, 8, 9, 19, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 43, 0, 9, 10, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 46, 0, 10, 11, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 49, 0, 11, 12, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 12, 13, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 13, 14, 34, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 58, 0, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 64, 0, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 70, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 79, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 85, 0, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 91, 0, 37, 40, 58, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 100, 0, 40, 43, 64, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 109, 0, 43, 46, 70, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_5(buffer, 121, 0, 46, 49, 79, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 133, 0, 49, 52, 85, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_8(buffer, 142, 0, 58, 64, 109, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_9(buffer, 154, 0, 64, 70, 121, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_10(buffer, 171, 0, 70, 79, 133, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_1(buffer, 187, 0, 91, 100, 142, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_2(buffer, 199, 0, 100, 109, 154, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_3(buffer, 211, 0, 109, 121, 171, ncols, alpha, beta, p);

            compute_prim_pp_electron_repulsion_2(buffer, 234, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 237, 3, 10, 22, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 240, 3, 11, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 243, 3, 12, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 246, 3, 13, 31, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 249, 3, 19, 43, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 252, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 255, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 258, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 261, 3, 31, 55, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 264, 3, 43, 64, ncols, p);

            compute_prim_fp_electron_repulsion_11(buffer, 267, 3, 46, 70, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 270, 3, 49, 79, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 273, 3, 52, 85, ncols, p);

            compute_prim_gp_electron_repulsion_12(buffer, 276, 3, 64, 109, ncols, p);

            compute_prim_gp_electron_repulsion_13(buffer, 279, 0, 3, 70, 270, 121, ncols, p);

            compute_prim_gp_electron_repulsion_14(buffer, 299, 3, 79, 133, ncols, p);

            compute_prim_hp_electron_repulsion_6(buffer, 308, 0, 3, 109, 279, 154, ncols, p);

            compute_prim_hp_electron_repulsion_7(buffer, 346, 0, 3, 121, 299, 171, ncols, p);

            compute_prim_ip_electron_repulsion_1(buffer, 376, 0, 3, 154, 346, 211, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 439, 3, 234, 37, 40, 249, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 442, 3, 237, 40, 43, 252, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 445, 3, 240, 43, 46, 255, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 448, 3, 243, 46, 49, 258, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 451, 3, 246, 49, 52, 261, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 454, 0, 3, 252, 445, 58, 64, 267, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_20(buffer, 463, 0, 3, 255, 448, 64, 70, 270, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_21(buffer, 472, 0, 3, 258, 451, 70, 79, 273, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 481, 0, 3, 439, 442, 264, 454, 91, 100, 276, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_21(buffer, 496, 0, 3, 442, 445, 267, 463, 100, 109, 279, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_22(buffer, 511, 0, 3, 445, 448, 270, 472, 109, 121, 299, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_9(buffer, 532, 0, 3, 454, 463, 279, 511, 142, 154, 346, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 586, 0, 3, 481, 496, 308, 532, 187, 199, 376, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 754, 586, 168, ncols);
        }
    }

    simdtrf::transform_id(values, nvalues, buffer, 754, nmax);
}

}  // namespace simdt2ceri
