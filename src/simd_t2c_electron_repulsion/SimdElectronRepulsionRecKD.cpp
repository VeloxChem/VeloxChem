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


#include "SimdElectronRepulsionRecKD.hpp"

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
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdTransformKD.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kd_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1297, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 43, 0, 7, 8, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 46, 0, 8, 9, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 49, 0, 9, 10, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 10, 11, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 11, 12, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 58, 0, 12, 13, 37, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 61, 0, 13, 14, 40, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 64, 0, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 70, 0, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 76, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 82, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 88, 0, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 94, 0, 31, 34, 58, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 100, 0, 34, 37, 61, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 106, 0, 43, 46, 76, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 115, 0, 46, 49, 82, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_5(buffer, 124, 0, 49, 52, 88, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 136, 0, 52, 55, 94, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 145, 0, 55, 58, 100, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 154, 0, 64, 70, 106, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 166, 0, 70, 76, 115, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_6(buffer, 178, 0, 76, 82, 124, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_7(buffer, 196, 0, 82, 88, 136, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 212, 0, 88, 94, 145, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_9(buffer, 225, 0, 106, 115, 178, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_10(buffer, 240, 0, 115, 124, 196, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_11(buffer, 263, 0, 124, 136, 212, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_1(buffer, 284, 0, 154, 166, 225, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_2(buffer, 299, 0, 166, 178, 240, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_3(buffer, 314, 0, 178, 196, 263, ncols, alpha, beta, p);

            compute_prim_pp_electron_repulsion_2(buffer, 344, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 347, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 350, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 353, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 356, 3, 13, 37, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 359, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 362, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 365, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 368, 3, 31, 55, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 371, 3, 34, 58, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 374, 3, 37, 61, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 377, 3, 46, 76, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 380, 3, 49, 82, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 383, 3, 52, 88, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 386, 3, 55, 94, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 389, 3, 58, 100, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 392, 3, 76, 115, ncols, p);

            compute_prim_gp_electron_repulsion_12(buffer, 395, 3, 82, 124, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 398, 3, 88, 136, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 401, 3, 94, 145, ncols, p);

            compute_prim_hp_electron_repulsion_12(buffer, 404, 3, 115, 178, ncols, p);

            compute_prim_hp_electron_repulsion_13(buffer, 407, 0, 3, 124, 398, 196, ncols, p);

            compute_prim_hp_electron_repulsion_14(buffer, 431, 3, 136, 212, ncols, p);

            compute_prim_ip_electron_repulsion_6(buffer, 440, 0, 3, 178, 407, 240, ncols, p);

            compute_prim_ip_electron_repulsion_7(buffer, 492, 0, 3, 196, 431, 263, ncols, p);

            compute_prim_kp_electron_repulsion_1(buffer, 529, 0, 3, 240, 492, 314, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 610, 3, 344, 43, 46, 362, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 613, 3, 347, 46, 49, 365, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 616, 3, 350, 49, 52, 368, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 619, 3, 353, 52, 55, 371, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 622, 3, 356, 55, 58, 374, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 625, 0, 3, 359, 610, 64, 70, 377, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 634, 0, 3, 362, 613, 70, 76, 380, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 643, 0, 3, 365, 616, 76, 82, 383, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 652, 0, 3, 368, 619, 82, 88, 386, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 661, 0, 3, 371, 622, 88, 94, 389, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 670, 0, 3, 610, 613, 380, 643, 106, 115, 395, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_31(buffer, 685, 0, 3, 613, 616, 383, 652, 115, 124, 398, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_32(buffer, 700, 0, 3, 616, 619, 386, 661, 124, 136, 401, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_22(buffer, 715, 0, 3, 625, 634, 392, 670, 154, 166, 404, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_23(buffer, 739, 0, 3, 634, 643, 395, 685, 166, 178, 407, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_24(buffer, 763, 0, 3, 643, 652, 398, 700, 178, 196, 431, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_9(buffer, 793, 0, 3, 670, 685, 407, 763, 225, 240, 492, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 865, 0, 3, 715, 739, 440, 793, 284, 299, 529, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1081, 865, 216, ncols);
        }
    }

    simdtrf::transform_kd(values, nvalues, buffer, 1081, nmax);
}

}  // namespace simdt2ceri
