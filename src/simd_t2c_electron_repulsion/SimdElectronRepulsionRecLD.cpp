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


#include "SimdElectronRepulsionRecLD.hpp"

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
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdTransformLD.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ld_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ld_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1787, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 45, 0, 7, 8, 18, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 48, 0, 8, 9, 21, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 9, 10, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 54, 0, 10, 11, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 57, 0, 11, 12, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 60, 0, 12, 13, 33, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 63, 0, 13, 14, 36, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 66, 0, 14, 15, 39, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 69, 0, 15, 16, 42, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 72, 0, 18, 21, 51, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 78, 0, 21, 24, 54, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 84, 0, 24, 27, 57, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 90, 0, 27, 30, 60, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 96, 0, 30, 33, 63, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 102, 0, 33, 36, 66, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 108, 0, 36, 39, 69, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 114, 0, 45, 48, 72, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 123, 0, 48, 51, 78, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 132, 0, 51, 54, 84, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 141, 0, 54, 57, 90, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 150, 0, 57, 60, 96, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 159, 0, 60, 63, 102, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 168, 0, 63, 66, 108, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_14(buffer, 177, 0, 72, 78, 132, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 189, 0, 78, 84, 141, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_7(buffer, 202, 0, 84, 90, 150, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 218, 0, 90, 96, 159, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 231, 0, 96, 102, 168, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_15(buffer, 244, 0, 114, 123, 177, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_17(buffer, 259, 0, 123, 132, 189, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_7(buffer, 274, 0, 132, 141, 202, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_8(buffer, 299, 0, 141, 150, 218, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_16(buffer, 320, 0, 150, 159, 231, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_9(buffer, 338, 0, 177, 189, 274, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_10(buffer, 356, 0, 189, 202, 299, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_11(buffer, 386, 0, 202, 218, 320, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_1(buffer, 413, 0, 244, 259, 338, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_2(buffer, 431, 0, 259, 274, 356, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_3(buffer, 449, 0, 274, 299, 386, ncols, alpha, beta, p);

            compute_prim_pp_electron_repulsion_2(buffer, 487, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 490, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 493, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 496, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 499, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 502, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 505, 3, 15, 39, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 508, 3, 21, 51, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 511, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 514, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 517, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 520, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 523, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 526, 3, 39, 69, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 529, 3, 51, 78, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 532, 3, 54, 84, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 535, 3, 57, 90, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 538, 3, 60, 96, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 541, 3, 63, 102, ncols, p);

            compute_prim_fp_electron_repulsion_10(buffer, 544, 3, 66, 108, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 547, 3, 78, 132, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 550, 3, 84, 141, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 553, 3, 90, 150, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 556, 3, 96, 159, ncols, p);

            compute_prim_gp_electron_repulsion_11(buffer, 559, 3, 102, 168, ncols, p);

            compute_prim_hp_electron_repulsion_11(buffer, 562, 3, 132, 189, ncols, p);

            compute_prim_hp_electron_repulsion_18(buffer, 565, 3, 141, 202, ncols, p);

            compute_prim_hp_electron_repulsion_11(buffer, 568, 3, 150, 218, ncols, p);

            compute_prim_hp_electron_repulsion_11(buffer, 571, 3, 159, 231, ncols, p);

            compute_prim_ip_electron_repulsion_12(buffer, 574, 3, 189, 274, ncols, p);

            compute_prim_ip_electron_repulsion_13(buffer, 577, 0, 3, 202, 568, 299, ncols, p);

            compute_prim_ip_electron_repulsion_14(buffer, 605, 3, 218, 320, ncols, p);

            compute_prim_kp_electron_repulsion_6(buffer, 614, 0, 3, 274, 577, 356, ncols, p);

            compute_prim_kp_electron_repulsion_7(buffer, 682, 0, 3, 299, 605, 386, ncols, p);

            compute_prim_lp_electron_repulsion_1(buffer, 726, 0, 3, 356, 682, 449, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 827, 3, 487, 45, 48, 508, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 830, 3, 490, 48, 51, 511, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 833, 3, 493, 51, 54, 514, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 836, 3, 496, 54, 57, 517, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 839, 3, 499, 57, 60, 520, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 842, 3, 502, 60, 63, 523, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_8(buffer, 845, 3, 505, 63, 66, 526, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 848, 0, 3, 511, 833, 72, 78, 532, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 857, 0, 3, 514, 836, 78, 84, 535, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 866, 0, 3, 517, 839, 84, 90, 538, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 875, 0, 3, 520, 842, 90, 96, 541, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_13(buffer, 884, 0, 3, 523, 845, 96, 102, 544, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 893, 0, 3, 827, 830, 529, 848, 114, 123, 547, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 908, 0, 3, 830, 833, 532, 857, 123, 132, 550, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 923, 0, 3, 833, 836, 535, 866, 132, 141, 553, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 938, 0, 3, 836, 839, 538, 875, 141, 150, 556, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_20(buffer, 953, 0, 3, 839, 842, 541, 884, 150, 159, 559, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_35(buffer, 968, 0, 3, 848, 857, 550, 923, 177, 189, 565, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_36(buffer, 992, 0, 3, 857, 866, 553, 938, 189, 202, 568, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_37(buffer, 1016, 0, 3, 866, 875, 556, 953, 202, 218, 571, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_22(buffer, 1040, 0, 3, 893, 908, 562, 968, 244, 259, 574, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_23(buffer, 1076, 0, 3, 908, 923, 565, 992, 259, 274, 577, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_24(buffer, 1112, 0, 3, 923, 938, 568, 1016, 274, 299, 605, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_9(buffer, 1154, 0, 3, 968, 992, 577, 1112, 338, 356, 682, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 1247, 0, 3, 1040, 1076, 614, 1154, 413, 431, 726, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1517, 1247, 270, ncols);
        }
    }

    simdtrf::transform_ld(values, nvalues, buffer, 1517, nmax);
}

}  // namespace simdt2ceri
