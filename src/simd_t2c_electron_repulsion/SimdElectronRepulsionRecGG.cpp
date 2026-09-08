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


#include "SimdElectronRepulsionRecGG.hpp"

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
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformGG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gg_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1254, nvalues);

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

            compute_prim_fs_electron_repulsion_4(buffer, 70, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 76, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 82, 0, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 88, 0, 37, 40, 58, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 94, 0, 40, 43, 64, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 100, 0, 43, 46, 70, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 106, 0, 46, 49, 76, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 112, 0, 49, 52, 82, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 118, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 121, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 124, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 127, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 130, 3, 14, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 133, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 136, 3, 10, 22, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 139, 3, 11, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 142, 3, 12, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 145, 3, 13, 31, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 148, 3, 19, 43, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 151, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 160, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 169, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 178, 3, 31, 55, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 187, 3, 43, 64, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 196, 3, 46, 70, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 205, 3, 49, 76, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 214, 3, 52, 82, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 223, 3, 64, 100, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 235, 3, 70, 106, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 247, 3, 76, 112, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 259, 3, 9, 10, 121, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 262, 3, 10, 11, 124, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 265, 3, 11, 12, 127, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 268, 3, 12, 13, 130, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 271, 0, 121, 262, ncols, p);

            compute_prim_pd_electron_repulsion_3(buffer, 274, 0, 124, 265, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 280, 0, 127, 268, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 283, 3, 133, 37, 40, 148, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_11(buffer, 286, 3, 136, 40, 43, 151, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_12(buffer, 289, 0, 3, 139, 274, 43, 46, 160, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_7(buffer, 301, 0, 3, 142, 280, 46, 49, 169, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 313, 3, 145, 49, 52, 178, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_7(buffer, 322, 0, 3, 151, 289, 58, 64, 196, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_8(buffer, 337, 0, 3, 160, 301, 64, 70, 205, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_9(buffer, 361, 0, 3, 169, 313, 70, 76, 214, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_2(buffer, 382, 0, 3, 283, 286, 187, 322, 88, 94, 223, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_3(buffer, 400, 0, 3, 286, 289, 196, 337, 94, 100, 235, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_4(buffer, 418, 0, 3, 289, 301, 205, 361, 100, 106, 247, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 454, 3, 118, 121, 262, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 457, 3, 121, 124, 265, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 460, 3, 124, 127, 268, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 463, 0, 259, 454, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 466, 0, 262, 457, ncols, p);

            compute_prim_pf_electron_repulsion_10(buffer, 469, 0, 265, 460, ncols, p);

            compute_prim_df_electron_repulsion_12(buffer, 472, 3, 271, 148, 151, 289, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_13(buffer, 475, 0, 3, 274, 469, 151, 160, 301, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 505, 3, 280, 160, 169, 313, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_5(buffer, 523, 0, 3, 289, 475, 187, 196, 337, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_6(buffer, 585, 0, 3, 301, 505, 196, 205, 361, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_1(buffer, 624, 0, 3, 472, 475, 337, 585, 223, 235, 418, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_9(buffer, 723, 3, 463, 283, 286, 472, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_10(buffer, 726, 3, 466, 286, 289, 475, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_11(buffer, 729, 3, 469, 289, 301, 505, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_4(buffer, 747, 0, 3, 475, 729, 322, 337, 585, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 804, 0, 3, 723, 726, 523, 747, 382, 400, 624, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1029, 804, 225, ncols);
        }
    }

    simdtrf::transform_gg_tri(values, nvalues, buffer, 1029, nmax);
}

}  // namespace simdt2ceri
