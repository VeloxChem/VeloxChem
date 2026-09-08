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


#include "SimdElectronRepulsionRecGH.hpp"

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
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformGH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gh_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1950, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 40, 0, 7, 8, 19, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 43, 0, 8, 9, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 46, 0, 9, 10, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 49, 0, 10, 11, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 11, 12, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 12, 13, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 58, 0, 13, 14, 37, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 61, 0, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 64, 0, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 70, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 76, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 82, 0, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 88, 0, 31, 34, 58, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 94, 0, 40, 43, 64, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 100, 0, 43, 46, 70, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 106, 0, 46, 49, 76, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 112, 0, 49, 52, 82, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 118, 0, 52, 55, 88, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 124, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 127, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 130, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 133, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 136, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 139, 3, 14, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 142, 3, 9, 22, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 145, 3, 10, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 148, 3, 11, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 151, 3, 12, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 154, 3, 13, 34, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 157, 3, 16, 40, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 160, 3, 19, 43, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 163, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 172, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 181, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 190, 3, 31, 55, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 199, 3, 34, 58, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 208, 3, 43, 64, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 217, 3, 46, 70, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 226, 3, 49, 76, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 235, 3, 52, 82, ncols, p);

            compute_prim_fp_electron_repulsion_7(buffer, 244, 3, 55, 88, ncols, p);

            compute_prim_gp_electron_repulsion_2(buffer, 253, 3, 61, 94, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 265, 3, 64, 100, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 277, 3, 70, 106, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 289, 3, 76, 112, ncols, p);

            compute_prim_gp_electron_repulsion_3(buffer, 301, 3, 82, 118, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 313, 3, 8, 9, 127, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 316, 3, 9, 10, 130, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 319, 3, 10, 11, 133, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 322, 3, 11, 12, 136, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 325, 3, 12, 13, 139, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 328, 0, 124, 313, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 331, 0, 127, 316, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 334, 0, 130, 319, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 337, 0, 133, 322, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 340, 0, 136, 325, ncols, p);

            compute_prim_dd_electron_repulsion_11(buffer, 343, 3, 142, 40, 43, 163, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 346, 3, 145, 43, 46, 172, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 355, 3, 148, 46, 49, 181, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 364, 3, 151, 49, 52, 190, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 373, 3, 154, 52, 55, 199, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_10(buffer, 382, 0, 3, 163, 346, 61, 64, 217, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 394, 0, 3, 172, 355, 64, 70, 226, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_12(buffer, 409, 0, 3, 181, 364, 70, 76, 235, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_11(buffer, 427, 0, 3, 190, 373, 76, 82, 244, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_5(buffer, 442, 0, 3, 343, 346, 217, 394, 94, 100, 277, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_6(buffer, 460, 0, 3, 346, 355, 226, 409, 100, 106, 289, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_7(buffer, 478, 0, 3, 355, 364, 235, 427, 106, 112, 301, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 496, 3, 124, 127, 316, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_2(buffer, 499, 3, 127, 130, 319, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_4(buffer, 502, 3, 130, 133, 322, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 508, 3, 133, 136, 325, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_10(buffer, 514, 0, 316, 499, ncols, p);

            compute_prim_pf_electron_repulsion_15(buffer, 517, 0, 3, 319, 502, 337, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 526, 0, 322, 508, ncols, p);

            compute_prim_df_electron_repulsion_9(buffer, 529, 3, 328, 157, 160, 343, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_14(buffer, 532, 3, 331, 160, 163, 346, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_15(buffer, 535, 0, 3, 334, 517, 163, 172, 355, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_16(buffer, 559, 0, 3, 337, 526, 172, 181, 364, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 583, 3, 340, 181, 190, 373, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_7(buffer, 601, 0, 3, 346, 535, 208, 217, 394, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_8(buffer, 625, 0, 3, 355, 559, 217, 226, 409, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_9(buffer, 667, 0, 3, 364, 583, 226, 235, 427, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_2(buffer, 697, 0, 3, 529, 532, 382, 601, 253, 265, 442, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_3(buffer, 733, 0, 3, 532, 535, 394, 625, 265, 277, 460, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_4(buffer, 769, 0, 3, 535, 559, 409, 667, 277, 289, 478, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_7(buffer, 823, 3, 313, 316, 499, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_8(buffer, 826, 3, 316, 319, 502, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_2(buffer, 829, 3, 319, 322, 508, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_17(buffer, 835, 0, 496, 823, ncols, p);

            compute_prim_pg_electron_repulsion_17(buffer, 838, 0, 499, 826, ncols, p);

            compute_prim_pg_electron_repulsion_10(buffer, 841, 0, 502, 829, ncols, p);

            compute_prim_dg_electron_repulsion_12(buffer, 847, 3, 514, 343, 346, 535, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_13(buffer, 850, 0, 3, 517, 841, 346, 355, 559, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 892, 3, 526, 355, 364, 583, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_5(buffer, 919, 0, 3, 535, 850, 382, 394, 625, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_6(buffer, 1011, 0, 3, 559, 892, 394, 409, 667, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_1(buffer, 1068, 0, 3, 847, 850, 625, 1011, 442, 460, 769, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_9(buffer, 1212, 3, 835, 529, 532, 847, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_10(buffer, 1215, 3, 838, 532, 535, 850, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_11(buffer, 1218, 3, 841, 535, 559, 892, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_4(buffer, 1245, 0, 3, 850, 1218, 601, 625, 1011, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 1320, 0, 3, 1212, 1215, 919, 1245, 697, 733, 1068, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1635, 1320, 315, ncols);
        }
    }

    simdtrf::transform_gh(values, nvalues, buffer, 1635, nmax);
}

}  // namespace simdt2ceri
