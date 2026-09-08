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


#include "SimdElectronRepulsionRecFH.hpp"

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
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformFH.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fh_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1170, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 15, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 36, 0, 7, 8, 18, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 39, 0, 8, 9, 21, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 42, 0, 9, 10, 24, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 45, 0, 10, 11, 27, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 48, 0, 11, 12, 30, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 51, 0, 12, 13, 33, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 54, 0, 15, 18, 39, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 57, 0, 18, 21, 42, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 60, 0, 21, 24, 45, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 63, 0, 24, 27, 48, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 66, 0, 27, 30, 51, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 69, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 72, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 75, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 78, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 81, 3, 14, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 84, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 87, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 90, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 93, 3, 12, 30, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 96, 3, 18, 39, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 105, 3, 21, 42, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 114, 3, 24, 45, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 123, 3, 27, 48, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 132, 3, 30, 51, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 141, 3, 36, 54, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 150, 3, 39, 57, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 159, 3, 42, 60, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 168, 3, 45, 63, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 177, 3, 48, 66, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 186, 3, 9, 10, 72, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 189, 3, 10, 11, 75, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 192, 3, 11, 12, 78, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 195, 3, 12, 13, 81, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 198, 0, 69, 186, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 201, 0, 72, 189, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 204, 0, 75, 192, ncols, p);

            compute_prim_dd_electron_repulsion_2(buffer, 207, 3, 84, 36, 39, 105, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 216, 3, 87, 39, 42, 114, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 225, 3, 90, 42, 45, 123, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 234, 3, 93, 45, 48, 132, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 243, 3, 105, 54, 57, 159, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 252, 3, 114, 57, 60, 168, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 261, 3, 123, 60, 63, 177, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_3(buffer, 270, 3, 69, 72, 189, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_3(buffer, 279, 3, 72, 75, 192, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 288, 3, 75, 78, 195, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_12(buffer, 294, 0, 3, 186, 270, 201, ncols, p);

            compute_prim_pf_electron_repulsion_12(buffer, 303, 0, 3, 189, 279, 204, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 312, 0, 192, 288, ncols, p);

            compute_prim_df_electron_repulsion_2(buffer, 315, 3, 198, 96, 105, 216, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_7(buffer, 333, 0, 3, 201, 303, 105, 114, 225, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_8(buffer, 360, 0, 3, 204, 312, 114, 123, 234, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 381, 3, 207, 141, 150, 243, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 399, 3, 216, 150, 159, 252, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_3(buffer, 417, 0, 3, 225, 360, 159, 168, 261, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_5(buffer, 450, 3, 186, 189, 279, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 459, 3, 189, 192, 288, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_9(buffer, 468, 0, 3, 270, 450, 303, ncols, p);

            compute_prim_pg_electron_repulsion_10(buffer, 489, 0, 279, 459, ncols, p);

            compute_prim_dg_electron_repulsion_5(buffer, 495, 0, 3, 294, 468, 207, 216, 333, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_6(buffer, 549, 0, 3, 303, 489, 216, 225, 360, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_1(buffer, 591, 0, 3, 333, 549, 243, 252, 417, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_4(buffer, 684, 3, 270, 279, 459, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_5(buffer, 693, 0, 450, 684, ncols, p);

            compute_prim_dh_electron_repulsion_2(buffer, 702, 0, 3, 468, 693, 315, 333, 549, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 750, 0, 3, 495, 702, 381, 399, 591, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 960, 750, 210, ncols);
        }
    }

    simdtrf::transform_fh(values, nvalues, buffer, 960, nmax);
}

}  // namespace simdt2ceri
