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


#include "SimdElectronRepulsionRecDL.hpp"

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
#include "SimdElectronRepulsionVrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDK.hpp"
#include "SimdElectronRepulsionVrrRecDL.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPL.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSL.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformDL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_dl_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_dl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1955, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 72, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 75, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 78, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 81, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 84, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 87, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 90, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 93, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 96, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 99, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 102, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 105, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 108, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 111, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 114, 3, 15, 39, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 117, 3, 21, 51, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 126, 3, 24, 54, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 135, 3, 27, 57, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 144, 3, 30, 60, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 153, 3, 33, 63, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 162, 3, 36, 66, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 171, 3, 39, 69, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 180, 3, 9, 10, 75, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 183, 3, 10, 11, 78, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 186, 3, 11, 12, 81, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 189, 3, 12, 13, 84, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 192, 3, 13, 14, 87, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 195, 3, 14, 15, 90, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 198, 3, 15, 16, 93, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 201, 0, 75, 183, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 204, 0, 78, 186, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 207, 0, 81, 189, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 210, 0, 84, 192, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 213, 0, 87, 195, ncols, p);

            compute_prim_dd_electron_repulsion_2(buffer, 216, 3, 96, 45, 48, 117, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 225, 3, 99, 48, 51, 126, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 234, 3, 102, 51, 54, 135, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 243, 3, 105, 54, 57, 144, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 252, 3, 108, 57, 60, 153, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 261, 3, 111, 60, 63, 162, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 270, 3, 114, 63, 66, 171, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 279, 3, 72, 75, 183, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 285, 3, 75, 78, 186, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 291, 3, 78, 81, 189, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 297, 3, 81, 84, 192, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 303, 3, 84, 87, 195, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 309, 3, 87, 90, 198, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_8(buffer, 315, 0, 180, 279, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 318, 0, 183, 285, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 321, 0, 186, 291, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 324, 0, 189, 297, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 327, 0, 192, 303, ncols, p);

            compute_prim_df_electron_repulsion_2(buffer, 330, 3, 201, 117, 126, 234, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 348, 3, 204, 126, 135, 243, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 366, 3, 207, 135, 144, 252, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 384, 3, 210, 144, 153, 261, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 402, 3, 213, 153, 162, 270, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 420, 3, 180, 183, 285, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 429, 3, 183, 186, 291, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 438, 3, 186, 189, 297, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 447, 3, 189, 192, 303, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 456, 3, 192, 195, 309, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_8(buffer, 465, 0, 285, 429, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 468, 0, 291, 438, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 471, 0, 297, 447, ncols, p);

            compute_prim_dg_electron_repulsion_2(buffer, 474, 3, 315, 216, 225, 330, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 501, 3, 318, 225, 234, 348, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 528, 3, 321, 234, 243, 366, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 555, 3, 324, 243, 252, 384, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 582, 3, 327, 252, 261, 402, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 609, 3, 279, 285, 429, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_4(buffer, 622, 3, 285, 291, 438, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 638, 3, 291, 297, 447, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 651, 3, 297, 303, 456, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_8(buffer, 664, 0, 420, 609, ncols, p);

            compute_prim_ph_electron_repulsion_9(buffer, 667, 0, 429, 622, ncols, p);

            compute_prim_ph_electron_repulsion_8(buffer, 670, 0, 438, 638, ncols, p);

            compute_prim_dh_electron_repulsion_2(buffer, 673, 3, 465, 330, 348, 528, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 712, 3, 468, 348, 366, 555, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 751, 3, 471, 366, 384, 582, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_3(buffer, 790, 3, 420, 429, 622, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_4(buffer, 815, 3, 429, 438, 638, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 836, 3, 438, 447, 651, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_6(buffer, 854, 0, 3, 622, 815, 670, ncols, p);

            compute_prim_pi_electron_repulsion_7(buffer, 878, 0, 638, 836, ncols, p);

            compute_prim_di_electron_repulsion_2(buffer, 884, 3, 664, 474, 501, 673, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_2(buffer, 938, 3, 667, 501, 528, 712, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_3(buffer, 992, 0, 3, 670, 878, 528, 555, 751, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_5(buffer, 1052, 3, 609, 622, 815, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_6(buffer, 1082, 3, 622, 638, 836, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_3(buffer, 1106, 0, 3, 790, 1052, 854, ncols, p);

            compute_prim_pk_electron_repulsion_4(buffer, 1172, 0, 3, 815, 1082, 878, ncols, p);

            compute_prim_dk_electron_repulsion_1(buffer, 1208, 0, 3, 854, 1172, 673, 712, 992, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_1(buffer, 1337, 3, 790, 815, 1082, ncols, alpha, beta, p);

            compute_prim_pl_electron_repulsion_1(buffer, 1361, 0, 3, 1052, 1337, 1172, ncols, p);

            compute_prim_dl_electron_repulsion_0(buffer, 1415, 0, 3, 1106, 1361, 884, 938, 1208, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1685, 1415, 270, ncols);
        }
    }

    simdtrf::transform_dl(values, nvalues, buffer, 1685, nmax);
}

}  // namespace simdt2ceri
