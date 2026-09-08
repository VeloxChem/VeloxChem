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


#include "SimdElectronRepulsionRecDK.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformDK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_dk_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_dk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1407, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 61, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 64, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 67, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 70, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 73, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 76, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 79, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 82, 3, 9, 22, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 85, 3, 10, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 88, 3, 11, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 91, 3, 12, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 94, 3, 13, 34, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 97, 3, 16, 40, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 106, 3, 19, 43, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 115, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 124, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 133, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 142, 3, 31, 55, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 151, 3, 34, 58, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 160, 3, 8, 9, 64, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 163, 3, 9, 10, 67, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 166, 3, 10, 11, 70, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 169, 3, 11, 12, 73, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 172, 3, 12, 13, 76, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 175, 3, 13, 14, 79, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 178, 0, 61, 160, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 181, 0, 64, 163, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 184, 0, 67, 166, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 187, 0, 70, 169, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 190, 0, 73, 172, ncols, p);

            compute_prim_dd_electron_repulsion_2(buffer, 193, 3, 82, 40, 43, 115, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 202, 3, 85, 43, 46, 124, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 211, 3, 88, 46, 49, 133, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 220, 3, 91, 49, 52, 142, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 229, 3, 94, 52, 55, 151, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 238, 3, 61, 64, 163, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 244, 3, 64, 67, 166, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 250, 3, 67, 70, 169, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 256, 3, 70, 73, 172, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 262, 3, 73, 76, 175, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_8(buffer, 268, 0, 163, 244, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 271, 0, 166, 250, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 274, 0, 169, 256, ncols, p);

            compute_prim_df_electron_repulsion_2(buffer, 277, 3, 178, 97, 106, 193, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 295, 3, 181, 106, 115, 202, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 313, 3, 184, 115, 124, 211, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 331, 3, 187, 124, 133, 220, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 349, 3, 190, 133, 142, 229, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 367, 3, 160, 163, 244, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_4(buffer, 376, 3, 163, 166, 250, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 388, 3, 166, 169, 256, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 397, 3, 169, 172, 262, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_8(buffer, 406, 0, 238, 367, ncols, p);

            compute_prim_pg_electron_repulsion_9(buffer, 409, 0, 244, 376, ncols, p);

            compute_prim_pg_electron_repulsion_8(buffer, 412, 0, 250, 388, ncols, p);

            compute_prim_dg_electron_repulsion_2(buffer, 415, 3, 268, 193, 202, 313, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 442, 3, 271, 202, 211, 331, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 469, 3, 274, 211, 220, 349, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_3(buffer, 496, 3, 238, 244, 376, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_4(buffer, 514, 3, 244, 250, 388, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 530, 3, 250, 256, 397, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_6(buffer, 543, 0, 3, 376, 514, 412, ncols, p);

            compute_prim_ph_electron_repulsion_7(buffer, 563, 0, 388, 530, ncols, p);

            compute_prim_dh_electron_repulsion_2(buffer, 569, 3, 406, 277, 295, 415, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_2(buffer, 608, 3, 409, 295, 313, 442, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_3(buffer, 647, 0, 3, 412, 563, 313, 331, 469, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_6(buffer, 692, 3, 367, 376, 514, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_7(buffer, 715, 3, 376, 388, 530, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_4(buffer, 733, 0, 3, 496, 692, 543, ncols, p);

            compute_prim_pi_electron_repulsion_5(buffer, 783, 0, 3, 514, 715, 563, ncols, p);

            compute_prim_di_electron_repulsion_1(buffer, 813, 0, 3, 543, 783, 415, 442, 647, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_4(buffer, 915, 3, 496, 514, 715, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_2(buffer, 933, 0, 3, 692, 915, 783, ncols, p);

            compute_prim_dk_electron_repulsion_0(buffer, 975, 0, 3, 733, 933, 569, 608, 813, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1191, 975, 216, ncols);
        }
    }

    simdtrf::transform_dk(values, nvalues, buffer, 1191, nmax);
}

}  // namespace simdt2ceri
