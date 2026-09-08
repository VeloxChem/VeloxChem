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


#include "SimdElectronRepulsionRecDI.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformDI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_di_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_di_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(986, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 58, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 61, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 64, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 67, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 70, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 73, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 76, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 79, 3, 10, 22, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 82, 3, 11, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 85, 3, 12, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 88, 3, 13, 31, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 91, 3, 19, 43, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 100, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 109, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 118, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 127, 3, 31, 55, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 136, 3, 9, 10, 61, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 139, 3, 10, 11, 64, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 142, 3, 11, 12, 67, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 145, 3, 12, 13, 70, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 148, 3, 13, 14, 73, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 151, 0, 61, 139, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 154, 0, 64, 142, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 157, 0, 67, 145, ncols, p);

            compute_prim_dd_electron_repulsion_2(buffer, 160, 3, 76, 37, 40, 91, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 169, 3, 79, 40, 43, 100, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 178, 3, 82, 43, 46, 109, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 187, 3, 85, 46, 49, 118, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 196, 3, 88, 49, 52, 127, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 205, 3, 58, 61, 139, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_3(buffer, 211, 3, 61, 64, 142, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 220, 3, 64, 67, 145, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 226, 3, 67, 70, 148, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_8(buffer, 232, 0, 136, 205, ncols, p);

            compute_prim_pf_electron_repulsion_9(buffer, 235, 0, 139, 211, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 238, 0, 142, 220, ncols, p);

            compute_prim_df_electron_repulsion_2(buffer, 241, 3, 151, 91, 100, 178, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 259, 3, 154, 100, 109, 187, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 277, 3, 157, 109, 118, 196, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_3(buffer, 295, 3, 136, 139, 211, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_4(buffer, 307, 3, 139, 142, 220, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 319, 3, 142, 145, 226, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_6(buffer, 328, 0, 3, 211, 307, 238, ncols, p);

            compute_prim_pg_electron_repulsion_7(buffer, 344, 0, 220, 319, ncols, p);

            compute_prim_dg_electron_repulsion_2(buffer, 350, 3, 232, 160, 169, 241, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_2(buffer, 377, 3, 235, 169, 178, 259, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_3(buffer, 404, 0, 3, 238, 344, 178, 187, 277, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_5(buffer, 437, 3, 205, 211, 307, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_6(buffer, 454, 3, 211, 220, 319, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_3(buffer, 467, 0, 3, 295, 437, 328, ncols, p);

            compute_prim_ph_electron_repulsion_4(buffer, 503, 0, 3, 307, 454, 344, ncols, p);

            compute_prim_dh_electron_repulsion_1(buffer, 527, 0, 3, 328, 503, 241, 259, 404, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 605, 3, 295, 307, 454, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_1(buffer, 618, 0, 3, 437, 605, 503, ncols, p);

            compute_prim_di_electron_repulsion_0(buffer, 650, 0, 3, 467, 618, 350, 377, 527, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 818, 650, 168, ncols);
        }
    }

    simdtrf::transform_di(values, nvalues, buffer, 818, nmax);
}

}  // namespace simdt2ceri
