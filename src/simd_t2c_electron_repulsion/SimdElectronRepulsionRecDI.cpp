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
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_di_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2151, nvalues);

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

            compute_prim_ds_electron_repulsion_0(buffer, 37, 0, 7, 8, 16, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 43, 0, 8, 9, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 49, 0, 9, 10, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 55, 0, 10, 11, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 61, 0, 11, 12, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 12, 13, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 13, 14, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 79, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 82, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 85, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 88, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 91, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 94, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 97, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 106, 3, 10, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 115, 3, 11, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 124, 3, 12, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 133, 3, 13, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 142, 3, 14, 34, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 151, 0, 3, 19, 106, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 169, 0, 3, 22, 115, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 187, 0, 3, 25, 124, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 205, 0, 3, 28, 133, 67, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 223, 0, 3, 31, 142, 73, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 241, 3, 9, 10, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 247, 3, 10, 11, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 253, 3, 11, 12, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 259, 3, 12, 13, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 265, 3, 13, 14, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 271, 0, 3, 79, 241, 106, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 289, 0, 3, 82, 247, 115, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 307, 0, 3, 85, 253, 124, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 325, 0, 3, 88, 259, 133, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 343, 0, 3, 91, 265, 142, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 361, 0, 3, 97, 271, 37, 43, 151, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 397, 0, 3, 106, 289, 43, 49, 169, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 433, 0, 3, 115, 307, 49, 55, 187, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 469, 0, 3, 124, 325, 55, 61, 205, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 505, 0, 3, 133, 343, 61, 67, 223, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 541, 3, 79, 82, 247, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 551, 3, 82, 85, 253, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 561, 3, 85, 88, 259, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 571, 3, 88, 91, 265, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 581, 0, 3, 241, 541, 289, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 611, 0, 3, 247, 551, 307, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 641, 0, 3, 253, 561, 325, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 671, 0, 3, 259, 571, 343, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 701, 0, 3, 289, 611, 151, 169, 433,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 761, 0, 3, 307, 641, 169, 187, 469,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 821, 0, 3, 325, 671, 187, 205, 505,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 881, 3, 241, 247, 551, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 896, 3, 247, 253, 561, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 911, 3, 253, 259, 571, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 926, 0, 3, 541, 881, 611, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 971, 0, 3, 551, 896, 641, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1016, 0, 3, 561, 911, 671, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1061, 0, 3, 581, 926, 361, 397, 701,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1151, 0, 3, 611, 971, 397, 433, 761,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1241, 0, 3, 641, 1016, 433, 469, 821,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1331, 3, 541, 551, 896, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1352, 3, 551, 561, 911, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 1373, 0, 3, 881, 1331, 971, ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 1436, 0, 3, 896, 1352, 1016, ncols, p);

            compute_prim_dh_electron_repulsion_0(buffer, 1499, 0, 3, 971, 1436, 701, 761, 1241,
                                                 ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 1625, 3, 881, 896, 1352, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 1653, 0, 3, 1331, 1625, 1436, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 1737, 0, 3, 1373, 1653, 1061, 1151,
                                                 1499, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1905, 1737, 168, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 2073, 1905, 6, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 2073, 13, nmax);
}

}  // namespace simdt2ceri
