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


#include "SimdElectronRepulsionRsRecPL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

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
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_pl_electron_repulsion(double               *values,
                                 const size_t          nvalues,
                                 const CBasisFunction &bra,
                                 const CBasisFunction &ket,
                                 const CSimdMatrix    &coordinates,
                                 CSimdMatrix          &buffer,
                                 const double          omega) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_pl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 2309, 1988, 270, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8,
                                                9}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 16, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 32, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 35, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 38, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 41, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 44, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 47, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 50, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 53, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 56, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 59, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 62, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 65, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 68, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 71, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 74, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 77, 3, 25, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 80, 3, 14, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 89, 3, 24, 29, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 98, 3, 7, 8, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 104, 3, 8, 9, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 110, 3, 9, 10, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 116, 3, 10, 11, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 122, 3, 11, 12, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 128, 3, 12, 13, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 134, 3, 13, 14, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 140, 3, 17, 18, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 146, 3, 18, 19, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 152, 3, 19, 20, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 158, 3, 20, 21, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 164, 3, 21, 22, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 170, 3, 22, 23, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 176, 3, 23, 24, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 182, 0, 3, 50, 134, 80, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 200, 0, 3, 74, 176, 89, ncols, p);

            compute_prim_sf_electron_repulsion_0(buffer, 218, 3, 32, 35, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 228, 3, 35, 38, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 238, 3, 38, 41, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 248, 3, 41, 44, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 258, 3, 44, 47, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 268, 3, 47, 50, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 278, 3, 56, 59, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 288, 3, 59, 62, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 298, 3, 62, 65, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 308, 3, 65, 68, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 318, 3, 68, 71, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 328, 3, 71, 74, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 338, 0, 3, 128, 268, 182, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 368, 0, 3, 170, 328, 200, ncols, p);

            compute_prim_sg_electron_repulsion_0(buffer, 398, 3, 98, 104, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 413, 3, 104, 110, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 428, 3, 110, 116, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 443, 3, 116, 122, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 458, 3, 122, 128, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 473, 3, 140, 146, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 488, 3, 146, 152, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 503, 3, 152, 158, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 518, 3, 158, 164, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 533, 3, 164, 170, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 548, 0, 3, 258, 458, 338, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 593, 0, 3, 318, 533, 368, ncols, p);

            compute_prim_sh_electron_repulsion_0(buffer, 638, 3, 218, 228, 413, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 659, 3, 228, 238, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 680, 3, 238, 248, 443, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 701, 3, 248, 258, 458, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 722, 3, 278, 288, 488, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 743, 3, 288, 298, 503, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 764, 3, 298, 308, 518, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 785, 3, 308, 318, 533, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 806, 0, 3, 443, 701, 548, ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 869, 0, 3, 518, 785, 593, ncols, p);

            compute_prim_si_electron_repulsion_0(buffer, 932, 3, 398, 413, 659, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 960, 3, 413, 428, 680, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 988, 3, 428, 443, 701, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 1016, 3, 473, 488, 743, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 1044, 3, 488, 503, 764, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 1072, 3, 503, 518, 785, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 1100, 0, 3, 680, 988, 806, ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 1184, 0, 3, 764, 1072, 869, ncols, p);

            compute_prim_sk_electron_repulsion_0(buffer, 1268, 3, 638, 659, 960, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 1304, 3, 659, 680, 988, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 1340, 3, 722, 743, 1044, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 1376, 3, 743, 764, 1072, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 1412, 0, 3, 960, 1304, 1100, ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 1520, 0, 3, 1044, 1376, 1184, ncols,
                                                 p);

            compute_prim_sl_electron_repulsion_0(buffer, 1628, 3, 932, 960, 1304, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 1673, 3, 1016, 1044, 1376, ncols, alpha,
                                                 beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 1718, 0, 3, 1268, 1628, 1412, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 1853, 0, 3, 1340, 1673, 1520, ncols,
                                                 p);

            simdfunc::contract_primitives(buffer, 1988, 1718, 270, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 2258, 2123, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 2258, 17, nmax);

    simdtrf::transform_l_inner(buffer, 2258, 1988, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 51 * nvalues, nvalues, buffer, 2258, 17, nmax);
}

}  // namespace simdt2ceri
