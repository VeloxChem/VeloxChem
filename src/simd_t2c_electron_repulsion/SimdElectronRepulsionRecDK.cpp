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
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"

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

    auto buffer = CSimdMatrix(3314, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 40, 0, 7, 8, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 46, 0, 8, 9, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 9, 10, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 10, 11, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 11, 12, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 12, 13, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 13, 14, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 82, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 85, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 88, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 91, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 94, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 97, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 100, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 103, 3, 8, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 112, 3, 9, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 121, 3, 10, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 130, 3, 11, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 139, 3, 12, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 148, 3, 13, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 157, 3, 14, 37, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 166, 0, 3, 16, 103, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 184, 0, 3, 19, 112, 46, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 202, 0, 3, 22, 121, 52, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 220, 0, 3, 25, 130, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 238, 0, 3, 28, 139, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 256, 0, 3, 31, 148, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 274, 0, 3, 34, 157, 76, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 292, 3, 8, 9, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 298, 3, 9, 10, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 304, 3, 10, 11, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 310, 3, 11, 12, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 316, 3, 12, 13, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 322, 3, 13, 14, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 328, 0, 3, 82, 292, 112, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 346, 0, 3, 85, 298, 121, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 364, 0, 3, 88, 304, 130, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 382, 0, 3, 91, 310, 139, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 400, 0, 3, 94, 316, 148, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 418, 0, 3, 97, 322, 157, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 436, 0, 3, 112, 346, 40, 46, 202, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 472, 0, 3, 121, 364, 46, 52, 220, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 508, 0, 3, 130, 382, 52, 58, 238, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 544, 0, 3, 139, 400, 58, 64, 256, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 580, 0, 3, 148, 418, 64, 70, 274, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 616, 3, 82, 85, 298, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 626, 3, 85, 88, 304, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 636, 3, 88, 91, 310, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 646, 3, 91, 94, 316, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 656, 3, 94, 97, 322, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 666, 0, 3, 292, 616, 346, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 696, 0, 3, 298, 626, 364, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 726, 0, 3, 304, 636, 382, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 756, 0, 3, 310, 646, 400, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 786, 0, 3, 316, 656, 418, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 816, 0, 3, 328, 666, 166, 184, 436,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 876, 0, 3, 346, 696, 184, 202, 472,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 936, 0, 3, 364, 726, 202, 220, 508,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 996, 0, 3, 382, 756, 220, 238, 544,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1056, 0, 3, 400, 786, 238, 256, 580,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1116, 3, 292, 298, 626, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1131, 3, 298, 304, 636, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1146, 3, 304, 310, 646, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1161, 3, 310, 316, 656, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1176, 0, 3, 616, 1116, 696, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1221, 0, 3, 626, 1131, 726, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1266, 0, 3, 636, 1146, 756, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1311, 0, 3, 646, 1161, 786, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1356, 0, 3, 696, 1221, 436, 472, 936,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1446, 0, 3, 726, 1266, 472, 508, 996,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1536, 0, 3, 756, 1311, 508, 544, 1056,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1626, 3, 616, 626, 1131, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1647, 3, 626, 636, 1146, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1668, 3, 636, 646, 1161, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 1689, 0, 3, 1116, 1626, 1221, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 1752, 0, 3, 1131, 1647, 1266, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 1815, 0, 3, 1146, 1668, 1311, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 1878, 0, 3, 1176, 1689, 816, 876, 1356,
                                                 ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2004, 0, 3, 1221, 1752, 876, 936, 1446,
                                                 ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2130, 0, 3, 1266, 1815, 936, 996, 1536,
                                                 ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2256, 3, 1116, 1131, 1647, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2284, 3, 1131, 1146, 1668, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 2312, 0, 3, 1626, 2256, 1752, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 2396, 0, 3, 1647, 2284, 1815, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 2480, 0, 3, 1752, 2396, 1356, 1446,
                                                 2130, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 2648, 3, 1626, 1647, 2284, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 2684, 0, 3, 2256, 2648, 2396, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 2792, 0, 3, 2312, 2684, 1878, 2004,
                                                 2480, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3008, 2792, 216, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 3224, 3008, 6, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 3224, 15, nmax);
}

}  // namespace simdt2ceri
