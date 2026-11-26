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


#include "SimdElectronRepulsionGeom10RecPK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_pk_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_pk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 3854, 3485, 324, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 82, 3, 7, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 85, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 88, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 91, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 94, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 97, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 100, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 103, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 106, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 109, 3, 8, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 118, 3, 9, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 127, 3, 10, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 136, 3, 11, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 145, 3, 12, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 154, 3, 13, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 163, 3, 14, 37, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 172, 0, 3, 16, 109, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 190, 0, 3, 19, 118, 46, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 208, 0, 3, 22, 127, 52, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 226, 0, 3, 25, 136, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 244, 0, 3, 28, 145, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 262, 0, 3, 31, 154, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 280, 0, 3, 34, 163, 76, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 298, 3, 7, 8, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 304, 3, 8, 9, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 310, 3, 9, 10, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 316, 3, 10, 11, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 322, 3, 11, 12, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 328, 3, 12, 13, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 334, 3, 13, 14, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 340, 0, 3, 88, 304, 118, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 358, 0, 3, 91, 310, 127, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 376, 0, 3, 94, 316, 136, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 394, 0, 3, 97, 322, 145, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 412, 0, 3, 100, 328, 154, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 430, 0, 3, 103, 334, 163, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 448, 0, 3, 118, 358, 40, 46, 208, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 484, 0, 3, 127, 376, 46, 52, 226, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 520, 0, 3, 136, 394, 52, 58, 244, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 556, 0, 3, 145, 412, 58, 64, 262, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 592, 0, 3, 154, 430, 64, 70, 280, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 628, 3, 82, 85, 298, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 638, 3, 85, 88, 304, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 648, 3, 88, 91, 310, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 658, 3, 91, 94, 316, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 668, 3, 94, 97, 322, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 678, 3, 97, 100, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 688, 3, 100, 103, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 698, 0, 3, 304, 648, 358, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 728, 0, 3, 310, 658, 376, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 758, 0, 3, 316, 668, 394, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 788, 0, 3, 322, 678, 412, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 818, 0, 3, 328, 688, 430, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 848, 0, 3, 340, 698, 172, 190, 448,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 908, 0, 3, 358, 728, 190, 208, 484,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 968, 0, 3, 376, 758, 208, 226, 520,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1028, 0, 3, 394, 788, 226, 244, 556,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1088, 0, 3, 412, 818, 244, 262, 592,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1148, 3, 298, 304, 648, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1163, 3, 304, 310, 658, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1178, 3, 310, 316, 668, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1193, 3, 316, 322, 678, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1208, 3, 322, 328, 688, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1223, 0, 3, 648, 1163, 728, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1268, 0, 3, 658, 1178, 758, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1313, 0, 3, 668, 1193, 788, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1358, 0, 3, 678, 1208, 818, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1403, 0, 3, 728, 1268, 448, 484, 968,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1493, 0, 3, 758, 1313, 484, 520, 1028,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1583, 0, 3, 788, 1358, 520, 556, 1088,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1673, 3, 628, 638, 1148, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1694, 3, 638, 648, 1163, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1715, 3, 648, 658, 1178, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1736, 3, 658, 668, 1193, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1757, 3, 668, 678, 1208, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 1778, 0, 3, 1163, 1715, 1268, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 1841, 0, 3, 1178, 1736, 1313, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 1904, 0, 3, 1193, 1757, 1358, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 1967, 0, 3, 1223, 1778, 848, 908, 1403,
                                                 ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2093, 0, 3, 1268, 1841, 908, 968, 1493,
                                                 ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2219, 0, 3, 1313, 1904, 968, 1028, 1583,
                                                 ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2345, 3, 1148, 1163, 1715, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2373, 3, 1163, 1178, 1736, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2401, 3, 1178, 1193, 1757, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 2429, 0, 3, 1715, 2373, 1841, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 2513, 0, 3, 1736, 2401, 1904, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 2597, 0, 3, 1841, 2513, 1403, 1493,
                                                 2219, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 2765, 3, 1673, 1694, 2345, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 2801, 3, 1715, 1736, 2401, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 2837, 0, 3, 2373, 2801, 2513, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 2945, 0, 3, 2429, 2837, 1967, 2093,
                                                 2597, ncols, alpha, beta, p);

            simdgeo::geom_p_x(buffer, 3161, 2765, 2945, 1, 36, ncols, alpha);

            simdgeo::geom_p_y(buffer, 3269, 2765, 2945, 1, 36, ncols, alpha);

            simdgeo::geom_p_z(buffer, 3377, 2765, 2945, 1, 36, ncols, alpha);

            simdfunc::contract_primitives(buffer, 3485, 3161, 324, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 3809, 3485, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 3809, 15, nmax);

    simdtrf::transform_k_inner(buffer, 3809, 3593, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 45 * nvalues, nvalues, buffer, 3809, 15, nmax);

    simdtrf::transform_k_inner(buffer, 3809, 3701, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 90 * nvalues, nvalues, buffer, 3809, 15, nmax);
}

}  // namespace simdt2ceri
