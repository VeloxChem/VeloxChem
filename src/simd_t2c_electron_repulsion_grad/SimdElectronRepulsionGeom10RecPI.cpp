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


#include "SimdElectronRepulsionGeom10RecPI.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_pi_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_pi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 2552, 2261, 252, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 79, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 82, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 85, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 88, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 91, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 94, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 97, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 100, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 109, 3, 10, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 118, 3, 11, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 127, 3, 12, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 136, 3, 13, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 145, 3, 14, 34, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 154, 0, 3, 19, 109, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 172, 0, 3, 22, 118, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 190, 0, 3, 25, 127, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 208, 0, 3, 28, 136, 67, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 226, 0, 3, 31, 145, 73, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 244, 3, 7, 8, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 250, 3, 8, 9, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 256, 3, 9, 10, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 262, 3, 10, 11, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 268, 3, 11, 12, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 274, 3, 12, 13, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 280, 3, 13, 14, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 286, 0, 3, 82, 256, 109, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 304, 0, 3, 85, 262, 118, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 322, 0, 3, 88, 268, 127, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 340, 0, 3, 91, 274, 136, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 358, 0, 3, 94, 280, 145, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 376, 0, 3, 100, 286, 37, 43, 154, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 412, 0, 3, 109, 304, 43, 49, 172, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 448, 0, 3, 118, 322, 49, 55, 190, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 484, 0, 3, 127, 340, 55, 61, 208, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 520, 0, 3, 136, 358, 61, 67, 226, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 556, 3, 79, 82, 256, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 566, 3, 82, 85, 262, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 576, 3, 85, 88, 268, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 586, 3, 88, 91, 274, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 596, 3, 91, 94, 280, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 606, 0, 3, 256, 566, 304, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 636, 0, 3, 262, 576, 322, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 666, 0, 3, 268, 586, 340, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 696, 0, 3, 274, 596, 358, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 726, 0, 3, 304, 636, 154, 172, 448,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 786, 0, 3, 322, 666, 172, 190, 484,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 846, 0, 3, 340, 696, 190, 208, 520,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 906, 3, 244, 250, 556, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 921, 3, 250, 256, 566, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 936, 3, 256, 262, 576, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 951, 3, 262, 268, 586, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 966, 3, 268, 274, 596, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 981, 0, 3, 566, 936, 636, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1026, 0, 3, 576, 951, 666, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1071, 0, 3, 586, 966, 696, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1116, 0, 3, 606, 981, 376, 412, 726,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1206, 0, 3, 636, 1026, 412, 448, 786,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1296, 0, 3, 666, 1071, 448, 484, 846,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1386, 3, 556, 566, 936, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1407, 3, 566, 576, 951, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 1428, 3, 576, 586, 966, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 1449, 0, 3, 936, 1407, 1026, ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 1512, 0, 3, 951, 1428, 1071, ncols, p);

            compute_prim_dh_electron_repulsion_0(buffer, 1575, 0, 3, 1026, 1512, 726, 786, 1296,
                                                 ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 1701, 3, 906, 921, 1386, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 1729, 3, 936, 951, 1428, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 1757, 0, 3, 1407, 1729, 1512, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 1841, 0, 3, 1449, 1757, 1116, 1206,
                                                 1575, ncols, alpha, beta, p);

            simdgeo::geom_p_x(buffer, 2009, 1701, 1841, 1, 28, ncols, alpha);

            simdgeo::geom_p_y(buffer, 2093, 1701, 1841, 1, 28, ncols, alpha);

            simdgeo::geom_p_z(buffer, 2177, 1701, 1841, 1, 28, ncols, alpha);

            simdfunc::contract_primitives(buffer, 2261, 2009, 252, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 2513, 2261, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 2513, 13, nmax);

    simdtrf::transform_i_inner(buffer, 2513, 2345, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 39 * nvalues, nvalues, buffer, 2513, 13, nmax);

    simdtrf::transform_i_inner(buffer, 2513, 2429, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 78 * nvalues, nvalues, buffer, 2513, 13, nmax);
}

}  // namespace simdt2ceri
