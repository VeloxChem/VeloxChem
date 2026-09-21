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


#include "SimdElectronRepulsionGeom10RsRecPG.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_pg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_pg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1903, 1606, 270, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 6, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 14, 6, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 15, 16, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 16, 17, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 17, 18, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 18, 19, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 19, 20, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 112, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 115, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 118, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 121, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 124, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 127, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 130, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 133, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 136, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 139, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 142, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 151, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 160, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 169, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 178, 3, 17, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 187, 3, 18, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 196, 3, 19, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 205, 3, 20, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 214, 0, 3, 25, 151, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 232, 0, 3, 28, 160, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 250, 0, 3, 31, 169, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 268, 0, 3, 40, 187, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 286, 0, 3, 43, 196, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 304, 0, 3, 46, 205, 106, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 322, 3, 7, 8, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 328, 3, 8, 9, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 334, 3, 9, 10, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 340, 3, 10, 11, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 346, 3, 11, 12, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 352, 3, 15, 16, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 358, 3, 16, 17, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 364, 3, 17, 18, 133, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 370, 3, 18, 19, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 376, 3, 19, 20, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 382, 0, 3, 115, 334, 151, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 400, 0, 3, 118, 340, 160, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 418, 0, 3, 121, 346, 169, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 436, 0, 3, 130, 364, 187, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 454, 0, 3, 133, 370, 196, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 472, 0, 3, 136, 376, 205, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 490, 0, 3, 142, 382, 52, 58, 214, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 526, 0, 3, 151, 400, 58, 64, 232, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 562, 0, 3, 160, 418, 64, 70, 250, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 598, 0, 3, 178, 436, 82, 88, 268, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 634, 0, 3, 187, 454, 88, 94, 286, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 670, 0, 3, 196, 472, 94, 100, 304,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 706, 3, 112, 115, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 716, 3, 115, 118, 340, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 726, 3, 118, 121, 346, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 736, 3, 127, 130, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 746, 3, 130, 133, 370, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 756, 3, 133, 136, 376, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 766, 0, 3, 334, 716, 400, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 796, 0, 3, 340, 726, 418, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 826, 0, 3, 364, 746, 454, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 856, 0, 3, 370, 756, 472, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 886, 0, 3, 400, 796, 214, 232, 562,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 946, 0, 3, 454, 856, 268, 286, 670,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1006, 3, 322, 328, 706, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1021, 3, 334, 340, 726, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1036, 3, 352, 358, 736, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1051, 3, 364, 370, 756, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1066, 0, 3, 716, 1021, 796, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1111, 0, 3, 746, 1051, 856, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1156, 0, 3, 766, 1066, 490, 526, 886,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1246, 0, 3, 826, 1111, 598, 634, 946,
                                                 ncols, alpha, beta, p);

            simdgeo::geom_p_x(buffer, 1336, 1036, 1246, 1, 15, ncols, alpha);

            simdgeo::geom_p_y(buffer, 1381, 1036, 1246, 1, 15, ncols, alpha);

            simdgeo::geom_p_z(buffer, 1426, 1036, 1246, 1, 15, ncols, alpha);

            simdgeo::geom_p_x(buffer, 1471, 1006, 1156, 1, 15, ncols, alpha);

            simdgeo::geom_p_y(buffer, 1516, 1006, 1156, 1, 15, ncols, alpha);

            simdgeo::geom_p_z(buffer, 1561, 1006, 1156, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 1606, 1471, 135, ncols);

            simdfunc::contract_primitives(buffer, 1741, 1336, 135, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 1876, 1741, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 1876, 9, nmax);

    simdtrf::transform_g_inner(buffer, 1876, 1786, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 27 * nvalues, nvalues, buffer, 1876, 9, nmax);

    simdtrf::transform_g_inner(buffer, 1876, 1831, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 54 * nvalues, nvalues, buffer, 1876, 9, nmax);

    simdtrf::transform_g_inner(buffer, 1876, 1606, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 81 * nvalues, nvalues, buffer, 1876, 9, nmax);

    simdtrf::transform_g_inner(buffer, 1876, 1651, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 108 * nvalues, nvalues, buffer, 1876, 9, nmax);

    simdtrf::transform_g_inner(buffer, 1876, 1696, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 135 * nvalues, nvalues, buffer, 1876, 9, nmax);
}

}  // namespace simdt2ceri
