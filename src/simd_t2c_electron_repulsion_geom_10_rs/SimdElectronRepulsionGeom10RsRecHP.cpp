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


#include "SimdElectronRepulsionGeom10RsRecHP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_hp_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_hp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1949, 1508, 378, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 14, {1, 2, 3, 4, 5, 6, 7},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 7, 8, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 8, 9, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 9, 10, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 10, 11, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 11, 12, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 15, 16, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 16, 17, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 17, 18, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 18, 19, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 19, 20, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 118, 0, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 128, 0, 25, 28, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 138, 0, 28, 31, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 148, 0, 31, 34, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 40, 43, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 43, 46, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 46, 49, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 49, 52, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 198, 0, 58, 64, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 213, 0, 64, 70, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 228, 0, 70, 76, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 243, 0, 88, 94, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 258, 0, 94, 100, 178, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 273, 0, 100, 106, 188, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 288, 0, 118, 128, 213, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 309, 0, 128, 138, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 330, 0, 158, 168, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 351, 0, 168, 178, 273, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 372, 0, 198, 213, 309, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 400, 0, 243, 258, 351, ncols, alpha,
                                                 beta, p);

            compute_prim_pp_electron_repulsion_0(buffer, 428, 3, 10, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 437, 3, 12, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 446, 3, 18, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 455, 3, 20, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 464, 0, 3, 28, 428, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 482, 0, 3, 34, 437, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 500, 0, 3, 46, 446, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 518, 0, 3, 52, 455, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 536, 0, 3, 64, 464, 128, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 566, 0, 3, 76, 482, 148, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 596, 0, 3, 94, 500, 168, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 626, 0, 3, 106, 518, 188, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 656, 0, 3, 118, 536, 198, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 701, 0, 3, 138, 566, 228, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 746, 0, 3, 158, 596, 243, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 791, 0, 3, 178, 626, 273, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 836, 0, 3, 213, 701, 309, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 899, 0, 3, 258, 791, 351, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 962, 0, 3, 288, 836, 372, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1046, 0, 3, 330, 899, 400, ncols, p);

            simdgeo::geom_h_x(buffer, 1130, 746, 1046, 1, 3, ncols, alpha);

            simdgeo::geom_h_y(buffer, 1193, 746, 1046, 1, 3, ncols, alpha);

            simdgeo::geom_h_z(buffer, 1256, 746, 1046, 1, 3, ncols, alpha);

            simdgeo::geom_h_x(buffer, 1319, 656, 962, 1, 3, ncols, alpha);

            simdgeo::geom_h_y(buffer, 1382, 656, 962, 1, 3, ncols, alpha);

            simdgeo::geom_h_z(buffer, 1445, 656, 962, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 1508, 1319, 189, ncols);

            simdfunc::contract_primitives(buffer, 1697, 1130, 189, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 1886, 1697, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 1886, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1886, 1760, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 33 * nvalues, nvalues, buffer, 1886, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1886, 1823, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 66 * nvalues, nvalues, buffer, 1886, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1886, 1508, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 99 * nvalues, nvalues, buffer, 1886, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1886, 1571, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 132 * nvalues, nvalues, buffer, 1886, 3, nmax);

    simdtrf::transform_p_inner(buffer, 1886, 1634, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 165 * nvalues, nvalues, buffer, 1886, 3, nmax);
}

}  // namespace simdt2ceri
