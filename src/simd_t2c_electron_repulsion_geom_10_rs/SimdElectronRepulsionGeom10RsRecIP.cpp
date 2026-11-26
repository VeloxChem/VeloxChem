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


#include "SimdElectronRepulsionGeom10RsRecIP.hpp"

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
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryI1.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ip_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ip_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 2804, 2216, 504, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 15, {1, 2, 3, 4, 5, 6, 7, 8},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 23, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 7, 8, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 8, 9, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 9, 10, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 10, 11, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 11, 12, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 12, 13, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 16, 17, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 17, 18, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 18, 19, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 19, 20, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 20, 21, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 21, 22, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 138, 0, 24, 27, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 148, 0, 27, 30, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 30, 33, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 33, 36, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 36, 39, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 45, 48, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 198, 0, 48, 51, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 208, 0, 51, 54, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 218, 0, 54, 57, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 228, 0, 57, 60, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 238, 0, 66, 72, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 253, 0, 72, 78, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 268, 0, 78, 84, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 283, 0, 84, 90, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 298, 0, 102, 108, 198, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 313, 0, 108, 114, 208, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 328, 0, 114, 120, 218, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 343, 0, 120, 126, 228, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 358, 0, 138, 148, 253, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 379, 0, 148, 158, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 400, 0, 158, 168, 283, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 421, 0, 188, 198, 313, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 442, 0, 198, 208, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 463, 0, 208, 218, 343, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 484, 0, 238, 253, 379, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 512, 0, 253, 268, 400, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 540, 0, 298, 313, 442, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 568, 0, 313, 328, 463, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 596, 0, 358, 379, 512, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 632, 0, 421, 442, 568, ncols, alpha,
                                                 beta, p);

            compute_prim_pp_electron_repulsion_0(buffer, 668, 3, 11, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 677, 3, 13, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 686, 3, 20, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 695, 3, 22, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 704, 0, 3, 33, 668, 84, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 722, 0, 3, 39, 677, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 740, 0, 3, 54, 686, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 758, 0, 3, 60, 695, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 776, 0, 3, 78, 704, 158, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 806, 0, 3, 90, 722, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 836, 0, 3, 114, 740, 208, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 866, 0, 3, 126, 758, 228, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 896, 0, 3, 148, 776, 253, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 941, 0, 3, 168, 806, 283, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 986, 0, 3, 198, 836, 313, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1031, 0, 3, 218, 866, 343, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1076, 0, 3, 238, 896, 358, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1139, 0, 3, 268, 941, 400, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1202, 0, 3, 298, 986, 421, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1265, 0, 3, 328, 1031, 463, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1328, 0, 3, 379, 1139, 512, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 1412, 0, 3, 442, 1265, 568, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 1496, 0, 3, 484, 1328, 596, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 1604, 0, 3, 540, 1412, 632, ncols, p);

            simdgeo::geom_i_x(buffer, 1712, 1202, 1604, 1, 3, ncols, alpha);

            simdgeo::geom_i_y(buffer, 1796, 1202, 1604, 1, 3, ncols, alpha);

            simdgeo::geom_i_z(buffer, 1880, 1202, 1604, 1, 3, ncols, alpha);

            simdgeo::geom_i_x(buffer, 1964, 1076, 1496, 1, 3, ncols, alpha);

            simdgeo::geom_i_y(buffer, 2048, 1076, 1496, 1, 3, ncols, alpha);

            simdgeo::geom_i_z(buffer, 2132, 1076, 1496, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 2216, 1964, 252, ncols);

            simdfunc::contract_primitives(buffer, 2468, 1712, 252, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 2720, 2468, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 2720, 3, nmax);

    simdtrf::transform_p_inner(buffer, 2720, 2552, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 39 * nvalues, nvalues, buffer, 2720, 3, nmax);

    simdtrf::transform_p_inner(buffer, 2720, 2636, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 78 * nvalues, nvalues, buffer, 2720, 3, nmax);

    simdtrf::transform_p_inner(buffer, 2720, 2216, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 117 * nvalues, nvalues, buffer, 2720, 3, nmax);

    simdtrf::transform_p_inner(buffer, 2720, 2300, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 156 * nvalues, nvalues, buffer, 2720, 3, nmax);

    simdtrf::transform_p_inner(buffer, 2720, 2384, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 195 * nvalues, nvalues, buffer, 2720, 3, nmax);
}

}  // namespace simdt2ceri
