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


#include "SimdElectronRepulsionGeom10RsRecDD.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_dd_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_dd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1104, 858, 216, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5}, ncols,
                                                fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 12, {1, 2, 3, 4, 5}, ncols, fj,
                                            mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 48, 0, 7, 8, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 54, 0, 8, 9, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 9, 10, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 13, 14, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 72, 0, 14, 15, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 78, 0, 15, 16, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 84, 0, 18, 21, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 94, 0, 21, 24, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 104, 0, 24, 27, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 114, 0, 33, 36, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 124, 0, 36, 39, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 134, 0, 39, 42, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 144, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 147, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 150, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 153, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 156, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 159, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 162, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 165, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 168, 3, 8, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 177, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 186, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 195, 3, 14, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 204, 3, 15, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 213, 3, 16, 45, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 222, 0, 3, 24, 177, 54, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 240, 0, 3, 27, 186, 60, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 258, 0, 3, 39, 204, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 276, 0, 3, 42, 213, 78, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 294, 0, 3, 54, 240, 104, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 324, 0, 3, 72, 276, 134, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 354, 3, 7, 8, 147, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 360, 3, 9, 10, 153, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 366, 3, 13, 14, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 372, 3, 15, 16, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 378, 0, 3, 144, 354, 168, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 396, 0, 3, 150, 360, 186, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 414, 0, 3, 156, 366, 195, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 432, 0, 3, 162, 372, 213, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 450, 0, 3, 177, 396, 48, 54, 240, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 486, 0, 3, 204, 432, 66, 72, 276, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 522, 0, 3, 222, 450, 84, 94, 294, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 582, 0, 3, 258, 486, 114, 124, 324,
                                                 ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 642, 414, 582, 1, 6, ncols, alpha);

            simdgeo::geom_d_y(buffer, 678, 414, 582, 1, 6, ncols, alpha);

            simdgeo::geom_d_z(buffer, 714, 414, 582, 1, 6, ncols, alpha);

            simdgeo::geom_d_x(buffer, 750, 378, 522, 1, 6, ncols, alpha);

            simdgeo::geom_d_y(buffer, 786, 378, 522, 1, 6, ncols, alpha);

            simdgeo::geom_d_z(buffer, 822, 378, 522, 1, 6, ncols, alpha);

            simdfunc::contract_primitives(buffer, 858, 750, 108, ncols);

            simdfunc::contract_primitives(buffer, 966, 642, 108, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 1074, 966, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 1074, 5, nmax);

    simdtrf::transform_d_inner(buffer, 1074, 1002, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 25 * nvalues, nvalues, buffer, 1074, 5, nmax);

    simdtrf::transform_d_inner(buffer, 1074, 1038, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 50 * nvalues, nvalues, buffer, 1074, 5, nmax);

    simdtrf::transform_d_inner(buffer, 1074, 858, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 75 * nvalues, nvalues, buffer, 1074, 5, nmax);

    simdtrf::transform_d_inner(buffer, 1074, 894, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 100 * nvalues, nvalues, buffer, 1074, 5, nmax);

    simdtrf::transform_d_inner(buffer, 1074, 930, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 125 * nvalues, nvalues, buffer, 1074, 5, nmax);
}

}  // namespace simdt2ceri
