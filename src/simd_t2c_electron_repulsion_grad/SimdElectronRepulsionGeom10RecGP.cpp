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


#include "SimdElectronRepulsionGeom10RecGP.hpp"

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
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_gp_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_gp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 670, 490, 135, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6}, ncols,
                                            fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 13, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 12, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 28, 0, 7, 8, 16, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 34, 0, 8, 9, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 40, 0, 9, 10, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 46, 0, 10, 11, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 52, 0, 13, 16, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 62, 0, 16, 19, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 72, 0, 19, 22, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 82, 0, 28, 34, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 97, 0, 34, 40, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 112, 0, 52, 62, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_pp_electron_repulsion_0(buffer, 133, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 142, 3, 11, 25, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 151, 0, 3, 16, 133, 34, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 169, 0, 3, 22, 142, 46, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 187, 0, 3, 28, 151, 52, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 217, 0, 3, 40, 169, 72, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 247, 0, 3, 62, 217, 97, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 292, 0, 3, 82, 247, 112, ncols, p);

            simdgeo::geom_g_x(buffer, 355, 187, 292, 1, 3, ncols, alpha);

            simdgeo::geom_g_y(buffer, 400, 187, 292, 1, 3, ncols, alpha);

            simdgeo::geom_g_z(buffer, 445, 187, 292, 1, 3, ncols, alpha);

            simdfunc::contract_primitives(buffer, 490, 355, 135, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 625, 490, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 625, 3, nmax);

    simdtrf::transform_p_inner(buffer, 625, 535, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 27 * nvalues, nvalues, buffer, 625, 3, nmax);

    simdtrf::transform_p_inner(buffer, 625, 580, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 54 * nvalues, nvalues, buffer, 625, 3, nmax);
}

}  // namespace simdt2ceri
