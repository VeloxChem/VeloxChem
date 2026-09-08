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


#include "SimdElectronRepulsionRecIP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdTransformIP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ip_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ip_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(355, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 32, 0, 7, 8, 17, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 35, 0, 8, 9, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 38, 0, 9, 10, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 41, 0, 10, 11, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 44, 0, 11, 12, 29, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 47, 0, 14, 17, 35, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 53, 0, 17, 20, 38, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 59, 0, 20, 23, 41, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 65, 0, 23, 26, 44, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 71, 0, 32, 35, 53, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_5(buffer, 80, 0, 35, 38, 59, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 92, 0, 38, 41, 65, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_6(buffer, 101, 0, 47, 53, 80, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_7(buffer, 119, 0, 53, 59, 92, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 135, 0, 71, 80, 119, ncols, alpha, beta, p);

            compute_prim_gp_electron_repulsion_11(buffer, 160, 3, 59, 92, ncols, p);

            compute_prim_hp_electron_repulsion_5(buffer, 163, 0, 3, 80, 160, 119, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 187, 0, 3, 101, 163, 135, ncols, p);

            simdfunc::contract_primitives(buffer, 271, 187, 84, ncols);
        }
    }

    simdtrf::transform_ip(values, nvalues, buffer, 271, nmax);
}

}  // namespace simdt2ceri
