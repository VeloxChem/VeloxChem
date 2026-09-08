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


#include "SimdElectronRepulsionRecKS.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionCtrVrrKS.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ks_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ks_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    // NOTE: the values are zeroed before anything writes them, the composed
    // step accumulating into them as the sum over primitives runs. The atom
    // pairs no pair of primitives reaches keep the zeros set here.

    std::fill(values, values + 15 * nvalues, 0.0);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(163, nvalues);

    buffer.zero();

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

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_boys_function(buffer, coordinates, 3, {1, 2, 3, 4, 5, 6, 7}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 11, 0, 4, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 5, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 6, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 10, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 32, 0, 4, 5, 17, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 35, 0, 5, 6, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 38, 0, 6, 7, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 41, 0, 7, 8, 26, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 44, 0, 8, 9, 29, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 47, 0, 11, 14, 32, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 53, 0, 14, 17, 35, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 59, 0, 17, 20, 38, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 65, 0, 20, 23, 41, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 71, 0, 23, 26, 44, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 77, 0, 32, 35, 59, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 86, 0, 35, 38, 65, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 95, 0, 38, 41, 71, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 104, 0, 47, 53, 77, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 117, 0, 53, 59, 86, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 130, 0, 59, 65, 95, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_6(buffer, 143, 0, 77, 86, 130, ncols, alpha, beta, p);

            compute_ctr_ks_electron_repulsion_0(values, nvalues, buffer, 0, 104, 117, 143, ncols, alpha, beta, p);
        }
    }
}

}  // namespace simdt2ceri
