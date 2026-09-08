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


#include "SimdElectronRepulsionRecSK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_sk_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_sk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(320, nvalues);

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

            simdfunc::compute_pb(buffer, coordinates, 0, ncols, fb);

            simdfunc::compute_boys_function(buffer, coordinates, 3, {1, 2, 3, 4, 5, 6, 7}, ncols,
                                            fj, mu);

            compute_prim_sp_electron_repulsion_0(buffer, 11, 0, 4, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 14, 0, 5, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 17, 0, 6, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 20, 0, 7, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 23, 0, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 26, 0, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 29, 0, 10, ncols);

            compute_prim_sd_electron_repulsion_0(buffer, 32, 0, 4, 5, 17, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 38, 0, 5, 6, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 44, 0, 6, 7, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 50, 0, 7, 8, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 56, 0, 8, 9, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 62, 0, 11, 14, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 72, 0, 14, 17, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 82, 0, 17, 20, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 92, 0, 20, 23, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 102, 0, 23, 26, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 112, 0, 32, 38, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 127, 0, 38, 44, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_sg_electron_repulsion_0(buffer, 142, 0, 44, 50, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_sh_electron_repulsion_0(buffer, 157, 0, 62, 72, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_sh_electron_repulsion_0(buffer, 178, 0, 72, 82, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_sh_electron_repulsion_0(buffer, 199, 0, 82, 92, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_si_electron_repulsion_0(buffer, 220, 0, 112, 127, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 248, 0, 157, 178, 220, ncols, alpha,
                                                 beta, p);

            simdfunc::contract_primitives(buffer, 284, 248, 36, ncols);
        }
    }

    simdtrf::transform_k_outer(values, nvalues, buffer, 284, 1, nmax);
}

}  // namespace simdt2ceri
