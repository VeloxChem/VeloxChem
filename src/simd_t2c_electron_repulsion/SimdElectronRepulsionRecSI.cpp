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


#include "SimdElectronRepulsionRecSI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionCtrVrrSI.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_si_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_si_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    // NOTE: the values are zeroed before anything writes them, the composed
    // step accumulating into them as the sum over primitives runs. The atom
    // pairs no pair of primitives reaches keep the zeros set here.

    std::fill(values, values + 13 * nvalues, 0.0);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(101, nvalues);

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

            const auto fb = a_exps[i] / p;

            simdfunc::compute_pb(buffer, coordinates, 0, ncols, fb);

            simdfunc::compute_full_boys_function(buffer, coordinates, 3, 6, ncols, fj, mu);

            compute_prim_sp_electron_repulsion_0(buffer, 11, 0, 6, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 14, 0, 7, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_sd_electron_repulsion_1(buffer, 26, 0, 4, 5, 11, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 29, 0, 5, 6, 14, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 32, 0, 6, 7, 17, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 35, 0, 7, 8, 20, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 38, 0, 8, 9, 23, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 41, 0, 11, 14, 32, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 47, 0, 14, 17, 35, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 53, 0, 17, 20, 38, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 59, 0, 26, 29, 41, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 68, 0, 29, 32, 47, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 77, 0, 32, 35, 53, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86, 0, 41, 47, 77, ncols, alpha, beta, p);

            compute_ctr_si_electron_repulsion_0(values, nvalues, buffer, 0, 59, 68, 86, ncols, alpha, beta, p);
        }
    }
}

}  // namespace simdt2ceri
