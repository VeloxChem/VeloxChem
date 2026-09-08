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


#include "SimdElectronRepulsionRecSL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionCtrVrrSL.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_sl_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_sl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    // NOTE: the values are zeroed before anything writes them, the composed
    // step accumulating into them as the sum over primitives runs. The atom
    // pairs no pair of primitives reaches keep the zeros set here.

    std::fill(values, values + 17 * nvalues, 0.0);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(249, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 3, 8, ncols, fj, mu);

            compute_prim_sp_electron_repulsion_0(buffer, 13, 0, 6, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 16, 0, 7, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_sd_electron_repulsion_1(buffer, 34, 0, 4, 5, 13, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 37, 0, 5, 6, 16, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 40, 0, 6, 7, 19, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 43, 0, 7, 8, 22, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 46, 0, 8, 9, 25, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 49, 0, 9, 10, 28, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 52, 0, 10, 11, 31, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 55, 0, 13, 16, 40, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 61, 0, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 67, 0, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 73, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 79, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 85, 0, 34, 37, 55, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 94, 0, 37, 40, 61, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 103, 0, 40, 43, 67, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 112, 0, 43, 46, 73, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 121, 0, 46, 49, 79, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 130, 0, 55, 61, 103, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 143, 0, 61, 67, 112, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 156, 0, 67, 73, 121, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 169, 0, 85, 94, 130, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 187, 0, 94, 103, 143, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 205, 0, 103, 112, 156, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 223, 0, 130, 143, 205, ncols, alpha, beta, p);

            compute_ctr_sl_electron_repulsion_0(values, nvalues, buffer, 0, 169, 187, 223, ncols, alpha, beta, p);
        }
    }
}

}  // namespace simdt2ceri
