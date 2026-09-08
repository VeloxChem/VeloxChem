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


#include "SimdElectronRepulsionRecPL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPL.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSL.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformPL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_pl_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_pl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(665, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9}, ncols, fj, mu);

            compute_prim_sp_electron_repulsion_0(buffer, 16, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 19, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 22, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 25, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 28, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 31, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 34, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 37, 3, 15, ncols);

            compute_prim_sd_electron_repulsion_1(buffer, 40, 3, 7, 8, 19, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 43, 3, 8, 9, 22, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 46, 3, 9, 10, 25, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 49, 3, 10, 11, 28, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 52, 3, 11, 12, 31, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 55, 3, 12, 13, 34, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 58, 3, 13, 14, 37, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 61, 3, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 67, 3, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 73, 3, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 79, 3, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 85, 3, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 91, 3, 31, 34, 58, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 97, 3, 40, 43, 67, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 106, 3, 43, 46, 73, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 115, 3, 46, 49, 79, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 124, 3, 49, 52, 85, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 133, 3, 52, 55, 91, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 142, 3, 61, 67, 106, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 155, 3, 67, 73, 115, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 168, 3, 73, 79, 124, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 181, 3, 79, 85, 133, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 194, 3, 97, 106, 155, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_4(buffer, 212, 3, 106, 115, 168, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 233, 3, 115, 124, 181, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_2(buffer, 251, 0, 168, 233, ncols, p);

            compute_prim_sk_electron_repulsion_2(buffer, 260, 3, 142, 155, 212, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_3(buffer, 293, 3, 155, 168, 233, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_1(buffer, 317, 0, 3, 212, 293, 251, ncols, p);

            compute_prim_sl_electron_repulsion_0(buffer, 371, 3, 194, 212, 293, ncols, alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 395, 0, 3, 260, 371, 317, ncols, p);

            simdfunc::contract_primitives(buffer, 530, 395, 135, ncols);
        }
    }

    simdtrf::transform_pl(values, nvalues, buffer, 530, nmax);
}

}  // namespace simdt2ceri
