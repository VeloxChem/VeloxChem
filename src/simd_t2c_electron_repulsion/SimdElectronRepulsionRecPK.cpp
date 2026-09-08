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


#include "SimdElectronRepulsionRecPK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformPK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_pk_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_pk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(493, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8}, ncols, fj, mu);

            compute_prim_sp_electron_repulsion_0(buffer, 15, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 18, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 21, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 24, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 27, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 30, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 33, 3, 14, ncols);

            compute_prim_sd_electron_repulsion_1(buffer, 36, 3, 7, 8, 18, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 39, 3, 8, 9, 21, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 42, 3, 9, 10, 24, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 45, 3, 10, 11, 27, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 48, 3, 11, 12, 30, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 51, 3, 12, 13, 33, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 54, 3, 15, 18, 39, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 60, 3, 18, 21, 42, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 66, 3, 21, 24, 45, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 72, 3, 24, 27, 48, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 78, 3, 27, 30, 51, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 84, 3, 36, 39, 60, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 93, 3, 39, 42, 66, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 102, 3, 42, 45, 72, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 111, 3, 45, 48, 78, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 120, 3, 54, 60, 93, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_4(buffer, 133, 3, 60, 66, 102, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_1(buffer, 149, 3, 66, 72, 111, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_2(buffer, 162, 0, 102, 149, ncols, p);

            compute_prim_si_electron_repulsion_3(buffer, 171, 3, 84, 93, 133, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 196, 3, 93, 102, 149, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_1(buffer, 214, 0, 3, 133, 196, 162, ncols, p);

            compute_prim_sk_electron_repulsion_1(buffer, 259, 3, 120, 133, 196, ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 277, 0, 3, 171, 259, 214, ncols, p);

            simdfunc::contract_primitives(buffer, 385, 277, 108, ncols);
        }
    }

    simdtrf::transform_pk(values, nvalues, buffer, 385, nmax);
}

}  // namespace simdt2ceri
