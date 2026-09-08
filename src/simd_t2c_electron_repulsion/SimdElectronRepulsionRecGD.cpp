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


#include "SimdElectronRepulsionRecGD.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdTransformGD.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_gd_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(412, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 6, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 29, 0, 7, 8, 14, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 32, 0, 8, 9, 17, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 35, 0, 9, 10, 20, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 38, 0, 10, 11, 23, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 41, 0, 11, 12, 26, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 44, 0, 14, 17, 35, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 50, 0, 17, 20, 38, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 59, 0, 20, 23, 41, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_1(buffer, 68, 0, 29, 32, 44, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_2(buffer, 74, 0, 32, 35, 50, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 80, 0, 35, 38, 59, ncols, alpha, beta, p);

            compute_prim_pp_electron_repulsion_2(buffer, 92, 3, 9, 17, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 95, 3, 10, 20, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 98, 3, 11, 23, ncols, p);

            compute_prim_dp_electron_repulsion_7(buffer, 101, 3, 17, 35, ncols, p);

            compute_prim_dp_electron_repulsion_5(buffer, 104, 0, 3, 20, 98, 38, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 116, 3, 23, 41, ncols, p);

            compute_prim_fp_electron_repulsion_5(buffer, 125, 0, 3, 35, 104, 50, ncols, p);

            compute_prim_fp_electron_repulsion_6(buffer, 141, 0, 3, 38, 116, 59, ncols, p);

            compute_prim_gp_electron_repulsion_1(buffer, 157, 0, 3, 50, 141, 80, ncols, p);

            compute_prim_dd_electron_repulsion_8(buffer, 190, 3, 92, 29, 32, 101, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_9(buffer, 193, 3, 95, 32, 35, 104, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 196, 3, 98, 35, 38, 116, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_4(buffer, 205, 0, 3, 104, 196, 44, 50, 141, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 232, 0, 3, 190, 193, 125, 205, 68, 74, 157, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 322, 232, 90, ncols);
        }
    }

    simdtrf::transform_gd(values, nvalues, buffer, 322, nmax);
}

}  // namespace simdt2ceri
