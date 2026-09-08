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


#include "SimdElectronRepulsionRecLP.hpp"

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
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdTransformLP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_lp_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_lp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(662, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 40, 0, 7, 8, 19, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 43, 0, 8, 9, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 46, 0, 9, 10, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 49, 0, 10, 11, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 11, 12, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 12, 13, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 58, 0, 13, 14, 37, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 61, 0, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 67, 0, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 73, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 79, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 85, 0, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_4(buffer, 91, 0, 31, 34, 58, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 97, 0, 40, 43, 67, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 106, 0, 43, 46, 73, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 115, 0, 46, 49, 79, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 124, 0, 49, 52, 85, ncols, alpha, beta, p);

            compute_prim_gs_electron_repulsion_8(buffer, 133, 0, 52, 55, 91, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 142, 0, 61, 67, 106, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 155, 0, 67, 73, 115, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 168, 0, 73, 79, 124, ncols, alpha, beta, p);

            compute_prim_hs_electron_repulsion_15(buffer, 181, 0, 79, 85, 133, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_16(buffer, 194, 0, 97, 106, 155, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_8(buffer, 212, 0, 106, 115, 168, ncols, alpha, beta, p);

            compute_prim_is_electron_repulsion_16(buffer, 233, 0, 115, 124, 181, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_7(buffer, 251, 0, 142, 155, 212, ncols, alpha, beta, p);

            compute_prim_ks_electron_repulsion_8(buffer, 284, 0, 155, 168, 233, ncols, alpha, beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 311, 0, 194, 212, 284, ncols, alpha, beta, p);

            compute_prim_ip_electron_repulsion_11(buffer, 353, 3, 168, 233, ncols, p);

            compute_prim_kp_electron_repulsion_5(buffer, 356, 0, 3, 212, 353, 284, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 392, 0, 3, 251, 356, 311, ncols, p);

            simdfunc::contract_primitives(buffer, 527, 392, 135, ncols);
        }
    }

    simdtrf::transform_lp(values, nvalues, buffer, 527, nmax);
}

}  // namespace simdt2ceri
