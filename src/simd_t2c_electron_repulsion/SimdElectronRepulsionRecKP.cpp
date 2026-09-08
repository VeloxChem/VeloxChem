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


#include "SimdElectronRepulsionRecKP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_kp_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_kp_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(910, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 15, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 36, 0, 7, 8, 18, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 42, 0, 8, 9, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 48, 0, 9, 10, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 54, 0, 10, 11, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 60, 0, 11, 12, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 66, 0, 12, 13, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 72, 0, 15, 18, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 82, 0, 18, 21, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 92, 0, 21, 24, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 102, 0, 24, 27, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 27, 30, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 122, 0, 36, 42, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 137, 0, 42, 48, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 152, 0, 48, 54, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 167, 0, 54, 60, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 182, 0, 72, 82, 137, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 203, 0, 82, 92, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 224, 0, 92, 102, 167, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 245, 0, 122, 137, 203, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 273, 0, 137, 152, 224, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 301, 0, 182, 203, 273, ncols, alpha,
                                                 beta, p);

            compute_prim_pp_electron_repulsion_0(buffer, 337, 3, 13, 33, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 346, 0, 3, 30, 337, 66, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 364, 0, 3, 60, 346, 112, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 394, 0, 3, 102, 364, 167, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 439, 0, 3, 152, 394, 224, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 502, 0, 3, 203, 439, 273, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 586, 0, 3, 245, 502, 301, ncols, p);

            simdfunc::contract_primitives(buffer, 694, 586, 108, ncols);
        }
    }

    simdtrf::transform_p_inner(buffer, 802, 694, 36, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 802, 3, nmax);
}

}  // namespace simdt2ceri
