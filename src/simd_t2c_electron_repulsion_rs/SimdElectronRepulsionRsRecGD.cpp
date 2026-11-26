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


#include "SimdElectronRepulsionRsRecGD.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_gd_electron_repulsion(double               *values,
                                 const size_t          nvalues,
                                 const CBasisFunction &bra,
                                 const CBasisFunction &ket,
                                 const CSimdMatrix    &coordinates,
                                 CSimdMatrix          &buffer,
                                 const double          omega) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_gd_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1591, 1336, 180, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 6, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 14, 6, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 52, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 15, 16, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 16, 17, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 17, 18, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 18, 19, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 19, 20, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 112, 0, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 122, 0, 25, 28, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 132, 0, 28, 31, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 37, 40, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 40, 43, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 43, 46, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 172, 0, 52, 58, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 187, 0, 58, 64, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 202, 0, 64, 70, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 217, 0, 82, 88, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 232, 0, 88, 94, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 247, 0, 94, 100, 162, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 262, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 265, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 268, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 271, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 274, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 277, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 280, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 283, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 286, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 295, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 304, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 313, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 322, 3, 17, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 331, 3, 18, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 340, 3, 19, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 349, 3, 20, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 358, 0, 3, 25, 295, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 376, 0, 3, 28, 304, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 394, 0, 3, 31, 313, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 412, 0, 3, 40, 331, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 430, 0, 3, 43, 340, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 448, 0, 3, 46, 349, 106, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 466, 0, 3, 64, 376, 122, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 496, 0, 3, 70, 394, 132, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 526, 0, 3, 94, 430, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 556, 0, 3, 100, 448, 162, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 586, 0, 3, 122, 496, 202, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 631, 0, 3, 152, 556, 247, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 676, 3, 9, 10, 265, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 682, 3, 10, 11, 268, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 688, 3, 11, 12, 271, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 694, 3, 17, 18, 277, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 700, 3, 18, 19, 280, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 706, 3, 19, 20, 283, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 712, 0, 3, 262, 676, 295, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 730, 0, 3, 265, 682, 304, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 748, 0, 3, 268, 688, 313, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 766, 0, 3, 274, 694, 331, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 784, 0, 3, 277, 700, 340, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 802, 0, 3, 280, 706, 349, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 820, 0, 3, 286, 712, 52, 58, 358, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 856, 0, 3, 295, 730, 58, 64, 376, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 892, 0, 3, 304, 748, 64, 70, 394, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 928, 0, 3, 322, 766, 82, 88, 412, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 964, 0, 3, 331, 784, 88, 94, 430, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1000, 0, 3, 340, 802, 94, 100, 448,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1036, 0, 3, 376, 892, 112, 122, 496,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1096, 0, 3, 430, 1000, 142, 152, 556,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1156, 0, 3, 820, 856, 466, 1036, 172,
                                                 187, 586, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1246, 0, 3, 928, 964, 526, 1096, 217,
                                                 232, 631, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1336, 1156, 180, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 1516, 1426, 15, 1, nmax);

    simdtrf::transform_g_outer(values, nvalues, buffer, 1516, 5, nmax);

    simdtrf::transform_d_inner(buffer, 1516, 1336, 15, 1, nmax);

    simdtrf::transform_g_outer(values + 45 * nvalues, nvalues, buffer, 1516, 5, nmax);
}

}  // namespace simdt2ceri
