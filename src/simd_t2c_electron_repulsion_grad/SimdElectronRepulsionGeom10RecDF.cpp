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


#include "SimdElectronRepulsionGeom10RecDF.hpp"

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
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_df_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_df_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 1120, 898, 180, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6}, ncols,
                                            fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 13, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 12, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 28, 0, 7, 8, 16, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 34, 0, 8, 9, 19, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 40, 0, 9, 10, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 46, 0, 10, 11, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 52, 0, 13, 16, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 62, 0, 16, 19, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 72, 0, 19, 22, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 82, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 85, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 88, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 91, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 94, 3, 12, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 97, 3, 9, 19, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 106, 3, 10, 22, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 115, 3, 11, 25, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 124, 0, 3, 16, 97, 34, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 142, 0, 3, 19, 106, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 160, 0, 3, 22, 115, 46, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 178, 0, 3, 28, 124, 52, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 208, 0, 3, 34, 142, 62, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 238, 0, 3, 40, 160, 72, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 268, 3, 7, 8, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 274, 3, 8, 9, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 280, 3, 9, 10, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 286, 3, 10, 11, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 292, 0, 3, 85, 274, 97, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 310, 0, 3, 88, 280, 106, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 328, 0, 3, 91, 286, 115, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 346, 0, 3, 97, 310, 28, 34, 142, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 382, 0, 3, 106, 328, 34, 40, 160, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 418, 0, 3, 142, 382, 52, 62, 238, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 478, 3, 82, 85, 274, ncols, alpha, beta,
                                                 p);

            compute_prim_sf_electron_repulsion_0(buffer, 488, 3, 88, 91, 286, ncols, alpha, beta,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 498, 0, 3, 268, 478, 292, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 528, 0, 3, 280, 488, 328, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 558, 0, 3, 310, 528, 124, 142, 382,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 618, 0, 3, 346, 558, 178, 208, 418,
                                                 ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 718, 498, 618, 1, 10, ncols, alpha);

            simdgeo::geom_d_y(buffer, 778, 498, 618, 1, 10, ncols, alpha);

            simdgeo::geom_d_z(buffer, 838, 498, 618, 1, 10, ncols, alpha);

            simdfunc::contract_primitives(buffer, 898, 718, 180, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 1078, 898, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 1078, 7, nmax);

    simdtrf::transform_f_inner(buffer, 1078, 958, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 35 * nvalues, nvalues, buffer, 1078, 7, nmax);

    simdtrf::transform_f_inner(buffer, 1078, 1018, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 70 * nvalues, nvalues, buffer, 1078, 7, nmax);
}

}  // namespace simdt2ceri
