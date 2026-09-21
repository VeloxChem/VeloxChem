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


#include "SimdElectronRepulsionGeom10RecDG.hpp"

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
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryD1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_dg_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_dg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 2007, 1683, 270, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7}, ncols,
                                            fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 14, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 17, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 20, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 23, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 26, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 29, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 32, 0, 13, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 35, 0, 7, 8, 20, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 41, 0, 8, 9, 23, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 47, 0, 9, 10, 26, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 53, 0, 10, 11, 29, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 59, 0, 11, 12, 32, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 65, 0, 14, 17, 35, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 75, 0, 17, 20, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 85, 0, 20, 23, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 95, 0, 23, 26, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 26, 29, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 115, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 118, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 121, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 124, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 127, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 130, 3, 13, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 133, 3, 9, 23, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 142, 3, 10, 26, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 151, 3, 11, 29, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 160, 3, 12, 32, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 169, 0, 3, 20, 133, 41, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 187, 0, 3, 23, 142, 47, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 205, 0, 3, 26, 151, 53, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 223, 0, 3, 29, 160, 59, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 241, 0, 3, 41, 187, 85, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 271, 0, 3, 47, 205, 95, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 301, 0, 3, 53, 223, 105, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 331, 3, 7, 8, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 337, 3, 8, 9, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 343, 3, 9, 10, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 349, 3, 10, 11, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 355, 3, 11, 12, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 361, 0, 3, 121, 343, 142, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 379, 0, 3, 124, 349, 151, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 397, 0, 3, 127, 355, 160, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 415, 0, 3, 133, 361, 35, 41, 187, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 451, 0, 3, 142, 379, 41, 47, 205, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 487, 0, 3, 151, 397, 47, 53, 223, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 523, 0, 3, 169, 415, 65, 75, 241, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 583, 0, 3, 187, 451, 75, 85, 271, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 643, 0, 3, 205, 487, 85, 95, 301, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 703, 3, 115, 118, 337, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 713, 3, 118, 121, 343, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 723, 3, 121, 124, 349, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 733, 3, 124, 127, 355, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 743, 0, 3, 337, 713, 361, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 773, 0, 3, 343, 723, 379, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 803, 0, 3, 349, 733, 397, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 833, 0, 3, 361, 773, 169, 187, 451,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 893, 0, 3, 379, 803, 187, 205, 487,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 953, 0, 3, 451, 893, 241, 271, 643,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1053, 3, 331, 337, 713, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1068, 3, 343, 349, 733, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1083, 0, 3, 703, 1053, 743, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1128, 0, 3, 723, 1068, 803, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1173, 0, 3, 773, 1128, 415, 451, 893,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 1263, 0, 3, 833, 1173, 523, 583, 953,
                                                 ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 1413, 1083, 1263, 1, 15, ncols, alpha);

            simdgeo::geom_d_y(buffer, 1503, 1083, 1263, 1, 15, ncols, alpha);

            simdgeo::geom_d_z(buffer, 1593, 1083, 1263, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 1683, 1413, 270, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 1953, 1683, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 1953, 9, nmax);

    simdtrf::transform_g_inner(buffer, 1953, 1773, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 45 * nvalues, nvalues, buffer, 1953, 9, nmax);

    simdtrf::transform_g_inner(buffer, 1953, 1863, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 90 * nvalues, nvalues, buffer, 1953, 9, nmax);
}

}  // namespace simdt2ceri
