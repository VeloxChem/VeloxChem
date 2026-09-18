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


#include "SimdElectronRepulsionGeom10RsRecDG.hpp"

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
compute_rs_geom_10_dg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_dg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 3954, 3360, 540, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 14, {1, 2, 3, 4, 5, 6, 7},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 7, 8, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 8, 9, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 9, 10, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 10, 11, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 11, 12, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 15, 16, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 16, 17, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 17, 18, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 18, 19, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 19, 20, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 124, 0, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 134, 0, 25, 28, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 144, 0, 28, 31, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 154, 0, 31, 34, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 164, 0, 34, 37, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 174, 0, 43, 46, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 184, 0, 46, 49, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 194, 0, 49, 52, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 52, 55, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 55, 58, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 224, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 227, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 230, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 233, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 236, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 239, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 242, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 245, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 248, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 251, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 254, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 257, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 260, 3, 9, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 269, 3, 10, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 278, 3, 11, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 287, 3, 12, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 296, 3, 17, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 305, 3, 18, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 314, 3, 19, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 323, 3, 20, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 332, 0, 3, 28, 260, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 350, 0, 3, 31, 269, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 368, 0, 3, 34, 278, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 386, 0, 3, 37, 287, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 404, 0, 3, 49, 296, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 422, 0, 3, 52, 305, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 440, 0, 3, 55, 314, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 458, 0, 3, 58, 323, 118, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 476, 0, 3, 70, 350, 144, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 506, 0, 3, 76, 368, 154, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 536, 0, 3, 82, 386, 164, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 566, 0, 3, 100, 422, 194, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 596, 0, 3, 106, 440, 204, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 626, 0, 3, 112, 458, 214, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 656, 3, 7, 8, 227, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 662, 3, 8, 9, 230, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 668, 3, 9, 10, 233, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 674, 3, 10, 11, 236, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 680, 3, 11, 12, 239, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 686, 3, 15, 16, 245, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 692, 3, 16, 17, 248, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 698, 3, 17, 18, 251, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 704, 3, 18, 19, 254, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 710, 3, 19, 20, 257, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 716, 0, 3, 230, 668, 269, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 734, 0, 3, 233, 674, 278, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 752, 0, 3, 236, 680, 287, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 770, 0, 3, 248, 698, 305, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 788, 0, 3, 251, 704, 314, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 806, 0, 3, 254, 710, 323, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 824, 0, 3, 260, 716, 64, 70, 350, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 860, 0, 3, 269, 734, 70, 76, 368, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 896, 0, 3, 278, 752, 76, 82, 386, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 932, 0, 3, 296, 770, 94, 100, 422,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 968, 0, 3, 305, 788, 100, 106, 440,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1004, 0, 3, 314, 806, 106, 112, 458,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1040, 0, 3, 332, 824, 124, 134, 476,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1100, 0, 3, 350, 860, 134, 144, 506,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1160, 0, 3, 368, 896, 144, 154, 536,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1220, 0, 3, 404, 932, 174, 184, 566,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1280, 0, 3, 422, 968, 184, 194, 596,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1340, 0, 3, 440, 1004, 194, 204, 626,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1400, 3, 224, 227, 662, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1410, 3, 227, 230, 668, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1420, 3, 230, 233, 674, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1430, 3, 233, 236, 680, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1440, 3, 242, 245, 692, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1450, 3, 245, 248, 698, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1460, 3, 248, 251, 704, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1470, 3, 251, 254, 710, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1480, 0, 3, 662, 1410, 716, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1510, 0, 3, 668, 1420, 734, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1540, 0, 3, 674, 1430, 752, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1570, 0, 3, 692, 1450, 770, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1600, 0, 3, 698, 1460, 788, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1630, 0, 3, 704, 1470, 806, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1660, 0, 3, 716, 1510, 332, 350, 860,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1720, 0, 3, 734, 1540, 350, 368, 896,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1780, 0, 3, 770, 1600, 404, 422, 968,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1840, 0, 3, 788, 1630, 422, 440, 1004,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1900, 0, 3, 860, 1720, 476, 506, 1160,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2000, 0, 3, 968, 1840, 566, 596, 1340,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2100, 3, 656, 662, 1410, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2115, 3, 668, 674, 1430, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2130, 3, 686, 692, 1450, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2145, 3, 698, 704, 1470, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2160, 0, 3, 1400, 2100, 1480, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2205, 0, 3, 1420, 2115, 1540, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2250, 0, 3, 1440, 2130, 1570, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2295, 0, 3, 1460, 2145, 1630, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2340, 0, 3, 1510, 2205, 824, 860, 1720,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2430, 0, 3, 1600, 2295, 932, 968, 1840,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2520, 0, 3, 1660, 2340, 1040, 1100,
                                                 1900, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2670, 0, 3, 1780, 2430, 1220, 1280,
                                                 2000, ncols, alpha, beta, p);

            simdgeo::geom_d_x(buffer, 2820, 2250, 2670, 1, 15, ncols, alpha);

            simdgeo::geom_d_y(buffer, 2910, 2250, 2670, 1, 15, ncols, alpha);

            simdgeo::geom_d_z(buffer, 3000, 2250, 2670, 1, 15, ncols, alpha);

            simdgeo::geom_d_x(buffer, 3090, 2160, 2520, 1, 15, ncols, alpha);

            simdgeo::geom_d_y(buffer, 3180, 2160, 2520, 1, 15, ncols, alpha);

            simdgeo::geom_d_z(buffer, 3270, 2160, 2520, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 3360, 3090, 270, ncols);

            simdfunc::contract_primitives(buffer, 3630, 2820, 270, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 3900, 3630, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 3900, 9, nmax);

    simdtrf::transform_g_inner(buffer, 3900, 3720, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 45 * nvalues, nvalues, buffer, 3900, 9, nmax);

    simdtrf::transform_g_inner(buffer, 3900, 3810, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 90 * nvalues, nvalues, buffer, 3900, 9, nmax);

    simdtrf::transform_g_inner(buffer, 3900, 3360, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 135 * nvalues, nvalues, buffer, 3900, 9, nmax);

    simdtrf::transform_g_inner(buffer, 3900, 3450, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 180 * nvalues, nvalues, buffer, 3900, 9, nmax);

    simdtrf::transform_g_inner(buffer, 3900, 3540, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 225 * nvalues, nvalues, buffer, 3900, 9, nmax);
}

}  // namespace simdt2ceri
