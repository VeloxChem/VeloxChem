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


#include "SimdElectronRepulsionGeom10RecDI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionGeom10VrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDH.hpp"
#include "SimdElectronRepulsionVrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_di_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_di_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 5246, 4664, 504, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 43, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 49, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 55, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 61, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 12, 13, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 79, 0, 13, 14, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 85, 0, 16, 19, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 95, 0, 19, 22, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 22, 25, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 115, 0, 25, 28, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 28, 31, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 31, 34, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 34, 37, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 155, 3, 8, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 170, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 173, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 176, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 179, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 188, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 197, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 206, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 215, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 224, 3, 14, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 233, 0, 3, 22, 179, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 251, 0, 3, 25, 188, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 269, 0, 3, 28, 197, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 287, 0, 3, 31, 206, 67, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 305, 0, 3, 34, 215, 73, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 323, 0, 3, 37, 224, 79, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 341, 0, 3, 49, 251, 105, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 371, 0, 3, 55, 269, 115, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 401, 0, 3, 61, 287, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 431, 0, 3, 67, 305, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 461, 0, 3, 73, 323, 145, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 491, 3, 7, 8, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 497, 3, 8, 9, 161, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 503, 3, 9, 10, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 509, 3, 10, 11, 167, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 515, 3, 11, 12, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 521, 3, 12, 13, 173, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 527, 3, 13, 14, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 533, 0, 3, 161, 503, 188, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 551, 0, 3, 164, 509, 197, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 569, 0, 3, 167, 515, 206, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 587, 0, 3, 170, 521, 215, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 605, 0, 3, 173, 527, 224, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 623, 0, 3, 179, 533, 43, 49, 251, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 659, 0, 3, 188, 551, 49, 55, 269, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 695, 0, 3, 197, 569, 55, 61, 287, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 731, 0, 3, 206, 587, 61, 67, 305, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 767, 0, 3, 215, 605, 67, 73, 323, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 803, 0, 3, 233, 623, 85, 95, 341, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 863, 0, 3, 251, 659, 95, 105, 371,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 923, 0, 3, 269, 695, 105, 115, 401,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 983, 0, 3, 287, 731, 115, 125, 431,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1043, 0, 3, 305, 767, 125, 135, 461,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1103, 3, 155, 158, 497, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1113, 3, 158, 161, 503, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1123, 3, 161, 164, 509, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1133, 3, 164, 167, 515, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1143, 3, 167, 170, 521, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1153, 3, 170, 173, 527, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1163, 0, 3, 503, 1123, 551, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1193, 0, 3, 509, 1133, 569, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1223, 0, 3, 515, 1143, 587, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1253, 0, 3, 521, 1153, 605, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1283, 0, 3, 533, 1163, 233, 251, 659,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1343, 0, 3, 551, 1193, 251, 269, 695,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1403, 0, 3, 569, 1223, 269, 287, 731,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1463, 0, 3, 587, 1253, 287, 305, 767,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1523, 0, 3, 659, 1343, 341, 371, 923,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1623, 0, 3, 695, 1403, 371, 401, 983,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1723, 0, 3, 731, 1463, 401, 431, 1043,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1823, 3, 491, 497, 1113, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1838, 3, 497, 503, 1123, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1853, 3, 503, 509, 1133, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1868, 3, 509, 515, 1143, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1883, 3, 515, 521, 1153, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1898, 0, 3, 1123, 1853, 1193, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 1943, 0, 3, 1133, 1868, 1223, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 1988, 0, 3, 1143, 1883, 1253, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2033, 0, 3, 1163, 1898, 623, 659, 1343,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2123, 0, 3, 1193, 1943, 659, 695, 1403,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2213, 0, 3, 1223, 1988, 695, 731, 1463,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2303, 0, 3, 1283, 2033, 803, 863, 1523,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2453, 0, 3, 1343, 2123, 863, 923, 1623,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2603, 0, 3, 1403, 2213, 923, 983, 1723,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2753, 3, 1103, 1113, 1838, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2774, 3, 1113, 1123, 1853, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2795, 3, 1123, 1133, 1868, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2816, 3, 1133, 1143, 1883, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 2837, 0, 3, 1838, 2774, 1898, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2900, 0, 3, 1853, 2795, 1943, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2963, 0, 3, 1868, 2816, 1988, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 3026, 0, 3, 1898, 2900, 1283, 1343,
                                                 2123, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 3152, 0, 3, 1943, 2963, 1343, 1403,
                                                 2213, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 3278, 0, 3, 2123, 3152, 1523, 1623,
                                                 2603, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3488, 3, 1823, 1838, 2774, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3516, 3, 1853, 1868, 2816, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 3544, 0, 3, 2753, 3488, 2837, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 3628, 0, 3, 2795, 3516, 2963, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 3712, 0, 3, 2900, 3628, 2033, 2123,
                                                 3152, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 3880, 0, 3, 3026, 3712, 2303, 2453,
                                                 3278, ncols, alpha, beta, p);

            compute_prim_geom_10_di_electron_repulsion_0(buffer, 4160, 3544, 3880, ncols,
                                                         alpha);

            compute_prim_geom_10_di_electron_repulsion_1(buffer, 4328, 3544, 3880, ncols,
                                                         alpha);

            compute_prim_geom_10_di_electron_repulsion_2(buffer, 4496, 3544, 3880, ncols,
                                                         alpha);

            simdfunc::contract_primitives(buffer, 4664, 4160, 504, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 5168, 4664, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 5168, 13, nmax);

    simdtrf::transform_i_inner(buffer, 5168, 4832, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 65 * nvalues, nvalues, buffer, 5168, 13, nmax);

    simdtrf::transform_i_inner(buffer, 5168, 5000, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 130 * nvalues, nvalues, buffer, 5168, 13, nmax);
}

}  // namespace simdt2ceri
