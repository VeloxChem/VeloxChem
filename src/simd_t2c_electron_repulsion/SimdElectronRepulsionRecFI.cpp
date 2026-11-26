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


#include "SimdElectronRepulsionRecFI.hpp"

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
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fi_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 4285, 3875, 280, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 155, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 158, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 161, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 164, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 167, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 170, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 173, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 182, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 191, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 200, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 209, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 218, 3, 14, 40, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 227, 0, 3, 22, 173, 49, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 245, 0, 3, 25, 182, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 263, 0, 3, 28, 191, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 281, 0, 3, 31, 200, 67, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 299, 0, 3, 34, 209, 73, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 317, 0, 3, 37, 218, 79, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 335, 0, 3, 49, 245, 105, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 365, 0, 3, 55, 263, 115, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 395, 0, 3, 61, 281, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 425, 0, 3, 67, 299, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 455, 0, 3, 73, 317, 145, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 485, 3, 9, 10, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 491, 3, 10, 11, 161, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 497, 3, 11, 12, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 503, 3, 12, 13, 167, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 509, 3, 13, 14, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 515, 0, 3, 155, 485, 182, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 533, 0, 3, 158, 491, 191, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 551, 0, 3, 161, 497, 200, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 569, 0, 3, 164, 503, 209, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 587, 0, 3, 167, 509, 218, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 605, 0, 3, 173, 515, 43, 49, 245, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 641, 0, 3, 182, 533, 49, 55, 263, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 677, 0, 3, 191, 551, 55, 61, 281, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 713, 0, 3, 200, 569, 61, 67, 299, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 749, 0, 3, 209, 587, 67, 73, 317, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 785, 0, 3, 227, 605, 85, 95, 335, ncols,
                                                 alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 845, 0, 3, 245, 641, 95, 105, 365,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 905, 0, 3, 263, 677, 105, 115, 395,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 965, 0, 3, 281, 713, 115, 125, 425,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1025, 0, 3, 299, 749, 125, 135, 455,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1085, 3, 155, 158, 491, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1095, 3, 158, 161, 497, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1105, 3, 161, 164, 503, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1115, 3, 164, 167, 509, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1125, 0, 3, 485, 1085, 533, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1155, 0, 3, 491, 1095, 551, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1185, 0, 3, 497, 1105, 569, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1215, 0, 3, 503, 1115, 587, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1245, 0, 3, 515, 1125, 227, 245, 641,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1305, 0, 3, 533, 1155, 245, 263, 677,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1365, 0, 3, 551, 1185, 263, 281, 713,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1425, 0, 3, 569, 1215, 281, 299, 749,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1485, 0, 3, 641, 1305, 335, 365, 905,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1585, 0, 3, 677, 1365, 365, 395, 965,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 1685, 0, 3, 713, 1425, 395, 425, 1025,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1785, 3, 485, 491, 1095, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1800, 3, 491, 497, 1105, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1815, 3, 497, 503, 1115, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1830, 0, 3, 1085, 1785, 1155, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 1875, 0, 3, 1095, 1800, 1185, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 1920, 0, 3, 1105, 1815, 1215, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 1965, 0, 3, 1125, 1830, 605, 641, 1305,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2055, 0, 3, 1155, 1875, 641, 677, 1365,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2145, 0, 3, 1185, 1920, 677, 713, 1425,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2235, 0, 3, 1245, 1965, 785, 845, 1485,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2385, 0, 3, 1305, 2055, 845, 905, 1585,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 2535, 0, 3, 1365, 2145, 905, 965, 1685,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2685, 3, 1085, 1095, 1800, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2706, 3, 1095, 1105, 1815, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 2727, 0, 3, 1785, 2685, 1875, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2790, 0, 3, 1800, 2706, 1920, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 2853, 0, 3, 1830, 2727, 1245, 1305,
                                                 2055, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2979, 0, 3, 1875, 2790, 1305, 1365,
                                                 2145, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 3105, 0, 3, 2055, 2979, 1485, 1585,
                                                 2535, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 3315, 3, 1785, 1800, 2706, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 3343, 0, 3, 2685, 3315, 2790, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 3427, 0, 3, 2727, 3343, 1965, 2055,
                                                 2979, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 3595, 0, 3, 2853, 3427, 2235, 2385,
                                                 3105, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 3875, 3595, 280, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 4155, 3875, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 4155, 13, nmax);
}

}  // namespace simdt2ceri
