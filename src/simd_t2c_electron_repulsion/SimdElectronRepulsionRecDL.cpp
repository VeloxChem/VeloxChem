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


#include "SimdElectronRepulsionRecDL.hpp"

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
#include "SimdElectronRepulsionVrrRecDK.hpp"
#include "SimdElectronRepulsionVrrRecDL.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPL.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSL.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_dl_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_dl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 4911, 4539, 270, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 45, 0, 7, 8, 18, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 8, 9, 21, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 9, 10, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 10, 11, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 11, 12, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 12, 13, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 13, 14, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 14, 15, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 15, 16, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 99, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 102, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 105, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 108, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 111, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 114, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 117, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 120, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 123, 3, 9, 21, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 132, 3, 10, 24, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 141, 3, 11, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 150, 3, 12, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 159, 3, 13, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 168, 3, 14, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 177, 3, 15, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 186, 3, 16, 42, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 195, 0, 3, 21, 132, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 213, 0, 3, 24, 141, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 231, 0, 3, 27, 150, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 249, 0, 3, 30, 159, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 267, 0, 3, 33, 168, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 285, 0, 3, 36, 177, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 303, 0, 3, 39, 186, 93, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 321, 3, 9, 10, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 327, 3, 10, 11, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 333, 3, 11, 12, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 339, 3, 12, 13, 111, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 345, 3, 13, 14, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 351, 3, 14, 15, 117, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 357, 3, 15, 16, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 363, 0, 3, 99, 321, 132, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 381, 0, 3, 102, 327, 141, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 399, 0, 3, 105, 333, 150, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 417, 0, 3, 108, 339, 159, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 435, 0, 3, 111, 345, 168, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 453, 0, 3, 114, 351, 177, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 471, 0, 3, 117, 357, 186, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 489, 0, 3, 123, 363, 45, 51, 195, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 525, 0, 3, 132, 381, 51, 57, 213, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 561, 0, 3, 141, 399, 57, 63, 231, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 597, 0, 3, 150, 417, 63, 69, 249, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 633, 0, 3, 159, 435, 69, 75, 267, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 669, 0, 3, 168, 453, 75, 81, 285, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 705, 0, 3, 177, 471, 81, 87, 303, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 741, 3, 99, 102, 327, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 751, 3, 102, 105, 333, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 761, 3, 105, 108, 339, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 771, 3, 108, 111, 345, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 781, 3, 111, 114, 351, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 791, 3, 114, 117, 357, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 801, 0, 3, 321, 741, 381, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 831, 0, 3, 327, 751, 399, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 861, 0, 3, 333, 761, 417, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 891, 0, 3, 339, 771, 435, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 921, 0, 3, 345, 781, 453, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 951, 0, 3, 351, 791, 471, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 981, 0, 3, 381, 831, 195, 213, 561,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1041, 0, 3, 399, 861, 213, 231, 597,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1101, 0, 3, 417, 891, 231, 249, 633,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1161, 0, 3, 435, 921, 249, 267, 669,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1221, 0, 3, 453, 951, 267, 285, 705,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1281, 3, 321, 327, 751, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1296, 3, 327, 333, 761, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1311, 3, 333, 339, 771, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1326, 3, 339, 345, 781, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 1341, 3, 345, 351, 791, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1356, 0, 3, 741, 1281, 831, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1401, 0, 3, 751, 1296, 861, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1446, 0, 3, 761, 1311, 891, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1491, 0, 3, 771, 1326, 921, ncols, p);

            compute_prim_pg_electron_repulsion_0(buffer, 1536, 0, 3, 781, 1341, 951, ncols, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1581, 0, 3, 801, 1356, 489, 525, 981,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1671, 0, 3, 831, 1401, 525, 561, 1041,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1761, 0, 3, 861, 1446, 561, 597, 1101,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1851, 0, 3, 891, 1491, 597, 633, 1161,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 1941, 0, 3, 921, 1536, 633, 669, 1221,
                                                 ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2031, 3, 741, 751, 1296, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2052, 3, 751, 761, 1311, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2073, 3, 761, 771, 1326, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 2094, 3, 771, 781, 1341, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 2115, 0, 3, 1281, 2031, 1401, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2178, 0, 3, 1296, 2052, 1446, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2241, 0, 3, 1311, 2073, 1491, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 2304, 0, 3, 1326, 2094, 1536, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 2367, 0, 3, 1401, 2178, 981, 1041, 1761,
                                                 ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2493, 0, 3, 1446, 2241, 1041, 1101,
                                                 1851, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 2619, 0, 3, 1491, 2304, 1101, 1161,
                                                 1941, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2745, 3, 1281, 1296, 2052, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2773, 3, 1296, 1311, 2073, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 2801, 3, 1311, 1326, 2094, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 2829, 0, 3, 2031, 2745, 2178, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 2913, 0, 3, 2052, 2773, 2241, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 2997, 0, 3, 2073, 2801, 2304, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 3081, 0, 3, 2115, 2829, 1581, 1671,
                                                 2367, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 3249, 0, 3, 2178, 2913, 1671, 1761,
                                                 2493, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 3417, 0, 3, 2241, 2997, 1761, 1851,
                                                 2619, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 3585, 3, 2031, 2052, 2773, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 3621, 3, 2052, 2073, 2801, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 3657, 0, 3, 2745, 3585, 2913, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 3765, 0, 3, 2773, 3621, 2997, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 3873, 0, 3, 2913, 3765, 2367, 2493,
                                                 3417, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 4089, 3, 2745, 2773, 3621, ncols, alpha,
                                                 beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 4134, 0, 3, 3585, 4089, 3765, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 4269, 0, 3, 3657, 4134, 3081, 3249,
                                                 3873, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4539, 4269, 270, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 4809, 4539, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 4809, 17, nmax);
}

}  // namespace simdt2ceri
