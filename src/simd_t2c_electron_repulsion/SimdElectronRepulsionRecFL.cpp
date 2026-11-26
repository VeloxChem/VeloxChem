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


#include "SimdElectronRepulsionRecFL.hpp"

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
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFL.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fl_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 9835, 9215, 450, nvalues);

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 18, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 21, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 17, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 51, 0, 7, 8, 24, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 57, 0, 8, 9, 27, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 63, 0, 9, 10, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 69, 0, 10, 11, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 11, 12, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 12, 13, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 13, 14, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 14, 15, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 99, 0, 15, 16, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 105, 0, 18, 21, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 115, 0, 21, 24, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 125, 0, 24, 27, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 135, 0, 27, 30, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 145, 0, 30, 33, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 155, 0, 33, 36, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 36, 39, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 39, 42, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 42, 45, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 195, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 198, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 201, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 204, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 207, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 210, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 213, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 216, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 219, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 228, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 237, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 246, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 255, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 264, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 273, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 282, 3, 16, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 291, 0, 3, 24, 219, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 309, 0, 3, 27, 228, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 327, 0, 3, 30, 237, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 345, 0, 3, 33, 246, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 363, 0, 3, 36, 255, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 381, 0, 3, 39, 264, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 399, 0, 3, 42, 273, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 417, 0, 3, 45, 282, 99, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 435, 0, 3, 57, 309, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 465, 0, 3, 63, 327, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 495, 0, 3, 69, 345, 145, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 525, 0, 3, 75, 363, 155, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 555, 0, 3, 81, 381, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 585, 0, 3, 87, 399, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 615, 0, 3, 93, 417, 185, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 645, 3, 9, 10, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 651, 3, 10, 11, 201, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 657, 3, 11, 12, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 663, 3, 12, 13, 207, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 669, 3, 13, 14, 210, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 675, 3, 14, 15, 213, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 681, 3, 15, 16, 216, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 687, 0, 3, 195, 645, 228, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 705, 0, 3, 198, 651, 237, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 723, 0, 3, 201, 657, 246, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 741, 0, 3, 204, 663, 255, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 759, 0, 3, 207, 669, 264, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 777, 0, 3, 210, 675, 273, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 795, 0, 3, 213, 681, 282, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 813, 0, 3, 219, 687, 51, 57, 309, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 849, 0, 3, 228, 705, 57, 63, 327, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 885, 0, 3, 237, 723, 63, 69, 345, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 921, 0, 3, 246, 741, 69, 75, 363, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 957, 0, 3, 255, 759, 75, 81, 381, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 993, 0, 3, 264, 777, 81, 87, 399, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1029, 0, 3, 273, 795, 87, 93, 417,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1065, 0, 3, 291, 813, 105, 115, 435,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1125, 0, 3, 309, 849, 115, 125, 465,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1185, 0, 3, 327, 885, 125, 135, 495,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1245, 0, 3, 345, 921, 135, 145, 525,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1305, 0, 3, 363, 957, 145, 155, 555,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1365, 0, 3, 381, 993, 155, 165, 585,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1425, 0, 3, 399, 1029, 165, 175, 615,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1485, 3, 195, 198, 651, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1495, 3, 198, 201, 657, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1505, 3, 201, 204, 663, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1515, 3, 204, 207, 669, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1525, 3, 207, 210, 675, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1535, 3, 210, 213, 681, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1545, 0, 3, 645, 1485, 705, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1575, 0, 3, 651, 1495, 723, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1605, 0, 3, 657, 1505, 741, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1635, 0, 3, 663, 1515, 759, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1665, 0, 3, 669, 1525, 777, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1695, 0, 3, 675, 1535, 795, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1725, 0, 3, 687, 1545, 291, 309, 849,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1785, 0, 3, 705, 1575, 309, 327, 885,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1845, 0, 3, 723, 1605, 327, 345, 921,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1905, 0, 3, 741, 1635, 345, 363, 957,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 1965, 0, 3, 759, 1665, 363, 381, 993,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2025, 0, 3, 777, 1695, 381, 399, 1029,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2085, 0, 3, 849, 1785, 435, 465, 1185,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2185, 0, 3, 885, 1845, 465, 495, 1245,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2285, 0, 3, 921, 1905, 495, 525, 1305,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2385, 0, 3, 957, 1965, 525, 555, 1365,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2485, 0, 3, 993, 2025, 555, 585, 1425,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2585, 3, 645, 651, 1495, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2600, 3, 651, 657, 1505, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2615, 3, 657, 663, 1515, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2630, 3, 663, 669, 1525, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2645, 3, 669, 675, 1535, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2660, 0, 3, 1485, 2585, 1575, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2705, 0, 3, 1495, 2600, 1605, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2750, 0, 3, 1505, 2615, 1635, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2795, 0, 3, 1515, 2630, 1665, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2840, 0, 3, 1525, 2645, 1695, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 2885, 0, 3, 1545, 2660, 813, 849, 1785,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 2975, 0, 3, 1575, 2705, 849, 885, 1845,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3065, 0, 3, 1605, 2750, 885, 921, 1905,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3155, 0, 3, 1635, 2795, 921, 957, 1965,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3245, 0, 3, 1665, 2840, 957, 993, 2025,
                                                 ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3335, 0, 3, 1725, 2885, 1065, 1125,
                                                 2085, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3485, 0, 3, 1785, 2975, 1125, 1185,
                                                 2185, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3635, 0, 3, 1845, 3065, 1185, 1245,
                                                 2285, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3785, 0, 3, 1905, 3155, 1245, 1305,
                                                 2385, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 3935, 0, 3, 1965, 3245, 1305, 1365,
                                                 2485, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4085, 3, 1485, 1495, 2600, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4106, 3, 1495, 1505, 2615, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4127, 3, 1505, 1515, 2630, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4148, 3, 1515, 1525, 2645, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 4169, 0, 3, 2585, 4085, 2705, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4232, 0, 3, 2600, 4106, 2750, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4295, 0, 3, 2615, 4127, 2795, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4358, 0, 3, 2630, 4148, 2840, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 4421, 0, 3, 2660, 4169, 1725, 1785,
                                                 2975, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4547, 0, 3, 2705, 4232, 1785, 1845,
                                                 3065, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4673, 0, 3, 2750, 4295, 1845, 1905,
                                                 3155, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4799, 0, 3, 2795, 4358, 1905, 1965,
                                                 3245, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 4925, 0, 3, 2975, 4547, 2085, 2185,
                                                 3635, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 5135, 0, 3, 3065, 4673, 2185, 2285,
                                                 3785, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 5345, 0, 3, 3155, 4799, 2285, 2385,
                                                 3935, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5555, 3, 2585, 2600, 4106, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5583, 3, 2600, 2615, 4127, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5611, 3, 2615, 2630, 4148, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 5639, 0, 3, 4085, 5555, 4232, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5723, 0, 3, 4106, 5583, 4295, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5807, 0, 3, 4127, 5611, 4358, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 5891, 0, 3, 4169, 5639, 2885, 2975,
                                                 4547, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6059, 0, 3, 4232, 5723, 2975, 3065,
                                                 4673, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6227, 0, 3, 4295, 5807, 3065, 3155,
                                                 4799, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 6395, 0, 3, 4421, 5891, 3335, 3485,
                                                 4925, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 6675, 0, 3, 4547, 6059, 3485, 3635,
                                                 5135, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 6955, 0, 3, 4673, 6227, 3635, 3785,
                                                 5345, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7235, 3, 4085, 4106, 5583, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7271, 3, 4106, 4127, 5611, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 7307, 0, 3, 5555, 7235, 5723, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7415, 0, 3, 5583, 7271, 5807, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 7523, 0, 3, 5639, 7307, 4421, 4547,
                                                 6059, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 7739, 0, 3, 5723, 7415, 4547, 4673,
                                                 6227, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 7955, 0, 3, 6059, 7739, 4925, 5135,
                                                 6955, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8315, 3, 5555, 5583, 7271, ncols, alpha,
                                                 beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 8360, 0, 3, 7235, 8315, 7415, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 8495, 0, 3, 7307, 8360, 5891, 6059,
                                                 7739, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 8765, 0, 3, 7523, 8495, 6395, 6675,
                                                 7955, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 9215, 8765, 450, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 9665, 9215, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 9665, 17, nmax);
}

}  // namespace simdt2ceri
