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


#include "SimdElectronRepulsionRecHI.hpp"

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
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_hi_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 14942, 14081, 588, nvalues);

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

            compute_prim_gs_electron_repulsion_0(buffer, 195, 0, 51, 57, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 210, 0, 57, 63, 135, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 225, 0, 63, 69, 145, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 240, 0, 69, 75, 155, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 255, 0, 75, 81, 165, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 270, 0, 81, 87, 175, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 285, 0, 87, 93, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_hs_electron_repulsion_0(buffer, 300, 0, 105, 115, 195, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 321, 0, 115, 125, 210, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 342, 0, 125, 135, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 363, 0, 135, 145, 240, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 384, 0, 145, 155, 255, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 405, 0, 155, 165, 270, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 426, 0, 165, 175, 285, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 447, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 450, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 453, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 456, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 459, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 462, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 465, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 468, 3, 17, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 471, 3, 9, 27, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 480, 3, 10, 30, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 489, 3, 11, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 498, 3, 12, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 507, 3, 13, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 516, 3, 14, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 525, 3, 15, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 534, 3, 16, 48, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 543, 0, 3, 24, 471, 57, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 561, 0, 3, 27, 480, 63, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 579, 0, 3, 30, 489, 69, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 597, 0, 3, 33, 498, 75, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 615, 0, 3, 36, 507, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 633, 0, 3, 39, 516, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 651, 0, 3, 42, 525, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 669, 0, 3, 45, 534, 99, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 687, 0, 3, 57, 561, 125, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 717, 0, 3, 63, 579, 135, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 747, 0, 3, 69, 597, 145, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 777, 0, 3, 75, 615, 155, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 807, 0, 3, 81, 633, 165, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 837, 0, 3, 87, 651, 175, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 867, 0, 3, 93, 669, 185, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 897, 0, 3, 125, 717, 210, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 942, 0, 3, 135, 747, 225, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 987, 0, 3, 145, 777, 240, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1032, 0, 3, 155, 807, 255, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1077, 0, 3, 165, 837, 270, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1122, 0, 3, 175, 867, 285, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1167, 0, 3, 210, 942, 342, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1230, 0, 3, 225, 987, 363, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1293, 0, 3, 240, 1032, 384, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1356, 0, 3, 255, 1077, 405, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 1419, 0, 3, 270, 1122, 426, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1482, 3, 9, 10, 450, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1488, 3, 10, 11, 453, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1494, 3, 11, 12, 456, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1500, 3, 12, 13, 459, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1506, 3, 13, 14, 462, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1512, 3, 14, 15, 465, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1518, 3, 15, 16, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1524, 0, 3, 447, 1482, 480, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1542, 0, 3, 450, 1488, 489, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1560, 0, 3, 453, 1494, 498, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1578, 0, 3, 456, 1500, 507, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1596, 0, 3, 459, 1506, 516, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1614, 0, 3, 462, 1512, 525, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1632, 0, 3, 465, 1518, 534, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1650, 0, 3, 471, 1524, 51, 57, 561,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1686, 0, 3, 480, 1542, 57, 63, 579,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1722, 0, 3, 489, 1560, 63, 69, 597,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1758, 0, 3, 498, 1578, 69, 75, 615,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1794, 0, 3, 507, 1596, 75, 81, 633,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1830, 0, 3, 516, 1614, 81, 87, 651,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1866, 0, 3, 525, 1632, 87, 93, 669,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1902, 0, 3, 543, 1650, 105, 115, 687,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1962, 0, 3, 561, 1686, 115, 125, 717,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2022, 0, 3, 579, 1722, 125, 135, 747,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2082, 0, 3, 597, 1758, 135, 145, 777,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2142, 0, 3, 615, 1794, 145, 155, 807,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2202, 0, 3, 633, 1830, 155, 165, 837,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 2262, 0, 3, 651, 1866, 165, 175, 867,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2322, 0, 3, 1650, 1686, 717, 2022, 195,
                                                 210, 942, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2412, 0, 3, 1686, 1722, 747, 2082, 210,
                                                 225, 987, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2502, 0, 3, 1722, 1758, 777, 2142, 225,
                                                 240, 1032, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2592, 0, 3, 1758, 1794, 807, 2202, 240,
                                                 255, 1077, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 2682, 0, 3, 1794, 1830, 837, 2262, 255,
                                                 270, 1122, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2772, 0, 3, 1902, 1962, 897, 2322, 300,
                                                 321, 1167, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 2898, 0, 3, 1962, 2022, 942, 2412, 321,
                                                 342, 1230, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3024, 0, 3, 2022, 2082, 987, 2502, 342,
                                                 363, 1293, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3150, 0, 3, 2082, 2142, 1032, 2592, 363,
                                                 384, 1356, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 3276, 0, 3, 2142, 2202, 1077, 2682, 384,
                                                 405, 1419, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3402, 3, 447, 450, 1488, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3412, 3, 450, 453, 1494, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3422, 3, 453, 456, 1500, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3432, 3, 456, 459, 1506, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3442, 3, 459, 462, 1512, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 3452, 3, 462, 465, 1518, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 3462, 0, 3, 1482, 3402, 1542, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3492, 0, 3, 1488, 3412, 1560, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3522, 0, 3, 1494, 3422, 1578, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3552, 0, 3, 1500, 3432, 1596, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3582, 0, 3, 1506, 3442, 1614, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 3612, 0, 3, 1512, 3452, 1632, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 3642, 0, 3, 1524, 3462, 543, 561, 1686,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3702, 0, 3, 1542, 3492, 561, 579, 1722,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3762, 0, 3, 1560, 3522, 579, 597, 1758,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3822, 0, 3, 1578, 3552, 597, 615, 1794,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3882, 0, 3, 1596, 3582, 615, 633, 1830,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 3942, 0, 3, 1614, 3612, 633, 651, 1866,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4002, 0, 3, 1686, 3702, 687, 717, 2022,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4102, 0, 3, 1722, 3762, 717, 747, 2082,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4202, 0, 3, 1758, 3822, 747, 777, 2142,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4302, 0, 3, 1794, 3882, 777, 807, 2202,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 4402, 0, 3, 1830, 3942, 807, 837, 2262,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4502, 0, 3, 3642, 3702, 2022, 4102, 897,
                                                 942, 2412, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4652, 0, 3, 3702, 3762, 2082, 4202, 942,
                                                 987, 2502, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4802, 0, 3, 3762, 3822, 2142, 4302, 987,
                                                 1032, 2592, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 4952, 0, 3, 3822, 3882, 2202, 4402,
                                                 1032, 1077, 2682, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5102, 0, 3, 4002, 4102, 2412, 4652,
                                                 1167, 1230, 3024, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5312, 0, 3, 4102, 4202, 2502, 4802,
                                                 1230, 1293, 3150, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 5522, 0, 3, 4202, 4302, 2592, 4952,
                                                 1293, 1356, 3276, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5732, 3, 1482, 1488, 3412, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5747, 3, 1488, 1494, 3422, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5762, 3, 1494, 1500, 3432, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5777, 3, 1500, 1506, 3442, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 5792, 3, 1506, 1512, 3452, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 5807, 0, 3, 3402, 5732, 3492, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5852, 0, 3, 3412, 5747, 3522, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5897, 0, 3, 3422, 5762, 3552, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5942, 0, 3, 3432, 5777, 3582, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 5987, 0, 3, 3442, 5792, 3612, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 6032, 0, 3, 3462, 5807, 1650, 1686,
                                                 3702, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6122, 0, 3, 3492, 5852, 1686, 1722,
                                                 3762, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6212, 0, 3, 3522, 5897, 1722, 1758,
                                                 3822, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6302, 0, 3, 3552, 5942, 1758, 1794,
                                                 3882, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 6392, 0, 3, 3582, 5987, 1794, 1830,
                                                 3942, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6482, 0, 3, 3642, 6032, 1902, 1962,
                                                 4002, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6632, 0, 3, 3702, 6122, 1962, 2022,
                                                 4102, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6782, 0, 3, 3762, 6212, 2022, 2082,
                                                 4202, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 6932, 0, 3, 3822, 6302, 2082, 2142,
                                                 4302, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 7082, 0, 3, 3882, 6392, 2142, 2202,
                                                 4402, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7232, 0, 3, 6032, 6122, 4102, 6782,
                                                 2322, 2412, 4652, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7457, 0, 3, 6122, 6212, 4202, 6932,
                                                 2412, 2502, 4802, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 7682, 0, 3, 6212, 6302, 4302, 7082,
                                                 2502, 2592, 4952, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 7907, 0, 3, 6482, 6632, 4502, 7232,
                                                 2772, 2898, 5102, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 8222, 0, 3, 6632, 6782, 4652, 7457,
                                                 2898, 3024, 5312, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 8537, 0, 3, 6782, 6932, 4802, 7682,
                                                 3024, 3150, 5522, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8852, 3, 3402, 3412, 5747, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8873, 3, 3412, 3422, 5762, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8894, 3, 3422, 3432, 5777, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 8915, 3, 3432, 3442, 5792, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 8936, 0, 3, 5732, 8852, 5852, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 8999, 0, 3, 5747, 8873, 5897, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 9062, 0, 3, 5762, 8894, 5942, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 9125, 0, 3, 5777, 8915, 5987, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 9188, 0, 3, 5807, 8936, 3642, 3702,
                                                 6122, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9314, 0, 3, 5852, 8999, 3702, 3762,
                                                 6212, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9440, 0, 3, 5897, 9062, 3762, 3822,
                                                 6302, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 9566, 0, 3, 5942, 9125, 3822, 3882,
                                                 6392, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9692, 0, 3, 6122, 9314, 4002, 4102,
                                                 6782, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 9902, 0, 3, 6212, 9440, 4102, 4202,
                                                 6932, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 10112, 0, 3, 6302, 9566, 4202, 4302,
                                                 7082, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 10322, 0, 3, 9188, 9314, 6782, 9902,
                                                 4502, 4652, 7457, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 10637, 0, 3, 9314, 9440, 6932, 10112,
                                                 4652, 4802, 7682, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 10952, 0, 3, 9692, 9902, 7457, 10637,
                                                 5102, 5312, 8537, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11393, 3, 5732, 5747, 8873, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11421, 3, 5747, 5762, 8894, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 11449, 3, 5762, 5777, 8915, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 11477, 0, 3, 8852, 11393, 8999, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11561, 0, 3, 8873, 11421, 9062, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 11645, 0, 3, 8894, 11449, 9125, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 11729, 0, 3, 8936, 11477, 6032, 6122,
                                                 9314, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 11897, 0, 3, 8999, 11561, 6122, 6212,
                                                 9440, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 12065, 0, 3, 9062, 11645, 6212, 6302,
                                                 9566, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 12233, 0, 3, 9188, 11729, 6482, 6632,
                                                 9692, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 12513, 0, 3, 9314, 11897, 6632, 6782,
                                                 9902, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 12793, 0, 3, 9440, 12065, 6782, 6932,
                                                 10112, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 13073, 0, 3, 11729, 11897, 9902, 12793,
                                                 7232, 7457, 10637, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 13493, 0, 3, 12233, 12513, 10322, 13073,
                                                 7907, 8222, 10952, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 14081, 13493, 588, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 14669, 14081, 21, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 14669, 13, nmax);
}

}  // namespace simdt2ceri
