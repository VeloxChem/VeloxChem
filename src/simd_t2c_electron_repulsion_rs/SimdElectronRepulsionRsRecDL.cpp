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


#include "SimdElectronRepulsionRsRecDL.hpp"

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
compute_rs_dl_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_dl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 9714, 9072, 540, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 10, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 18, 10, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 29, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 84, 0, 7, 8, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 90, 0, 8, 9, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 9, 10, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 10, 11, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 11, 12, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 12, 13, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 13, 14, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 14, 15, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 15, 16, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 19, 20, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 20, 21, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 21, 22, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 22, 23, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 23, 24, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 24, 25, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 25, 26, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 26, 27, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 27, 28, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_sp_electron_repulsion_0(buffer, 192, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 195, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 198, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 201, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 204, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 207, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 210, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 213, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 216, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 219, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 222, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 225, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 228, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 231, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 234, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 237, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 240, 3, 9, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 249, 3, 10, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 258, 3, 11, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 267, 3, 12, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 276, 3, 13, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 285, 3, 14, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 294, 3, 15, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 303, 3, 16, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 312, 3, 21, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 321, 3, 22, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 330, 3, 23, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 339, 3, 24, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 348, 3, 25, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 357, 3, 26, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 366, 3, 27, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 375, 3, 28, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 384, 0, 3, 33, 249, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 402, 0, 3, 36, 258, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 420, 0, 3, 39, 267, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 438, 0, 3, 42, 276, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 456, 0, 3, 45, 285, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 474, 0, 3, 48, 294, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 492, 0, 3, 51, 303, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 510, 0, 3, 60, 321, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 528, 0, 3, 63, 330, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 546, 0, 3, 66, 339, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 564, 0, 3, 69, 348, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 582, 0, 3, 72, 357, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 600, 0, 3, 75, 366, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 618, 0, 3, 78, 375, 186, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 636, 3, 9, 10, 195, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 642, 3, 10, 11, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 648, 3, 11, 12, 201, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 654, 3, 12, 13, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 660, 3, 13, 14, 207, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 666, 3, 14, 15, 210, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 672, 3, 15, 16, 213, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 678, 3, 21, 22, 219, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 684, 3, 22, 23, 222, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 690, 3, 23, 24, 225, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 696, 3, 24, 25, 228, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 702, 3, 25, 26, 231, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 708, 3, 26, 27, 234, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 714, 3, 27, 28, 237, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 720, 0, 3, 192, 636, 249, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 738, 0, 3, 195, 642, 258, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 756, 0, 3, 198, 648, 267, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 774, 0, 3, 201, 654, 276, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 792, 0, 3, 204, 660, 285, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 810, 0, 3, 207, 666, 294, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 828, 0, 3, 210, 672, 303, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 846, 0, 3, 216, 678, 321, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 864, 0, 3, 219, 684, 330, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 882, 0, 3, 222, 690, 339, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 900, 0, 3, 225, 696, 348, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 918, 0, 3, 228, 702, 357, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 936, 0, 3, 231, 708, 366, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 954, 0, 3, 234, 714, 375, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 972, 0, 3, 240, 720, 84, 90, 384, ncols,
                                                 alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1008, 0, 3, 249, 738, 90, 96, 402,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1044, 0, 3, 258, 756, 96, 102, 420,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1080, 0, 3, 267, 774, 102, 108, 438,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1116, 0, 3, 276, 792, 108, 114, 456,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1152, 0, 3, 285, 810, 114, 120, 474,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1188, 0, 3, 294, 828, 120, 126, 492,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1224, 0, 3, 312, 846, 138, 144, 510,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1260, 0, 3, 321, 864, 144, 150, 528,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1296, 0, 3, 330, 882, 150, 156, 546,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1332, 0, 3, 339, 900, 156, 162, 564,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1368, 0, 3, 348, 918, 162, 168, 582,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1404, 0, 3, 357, 936, 168, 174, 600,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1440, 0, 3, 366, 954, 174, 180, 618,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1476, 3, 192, 195, 642, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1486, 3, 195, 198, 648, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1496, 3, 198, 201, 654, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1506, 3, 201, 204, 660, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1516, 3, 204, 207, 666, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1526, 3, 207, 210, 672, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1536, 3, 216, 219, 684, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1546, 3, 219, 222, 690, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1556, 3, 222, 225, 696, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1566, 3, 225, 228, 702, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1576, 3, 228, 231, 708, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1586, 3, 231, 234, 714, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1596, 0, 3, 636, 1476, 738, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1626, 0, 3, 642, 1486, 756, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1656, 0, 3, 648, 1496, 774, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1686, 0, 3, 654, 1506, 792, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1716, 0, 3, 660, 1516, 810, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1746, 0, 3, 666, 1526, 828, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1776, 0, 3, 678, 1536, 864, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1806, 0, 3, 684, 1546, 882, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1836, 0, 3, 690, 1556, 900, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1866, 0, 3, 696, 1566, 918, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1896, 0, 3, 702, 1576, 936, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1926, 0, 3, 708, 1586, 954, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 1956, 0, 3, 738, 1626, 384, 402, 1044,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2016, 0, 3, 756, 1656, 402, 420, 1080,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2076, 0, 3, 774, 1686, 420, 438, 1116,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2136, 0, 3, 792, 1716, 438, 456, 1152,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2196, 0, 3, 810, 1746, 456, 474, 1188,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2256, 0, 3, 864, 1806, 510, 528, 1296,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2316, 0, 3, 882, 1836, 528, 546, 1332,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2376, 0, 3, 900, 1866, 546, 564, 1368,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2436, 0, 3, 918, 1896, 564, 582, 1404,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2496, 0, 3, 936, 1926, 582, 600, 1440,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2556, 3, 636, 642, 1486, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2571, 3, 642, 648, 1496, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2586, 3, 648, 654, 1506, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2601, 3, 654, 660, 1516, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2616, 3, 660, 666, 1526, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2631, 3, 678, 684, 1546, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2646, 3, 684, 690, 1556, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2661, 3, 690, 696, 1566, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2676, 3, 696, 702, 1576, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2691, 3, 702, 708, 1586, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2706, 0, 3, 1476, 2556, 1626, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2751, 0, 3, 1486, 2571, 1656, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2796, 0, 3, 1496, 2586, 1686, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2841, 0, 3, 1506, 2601, 1716, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2886, 0, 3, 1516, 2616, 1746, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2931, 0, 3, 1536, 2631, 1806, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2976, 0, 3, 1546, 2646, 1836, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3021, 0, 3, 1556, 2661, 1866, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3066, 0, 3, 1566, 2676, 1896, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3111, 0, 3, 1576, 2691, 1926, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 3156, 0, 3, 1596, 2706, 972, 1008, 1956,
                                                 ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3246, 0, 3, 1626, 2751, 1008, 1044,
                                                 2016, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3336, 0, 3, 1656, 2796, 1044, 1080,
                                                 2076, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3426, 0, 3, 1686, 2841, 1080, 1116,
                                                 2136, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3516, 0, 3, 1716, 2886, 1116, 1152,
                                                 2196, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3606, 0, 3, 1776, 2931, 1224, 1260,
                                                 2256, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3696, 0, 3, 1806, 2976, 1260, 1296,
                                                 2316, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3786, 0, 3, 1836, 3021, 1296, 1332,
                                                 2376, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3876, 0, 3, 1866, 3066, 1332, 1368,
                                                 2436, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3966, 0, 3, 1896, 3111, 1368, 1404,
                                                 2496, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4056, 3, 1476, 1486, 2571, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4077, 3, 1486, 1496, 2586, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4098, 3, 1496, 1506, 2601, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4119, 3, 1506, 1516, 2616, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4140, 3, 1536, 1546, 2646, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4161, 3, 1546, 1556, 2661, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4182, 3, 1556, 1566, 2676, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4203, 3, 1566, 1576, 2691, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 4224, 0, 3, 2556, 4056, 2751, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4287, 0, 3, 2571, 4077, 2796, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4350, 0, 3, 2586, 4098, 2841, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4413, 0, 3, 2601, 4119, 2886, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4476, 0, 3, 2631, 4140, 2976, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4539, 0, 3, 2646, 4161, 3021, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4602, 0, 3, 2661, 4182, 3066, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4665, 0, 3, 2676, 4203, 3111, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 4728, 0, 3, 2751, 4287, 1956, 2016,
                                                 3336, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4854, 0, 3, 2796, 4350, 2016, 2076,
                                                 3426, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 4980, 0, 3, 2841, 4413, 2076, 2136,
                                                 3516, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5106, 0, 3, 2976, 4539, 2256, 2316,
                                                 3786, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5232, 0, 3, 3021, 4602, 2316, 2376,
                                                 3876, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5358, 0, 3, 3066, 4665, 2376, 2436,
                                                 3966, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5484, 3, 2556, 2571, 4077, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5512, 3, 2571, 2586, 4098, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5540, 3, 2586, 2601, 4119, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5568, 3, 2631, 2646, 4161, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5596, 3, 2646, 2661, 4182, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5624, 3, 2661, 2676, 4203, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 5652, 0, 3, 4056, 5484, 4287, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5736, 0, 3, 4077, 5512, 4350, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5820, 0, 3, 4098, 5540, 4413, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5904, 0, 3, 4140, 5568, 4539, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 5988, 0, 3, 4161, 5596, 4602, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 6072, 0, 3, 4182, 5624, 4665, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 6156, 0, 3, 4224, 5652, 3156, 3246,
                                                 4728, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6324, 0, 3, 4287, 5736, 3246, 3336,
                                                 4854, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6492, 0, 3, 4350, 5820, 3336, 3426,
                                                 4980, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6660, 0, 3, 4476, 5904, 3606, 3696,
                                                 5106, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6828, 0, 3, 4539, 5988, 3696, 3786,
                                                 5232, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6996, 0, 3, 4602, 6072, 3786, 3876,
                                                 5358, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7164, 3, 4056, 4077, 5512, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7200, 3, 4077, 4098, 5540, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7236, 3, 4140, 4161, 5596, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7272, 3, 4161, 4182, 5624, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 7308, 0, 3, 5484, 7164, 5736, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7416, 0, 3, 5512, 7200, 5820, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7524, 0, 3, 5568, 7236, 5988, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7632, 0, 3, 5596, 7272, 6072, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 7740, 0, 3, 5736, 7416, 4728, 4854,
                                                 6492, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 7956, 0, 3, 5988, 7632, 5106, 5232,
                                                 6996, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8172, 3, 5484, 5512, 7200, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8217, 3, 5568, 5596, 7272, ncols, alpha,
                                                 beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 8262, 0, 3, 7164, 8172, 7416, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 8397, 0, 3, 7236, 8217, 7632, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 8532, 0, 3, 7308, 8262, 6156, 6324,
                                                 7740, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 8802, 0, 3, 7524, 8397, 6660, 6828,
                                                 7956, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 9072, 8532, 540, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 9612, 9342, 6, 1, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 9612, 17, nmax);

    simdtrf::transform_l_inner(buffer, 9612, 9072, 6, 1, nmax);

    simdtrf::transform_d_outer(values + 85 * nvalues, nvalues, buffer, 9612, 17, nmax);
}

}  // namespace simdt2ceri
