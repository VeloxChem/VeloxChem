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


#include "SimdElectronRepulsionGeom10RsRecPL.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdTransformL.hpp"
#include "SimdTransformP.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_pl_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_pl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 11169, 10308, 810, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 192, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 195, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 198, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 201, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 204, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 207, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 210, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 213, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 216, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 219, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 222, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 225, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 228, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 231, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 234, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 237, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 240, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 243, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 246, 3, 9, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 255, 3, 10, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 264, 3, 11, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 273, 3, 12, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 282, 3, 13, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 291, 3, 14, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 300, 3, 15, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 309, 3, 16, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 318, 3, 21, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 327, 3, 22, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 336, 3, 23, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 345, 3, 24, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 354, 3, 25, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 363, 3, 26, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 372, 3, 27, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 381, 3, 28, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 390, 0, 3, 33, 255, 96, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 408, 0, 3, 36, 264, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 426, 0, 3, 39, 273, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 444, 0, 3, 42, 282, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 462, 0, 3, 45, 291, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 480, 0, 3, 48, 300, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 498, 0, 3, 51, 309, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 516, 0, 3, 60, 327, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 534, 0, 3, 63, 336, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 552, 0, 3, 66, 345, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 570, 0, 3, 69, 354, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 588, 0, 3, 72, 363, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 606, 0, 3, 75, 372, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 624, 0, 3, 78, 381, 186, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 642, 3, 7, 8, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 648, 3, 8, 9, 195, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 654, 3, 9, 10, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 660, 3, 10, 11, 201, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 666, 3, 11, 12, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 672, 3, 12, 13, 207, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 678, 3, 13, 14, 210, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 684, 3, 14, 15, 213, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 690, 3, 15, 16, 216, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 696, 3, 19, 20, 219, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 702, 3, 20, 21, 222, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 708, 3, 21, 22, 225, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 714, 3, 22, 23, 228, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 720, 3, 23, 24, 231, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 726, 3, 24, 25, 234, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 732, 3, 25, 26, 237, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 738, 3, 26, 27, 240, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 744, 3, 27, 28, 243, ncols, alpha, beta,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 750, 0, 3, 195, 654, 255, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 768, 0, 3, 198, 660, 264, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 786, 0, 3, 201, 666, 273, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 804, 0, 3, 204, 672, 282, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 822, 0, 3, 207, 678, 291, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 840, 0, 3, 210, 684, 300, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 858, 0, 3, 213, 690, 309, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 876, 0, 3, 222, 708, 327, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 894, 0, 3, 225, 714, 336, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 912, 0, 3, 228, 720, 345, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 930, 0, 3, 231, 726, 354, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 948, 0, 3, 234, 732, 363, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 966, 0, 3, 237, 738, 372, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 984, 0, 3, 240, 744, 381, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1002, 0, 3, 246, 750, 84, 90, 390,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1038, 0, 3, 255, 768, 90, 96, 408,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1074, 0, 3, 264, 786, 96, 102, 426,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1110, 0, 3, 273, 804, 102, 108, 444,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1146, 0, 3, 282, 822, 108, 114, 462,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1182, 0, 3, 291, 840, 114, 120, 480,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1218, 0, 3, 300, 858, 120, 126, 498,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1254, 0, 3, 318, 876, 138, 144, 516,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1290, 0, 3, 327, 894, 144, 150, 534,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1326, 0, 3, 336, 912, 150, 156, 552,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1362, 0, 3, 345, 930, 156, 162, 570,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1398, 0, 3, 354, 948, 162, 168, 588,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1434, 0, 3, 363, 966, 168, 174, 606,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1470, 0, 3, 372, 984, 174, 180, 624,
                                                 ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1506, 3, 192, 195, 654, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1516, 3, 195, 198, 660, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1526, 3, 198, 201, 666, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1536, 3, 201, 204, 672, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1546, 3, 204, 207, 678, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1556, 3, 207, 210, 684, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1566, 3, 210, 213, 690, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1576, 3, 219, 222, 708, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1586, 3, 222, 225, 714, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1596, 3, 225, 228, 720, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1606, 3, 228, 231, 726, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1616, 3, 231, 234, 732, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1626, 3, 234, 237, 738, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1636, 3, 237, 240, 744, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1646, 0, 3, 654, 1516, 768, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1676, 0, 3, 660, 1526, 786, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1706, 0, 3, 666, 1536, 804, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1736, 0, 3, 672, 1546, 822, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1766, 0, 3, 678, 1556, 840, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1796, 0, 3, 684, 1566, 858, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1826, 0, 3, 708, 1586, 894, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1856, 0, 3, 714, 1596, 912, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1886, 0, 3, 720, 1606, 930, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1916, 0, 3, 726, 1616, 948, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1946, 0, 3, 732, 1626, 966, ncols, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1976, 0, 3, 738, 1636, 984, ncols, p);

            compute_prim_df_electron_repulsion_0(buffer, 2006, 0, 3, 768, 1676, 390, 408, 1074,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2066, 0, 3, 786, 1706, 408, 426, 1110,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2126, 0, 3, 804, 1736, 426, 444, 1146,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2186, 0, 3, 822, 1766, 444, 462, 1182,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2246, 0, 3, 840, 1796, 462, 480, 1218,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2306, 0, 3, 894, 1856, 516, 534, 1326,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2366, 0, 3, 912, 1886, 534, 552, 1362,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2426, 0, 3, 930, 1916, 552, 570, 1398,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2486, 0, 3, 948, 1946, 570, 588, 1434,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2546, 0, 3, 966, 1976, 588, 606, 1470,
                                                 ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2606, 3, 642, 648, 1506, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2621, 3, 648, 654, 1516, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2636, 3, 654, 660, 1526, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2651, 3, 660, 666, 1536, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2666, 3, 666, 672, 1546, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2681, 3, 672, 678, 1556, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2696, 3, 678, 684, 1566, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2711, 3, 696, 702, 1576, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2726, 3, 702, 708, 1586, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2741, 3, 708, 714, 1596, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2756, 3, 714, 720, 1606, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2771, 3, 720, 726, 1616, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2786, 3, 726, 732, 1626, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 2801, 3, 732, 738, 1636, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 2816, 0, 3, 1516, 2636, 1676, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2861, 0, 3, 1526, 2651, 1706, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2906, 0, 3, 1536, 2666, 1736, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2951, 0, 3, 1546, 2681, 1766, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 2996, 0, 3, 1556, 2696, 1796, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3041, 0, 3, 1586, 2741, 1856, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3086, 0, 3, 1596, 2756, 1886, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3131, 0, 3, 1606, 2771, 1916, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3176, 0, 3, 1616, 2786, 1946, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 3221, 0, 3, 1626, 2801, 1976, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 3266, 0, 3, 1646, 2816, 1002, 1038,
                                                 2006, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3356, 0, 3, 1676, 2861, 1038, 1074,
                                                 2066, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3446, 0, 3, 1706, 2906, 1074, 1110,
                                                 2126, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3536, 0, 3, 1736, 2951, 1110, 1146,
                                                 2186, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3626, 0, 3, 1766, 2996, 1146, 1182,
                                                 2246, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3716, 0, 3, 1826, 3041, 1254, 1290,
                                                 2306, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3806, 0, 3, 1856, 3086, 1290, 1326,
                                                 2366, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3896, 0, 3, 1886, 3131, 1326, 1362,
                                                 2426, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 3986, 0, 3, 1916, 3176, 1362, 1398,
                                                 2486, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 4076, 0, 3, 1946, 3221, 1398, 1434,
                                                 2546, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4166, 3, 1506, 1516, 2636, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4187, 3, 1516, 1526, 2651, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4208, 3, 1526, 1536, 2666, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4229, 3, 1536, 1546, 2681, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4250, 3, 1546, 1556, 2696, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4271, 3, 1576, 1586, 2741, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4292, 3, 1586, 1596, 2756, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4313, 3, 1596, 1606, 2771, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4334, 3, 1606, 1616, 2786, ncols, alpha,
                                                 beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 4355, 3, 1616, 1626, 2801, ncols, alpha,
                                                 beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 4376, 0, 3, 2636, 4187, 2861, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4439, 0, 3, 2651, 4208, 2906, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4502, 0, 3, 2666, 4229, 2951, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4565, 0, 3, 2681, 4250, 2996, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4628, 0, 3, 2741, 4292, 3086, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4691, 0, 3, 2756, 4313, 3131, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4754, 0, 3, 2771, 4334, 3176, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 4817, 0, 3, 2786, 4355, 3221, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 4880, 0, 3, 2861, 4439, 2006, 2066,
                                                 3446, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5006, 0, 3, 2906, 4502, 2066, 2126,
                                                 3536, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5132, 0, 3, 2951, 4565, 2126, 2186,
                                                 3626, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5258, 0, 3, 3086, 4691, 2306, 2366,
                                                 3896, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5384, 0, 3, 3131, 4754, 2366, 2426,
                                                 3986, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 5510, 0, 3, 3176, 4817, 2426, 2486,
                                                 4076, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5636, 3, 2606, 2621, 4166, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5664, 3, 2621, 2636, 4187, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5692, 3, 2636, 2651, 4208, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5720, 3, 2651, 2666, 4229, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5748, 3, 2666, 2681, 4250, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5776, 3, 2711, 2726, 4271, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5804, 3, 2726, 2741, 4292, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5832, 3, 2741, 2756, 4313, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5860, 3, 2756, 2771, 4334, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 5888, 3, 2771, 2786, 4355, ncols, alpha,
                                                 beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 5916, 0, 3, 4187, 5692, 4439, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 6000, 0, 3, 4208, 5720, 4502, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 6084, 0, 3, 4229, 5748, 4565, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 6168, 0, 3, 4292, 5832, 4691, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 6252, 0, 3, 4313, 5860, 4754, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 6336, 0, 3, 4334, 5888, 4817, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 6420, 0, 3, 4376, 5916, 3266, 3356,
                                                 4880, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6588, 0, 3, 4439, 6000, 3356, 3446,
                                                 5006, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6756, 0, 3, 4502, 6084, 3446, 3536,
                                                 5132, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 6924, 0, 3, 4628, 6168, 3716, 3806,
                                                 5258, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 7092, 0, 3, 4691, 6252, 3806, 3896,
                                                 5384, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 7260, 0, 3, 4754, 6336, 3896, 3986,
                                                 5510, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7428, 3, 4166, 4187, 5692, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7464, 3, 4187, 4208, 5720, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7500, 3, 4208, 4229, 5748, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7536, 3, 4271, 4292, 5832, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7572, 3, 4292, 4313, 5860, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 7608, 3, 4313, 4334, 5888, ncols, alpha,
                                                 beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 7644, 0, 3, 5692, 7464, 6000, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7752, 0, 3, 5720, 7500, 6084, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7860, 0, 3, 5832, 7572, 6252, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 7968, 0, 3, 5860, 7608, 6336, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 8076, 0, 3, 6000, 7752, 4880, 5006,
                                                 6756, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 8292, 0, 3, 6252, 7968, 5258, 5384,
                                                 7260, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8508, 3, 5636, 5664, 7428, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8553, 3, 5692, 5720, 7500, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8598, 3, 5776, 5804, 7536, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 8643, 3, 5832, 5860, 7608, ncols, alpha,
                                                 beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 8688, 0, 3, 7464, 8553, 7752, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 8823, 0, 3, 7572, 8643, 7968, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 8958, 0, 3, 7644, 8688, 6420, 6588,
                                                 8076, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 9228, 0, 3, 7860, 8823, 6924, 7092,
                                                 8292, ncols, alpha, beta, p);

            simdgeo::geom_p_x(buffer, 9498, 8598, 9228, 1, 45, ncols, alpha);

            simdgeo::geom_p_y(buffer, 9633, 8598, 9228, 1, 45, ncols, alpha);

            simdgeo::geom_p_z(buffer, 9768, 8598, 9228, 1, 45, ncols, alpha);

            simdgeo::geom_p_x(buffer, 9903, 8508, 8958, 1, 45, ncols, alpha);

            simdgeo::geom_p_y(buffer, 10038, 8508, 8958, 1, 45, ncols, alpha);

            simdgeo::geom_p_z(buffer, 10173, 8508, 8958, 1, 45, ncols, alpha);

            simdfunc::contract_primitives(buffer, 10308, 9903, 405, ncols);

            simdfunc::contract_primitives(buffer, 10713, 9498, 405, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 11118, 10713, 3, 1, nmax);

    simdtrf::transform_p_outer(values, nvalues, buffer, 11118, 17, nmax);

    simdtrf::transform_l_inner(buffer, 11118, 10848, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 51 * nvalues, nvalues, buffer, 11118, 17, nmax);

    simdtrf::transform_l_inner(buffer, 11118, 10983, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 102 * nvalues, nvalues, buffer, 11118, 17, nmax);

    simdtrf::transform_l_inner(buffer, 11118, 10308, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 153 * nvalues, nvalues, buffer, 11118, 17, nmax);

    simdtrf::transform_l_inner(buffer, 11118, 10443, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 204 * nvalues, nvalues, buffer, 11118, 17, nmax);

    simdtrf::transform_l_inner(buffer, 11118, 10578, 3, 1, nmax);

    simdtrf::transform_p_outer(values + 255 * nvalues, nvalues, buffer, 11118, 17, nmax);
}

}  // namespace simdt2ceri
