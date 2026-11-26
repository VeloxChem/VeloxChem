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


#include "SimdElectronRepulsionGeom10RsRecIG.hpp"

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
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryI1.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ig_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ig_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 31014, 28242, 2520, nvalues);

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

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8,
                                                9, 10, 11}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 18, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 84, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 87, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 90, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 93, 0, 29, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 96, 0, 7, 8, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 102, 0, 8, 9, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 108, 0, 9, 10, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 114, 0, 10, 11, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 120, 0, 11, 12, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 126, 0, 12, 13, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 13, 14, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 14, 15, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 15, 16, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 19, 20, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 20, 21, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 21, 22, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 22, 23, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 23, 24, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 24, 25, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 25, 26, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 192, 0, 26, 27, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 198, 0, 27, 28, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 204, 0, 30, 33, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 214, 0, 33, 36, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 224, 0, 36, 39, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 234, 0, 39, 42, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 244, 0, 42, 45, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 254, 0, 45, 48, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 264, 0, 48, 51, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 274, 0, 51, 54, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 284, 0, 54, 57, 144, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 294, 0, 63, 66, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 304, 0, 66, 69, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 314, 0, 69, 72, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 324, 0, 72, 75, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 334, 0, 75, 78, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 344, 0, 78, 81, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 354, 0, 81, 84, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 364, 0, 84, 87, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 374, 0, 87, 90, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 384, 0, 96, 102, 224, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 399, 0, 102, 108, 234, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 414, 0, 108, 114, 244, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 429, 0, 114, 120, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 444, 0, 120, 126, 264, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 459, 0, 126, 132, 274, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 474, 0, 132, 138, 284, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 489, 0, 150, 156, 314, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 504, 0, 156, 162, 324, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 519, 0, 162, 168, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 534, 0, 168, 174, 344, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 549, 0, 174, 180, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 564, 0, 180, 186, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 579, 0, 186, 192, 374, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 594, 0, 204, 214, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 615, 0, 214, 224, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 636, 0, 224, 234, 414, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 657, 0, 234, 244, 429, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 678, 0, 244, 254, 444, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 699, 0, 254, 264, 459, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 720, 0, 264, 274, 474, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 741, 0, 294, 304, 489, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 762, 0, 304, 314, 504, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 783, 0, 314, 324, 519, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 804, 0, 324, 334, 534, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 825, 0, 334, 344, 549, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 846, 0, 344, 354, 564, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 867, 0, 354, 364, 579, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 888, 0, 384, 399, 636, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 916, 0, 399, 414, 657, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 944, 0, 414, 429, 678, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 972, 0, 429, 444, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1000, 0, 444, 459, 720, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1028, 0, 489, 504, 783, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1056, 0, 504, 519, 804, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1084, 0, 519, 534, 825, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1112, 0, 534, 549, 846, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1140, 0, 549, 564, 867, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1168, 0, 594, 615, 888, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1204, 0, 615, 636, 916, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1240, 0, 636, 657, 944, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1276, 0, 657, 678, 972, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1312, 0, 678, 699, 1000, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1348, 0, 741, 762, 1028, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1384, 0, 762, 783, 1056, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1420, 0, 783, 804, 1084, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1456, 0, 804, 825, 1112, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1492, 0, 825, 846, 1140, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1528, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1531, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1534, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1537, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1540, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1543, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1546, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1549, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1552, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1555, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1558, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1561, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1564, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1567, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1570, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1573, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1576, 3, 9, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1585, 3, 10, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1594, 3, 11, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1603, 3, 12, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1612, 3, 13, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1621, 3, 14, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1630, 3, 15, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1639, 3, 16, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1648, 3, 21, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1657, 3, 22, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1666, 3, 23, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1675, 3, 24, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1684, 3, 25, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1693, 3, 26, 87, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1702, 3, 27, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1711, 3, 28, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1720, 0, 3, 36, 1576, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1738, 0, 3, 39, 1585, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1756, 0, 3, 42, 1594, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1774, 0, 3, 45, 1603, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1792, 0, 3, 48, 1612, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1810, 0, 3, 51, 1621, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1828, 0, 3, 54, 1630, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1846, 0, 3, 57, 1639, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1864, 0, 3, 69, 1648, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1882, 0, 3, 72, 1657, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1900, 0, 3, 75, 1666, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1918, 0, 3, 78, 1675, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1936, 0, 3, 81, 1684, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1954, 0, 3, 84, 1693, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1972, 0, 3, 87, 1702, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1990, 0, 3, 90, 1711, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2008, 0, 3, 102, 1738, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2038, 0, 3, 108, 1756, 234, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2068, 0, 3, 114, 1774, 244, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2098, 0, 3, 120, 1792, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2128, 0, 3, 126, 1810, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2158, 0, 3, 132, 1828, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2188, 0, 3, 138, 1846, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2218, 0, 3, 156, 1882, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2248, 0, 3, 162, 1900, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2278, 0, 3, 168, 1918, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2308, 0, 3, 174, 1936, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2338, 0, 3, 180, 1954, 354, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2368, 0, 3, 186, 1972, 364, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2398, 0, 3, 192, 1990, 374, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2428, 0, 3, 224, 2038, 399, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2473, 0, 3, 234, 2068, 414, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2518, 0, 3, 244, 2098, 429, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2563, 0, 3, 254, 2128, 444, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2608, 0, 3, 264, 2158, 459, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2653, 0, 3, 274, 2188, 474, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2698, 0, 3, 314, 2248, 504, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2743, 0, 3, 324, 2278, 519, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2788, 0, 3, 334, 2308, 534, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2833, 0, 3, 344, 2338, 549, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2878, 0, 3, 354, 2368, 564, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2923, 0, 3, 364, 2398, 579, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2968, 0, 3, 399, 2473, 636, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3031, 0, 3, 414, 2518, 657, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3094, 0, 3, 429, 2563, 678, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3157, 0, 3, 444, 2608, 699, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3220, 0, 3, 459, 2653, 720, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3283, 0, 3, 504, 2743, 783, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3346, 0, 3, 519, 2788, 804, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3409, 0, 3, 534, 2833, 825, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3472, 0, 3, 549, 2878, 846, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3535, 0, 3, 564, 2923, 867, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3598, 0, 3, 636, 3031, 916, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3682, 0, 3, 657, 3094, 944, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3766, 0, 3, 678, 3157, 972, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3850, 0, 3, 699, 3220, 1000, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3934, 0, 3, 783, 3346, 1056, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4018, 0, 3, 804, 3409, 1084, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4102, 0, 3, 825, 3472, 1112, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4186, 0, 3, 846, 3535, 1140, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4270, 0, 3, 916, 3682, 1240, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4378, 0, 3, 944, 3766, 1276, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4486, 0, 3, 972, 3850, 1312, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4594, 0, 3, 1056, 4018, 1420, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 4702, 0, 3, 1084, 4102, 1456, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 4810, 0, 3, 1112, 4186, 1492, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 4918, 3, 9, 10, 1531, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4924, 3, 10, 11, 1534, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4930, 3, 11, 12, 1537, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4936, 3, 12, 13, 1540, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4942, 3, 13, 14, 1543, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4948, 3, 14, 15, 1546, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4954, 3, 15, 16, 1549, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4960, 3, 21, 22, 1555, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4966, 3, 22, 23, 1558, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4972, 3, 23, 24, 1561, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4978, 3, 24, 25, 1564, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4984, 3, 25, 26, 1567, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4990, 3, 26, 27, 1570, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4996, 3, 27, 28, 1573, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 5002, 0, 3, 1528, 4918, 1585, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5020, 0, 3, 1531, 4924, 1594, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5038, 0, 3, 1534, 4930, 1603, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5056, 0, 3, 1537, 4936, 1612, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5074, 0, 3, 1540, 4942, 1621, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5092, 0, 3, 1543, 4948, 1630, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5110, 0, 3, 1546, 4954, 1639, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5128, 0, 3, 1552, 4960, 1657, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5146, 0, 3, 1555, 4966, 1666, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5164, 0, 3, 1558, 4972, 1675, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5182, 0, 3, 1561, 4978, 1684, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5200, 0, 3, 1564, 4984, 1693, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5218, 0, 3, 1567, 4990, 1702, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5236, 0, 3, 1570, 4996, 1711, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 5254, 0, 3, 1576, 5002, 96, 102, 1738,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5290, 0, 3, 1585, 5020, 102, 108, 1756,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5326, 0, 3, 1594, 5038, 108, 114, 1774,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5362, 0, 3, 1603, 5056, 114, 120, 1792,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5398, 0, 3, 1612, 5074, 120, 126, 1810,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5434, 0, 3, 1621, 5092, 126, 132, 1828,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5470, 0, 3, 1630, 5110, 132, 138, 1846,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5506, 0, 3, 1648, 5128, 150, 156, 1882,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5542, 0, 3, 1657, 5146, 156, 162, 1900,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5578, 0, 3, 1666, 5164, 162, 168, 1918,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5614, 0, 3, 1675, 5182, 168, 174, 1936,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5650, 0, 3, 1684, 5200, 174, 180, 1954,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5686, 0, 3, 1693, 5218, 180, 186, 1972,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5722, 0, 3, 1702, 5236, 186, 192, 1990,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5758, 0, 3, 1720, 5254, 204, 214, 2008,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5818, 0, 3, 1738, 5290, 214, 224, 2038,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5878, 0, 3, 1756, 5326, 224, 234, 2068,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5938, 0, 3, 1774, 5362, 234, 244, 2098,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5998, 0, 3, 1792, 5398, 244, 254, 2128,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6058, 0, 3, 1810, 5434, 254, 264, 2158,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6118, 0, 3, 1828, 5470, 264, 274, 2188,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6178, 0, 3, 1864, 5506, 294, 304, 2218,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6238, 0, 3, 1882, 5542, 304, 314, 2248,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6298, 0, 3, 1900, 5578, 314, 324, 2278,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6358, 0, 3, 1918, 5614, 324, 334, 2308,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6418, 0, 3, 1936, 5650, 334, 344, 2338,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6478, 0, 3, 1954, 5686, 344, 354, 2368,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6538, 0, 3, 1972, 5722, 354, 364, 2398,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6598, 0, 3, 5254, 5290, 2038, 5878, 384,
                                                 399, 2473, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6688, 0, 3, 5290, 5326, 2068, 5938, 399,
                                                 414, 2518, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6778, 0, 3, 5326, 5362, 2098, 5998, 414,
                                                 429, 2563, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6868, 0, 3, 5362, 5398, 2128, 6058, 429,
                                                 444, 2608, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6958, 0, 3, 5398, 5434, 2158, 6118, 444,
                                                 459, 2653, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7048, 0, 3, 5506, 5542, 2248, 6298, 489,
                                                 504, 2743, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7138, 0, 3, 5542, 5578, 2278, 6358, 504,
                                                 519, 2788, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7228, 0, 3, 5578, 5614, 2308, 6418, 519,
                                                 534, 2833, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7318, 0, 3, 5614, 5650, 2338, 6478, 534,
                                                 549, 2878, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7408, 0, 3, 5650, 5686, 2368, 6538, 549,
                                                 564, 2923, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7498, 0, 3, 5758, 5818, 2428, 6598, 594,
                                                 615, 2968, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7624, 0, 3, 5818, 5878, 2473, 6688, 615,
                                                 636, 3031, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7750, 0, 3, 5878, 5938, 2518, 6778, 636,
                                                 657, 3094, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7876, 0, 3, 5938, 5998, 2563, 6868, 657,
                                                 678, 3157, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8002, 0, 3, 5998, 6058, 2608, 6958, 678,
                                                 699, 3220, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8128, 0, 3, 6178, 6238, 2698, 7048, 741,
                                                 762, 3283, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8254, 0, 3, 6238, 6298, 2743, 7138, 762,
                                                 783, 3346, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8380, 0, 3, 6298, 6358, 2788, 7228, 783,
                                                 804, 3409, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8506, 0, 3, 6358, 6418, 2833, 7318, 804,
                                                 825, 3472, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8632, 0, 3, 6418, 6478, 2878, 7408, 825,
                                                 846, 3535, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8758, 0, 3, 6598, 6688, 3031, 7750, 888,
                                                 916, 3682, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8926, 0, 3, 6688, 6778, 3094, 7876, 916,
                                                 944, 3766, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9094, 0, 3, 6778, 6868, 3157, 8002, 944,
                                                 972, 3850, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9262, 0, 3, 7048, 7138, 3346, 8380,
                                                 1028, 1056, 4018, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9430, 0, 3, 7138, 7228, 3409, 8506,
                                                 1056, 1084, 4102, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9598, 0, 3, 7228, 7318, 3472, 8632,
                                                 1084, 1112, 4186, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9766, 0, 3, 7498, 7624, 3598, 8758,
                                                 1168, 1204, 4270, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9982, 0, 3, 7624, 7750, 3682, 8926,
                                                 1204, 1240, 4378, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10198, 0, 3, 7750, 7876, 3766, 9094,
                                                 1240, 1276, 4486, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10414, 0, 3, 8128, 8254, 3934, 9262,
                                                 1348, 1384, 4594, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10630, 0, 3, 8254, 8380, 4018, 9430,
                                                 1384, 1420, 4702, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10846, 0, 3, 8380, 8506, 4102, 9598,
                                                 1420, 1456, 4810, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11062, 3, 1528, 1531, 4924, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11072, 3, 1531, 1534, 4930, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11082, 3, 1534, 1537, 4936, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11092, 3, 1537, 1540, 4942, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11102, 3, 1540, 1543, 4948, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11112, 3, 1543, 1546, 4954, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11122, 3, 1552, 1555, 4966, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11132, 3, 1555, 1558, 4972, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11142, 3, 1558, 1561, 4978, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11152, 3, 1561, 1564, 4984, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11162, 3, 1564, 1567, 4990, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11172, 3, 1567, 1570, 4996, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 11182, 0, 3, 4918, 11062, 5020, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11212, 0, 3, 4924, 11072, 5038, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11242, 0, 3, 4930, 11082, 5056, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11272, 0, 3, 4936, 11092, 5074, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11302, 0, 3, 4942, 11102, 5092, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11332, 0, 3, 4948, 11112, 5110, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11362, 0, 3, 4960, 11122, 5146, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11392, 0, 3, 4966, 11132, 5164, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11422, 0, 3, 4972, 11142, 5182, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11452, 0, 3, 4978, 11152, 5200, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11482, 0, 3, 4984, 11162, 5218, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11512, 0, 3, 4990, 11172, 5236, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 11542, 0, 3, 5002, 11182, 1720, 1738,
                                                 5290, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11602, 0, 3, 5020, 11212, 1738, 1756,
                                                 5326, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11662, 0, 3, 5038, 11242, 1756, 1774,
                                                 5362, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11722, 0, 3, 5056, 11272, 1774, 1792,
                                                 5398, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11782, 0, 3, 5074, 11302, 1792, 1810,
                                                 5434, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11842, 0, 3, 5092, 11332, 1810, 1828,
                                                 5470, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11902, 0, 3, 5128, 11362, 1864, 1882,
                                                 5542, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11962, 0, 3, 5146, 11392, 1882, 1900,
                                                 5578, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12022, 0, 3, 5164, 11422, 1900, 1918,
                                                 5614, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12082, 0, 3, 5182, 11452, 1918, 1936,
                                                 5650, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12142, 0, 3, 5200, 11482, 1936, 1954,
                                                 5686, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12202, 0, 3, 5218, 11512, 1954, 1972,
                                                 5722, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12262, 0, 3, 5290, 11602, 2008, 2038,
                                                 5878, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12362, 0, 3, 5326, 11662, 2038, 2068,
                                                 5938, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12462, 0, 3, 5362, 11722, 2068, 2098,
                                                 5998, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12562, 0, 3, 5398, 11782, 2098, 2128,
                                                 6058, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12662, 0, 3, 5434, 11842, 2128, 2158,
                                                 6118, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12762, 0, 3, 5542, 11962, 2218, 2248,
                                                 6298, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12862, 0, 3, 5578, 12022, 2248, 2278,
                                                 6358, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12962, 0, 3, 5614, 12082, 2278, 2308,
                                                 6418, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13062, 0, 3, 5650, 12142, 2308, 2338,
                                                 6478, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13162, 0, 3, 5686, 12202, 2338, 2368,
                                                 6538, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13262, 0, 3, 11542, 11602, 5878, 12362,
                                                 2428, 2473, 6688, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13412, 0, 3, 11602, 11662, 5938, 12462,
                                                 2473, 2518, 6778, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13562, 0, 3, 11662, 11722, 5998, 12562,
                                                 2518, 2563, 6868, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13712, 0, 3, 11722, 11782, 6058, 12662,
                                                 2563, 2608, 6958, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13862, 0, 3, 11902, 11962, 6298, 12862,
                                                 2698, 2743, 7138, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14012, 0, 3, 11962, 12022, 6358, 12962,
                                                 2743, 2788, 7228, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14162, 0, 3, 12022, 12082, 6418, 13062,
                                                 2788, 2833, 7318, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14312, 0, 3, 12082, 12142, 6478, 13162,
                                                 2833, 2878, 7408, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14462, 0, 3, 12262, 12362, 6688, 13412,
                                                 2968, 3031, 7750, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14672, 0, 3, 12362, 12462, 6778, 13562,
                                                 3031, 3094, 7876, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14882, 0, 3, 12462, 12562, 6868, 13712,
                                                 3094, 3157, 8002, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15092, 0, 3, 12762, 12862, 7138, 14012,
                                                 3283, 3346, 8380, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15302, 0, 3, 12862, 12962, 7228, 14162,
                                                 3346, 3409, 8506, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15512, 0, 3, 12962, 13062, 7318, 14312,
                                                 3409, 3472, 8632, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15722, 0, 3, 13262, 13412, 7750, 14672,
                                                 3598, 3682, 8926, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16002, 0, 3, 13412, 13562, 7876, 14882,
                                                 3682, 3766, 9094, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16282, 0, 3, 13862, 14012, 8380, 15302,
                                                 3934, 4018, 9430, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16562, 0, 3, 14012, 14162, 8506, 15512,
                                                 4018, 4102, 9598, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 16842, 0, 3, 14462, 14672, 8926, 16002,
                                                 4270, 4378, 10198, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 17202, 0, 3, 15092, 15302, 9430, 16562,
                                                 4594, 4702, 10846, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17562, 3, 4918, 4924, 11072, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17577, 3, 4924, 4930, 11082, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17592, 3, 4930, 4936, 11092, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17607, 3, 4936, 4942, 11102, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17622, 3, 4942, 4948, 11112, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17637, 3, 4960, 4966, 11132, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17652, 3, 4966, 4972, 11142, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17667, 3, 4972, 4978, 11152, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17682, 3, 4978, 4984, 11162, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 17697, 3, 4984, 4990, 11172, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 17712, 0, 3, 11062, 17562, 11212, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17757, 0, 3, 11072, 17577, 11242, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17802, 0, 3, 11082, 17592, 11272, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17847, 0, 3, 11092, 17607, 11302, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17892, 0, 3, 11102, 17622, 11332, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17937, 0, 3, 11122, 17637, 11392, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 17982, 0, 3, 11132, 17652, 11422, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18027, 0, 3, 11142, 17667, 11452, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18072, 0, 3, 11152, 17682, 11482, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 18117, 0, 3, 11162, 17697, 11512, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 18162, 0, 3, 11182, 17712, 5254, 5290,
                                                 11602, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18252, 0, 3, 11212, 17757, 5290, 5326,
                                                 11662, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18342, 0, 3, 11242, 17802, 5326, 5362,
                                                 11722, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18432, 0, 3, 11272, 17847, 5362, 5398,
                                                 11782, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18522, 0, 3, 11302, 17892, 5398, 5434,
                                                 11842, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18612, 0, 3, 11362, 17937, 5506, 5542,
                                                 11962, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18702, 0, 3, 11392, 17982, 5542, 5578,
                                                 12022, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18792, 0, 3, 11422, 18027, 5578, 5614,
                                                 12082, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18882, 0, 3, 11452, 18072, 5614, 5650,
                                                 12142, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 18972, 0, 3, 11482, 18117, 5650, 5686,
                                                 12202, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19062, 0, 3, 11542, 18162, 5758, 5818,
                                                 12262, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19212, 0, 3, 11602, 18252, 5818, 5878,
                                                 12362, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19362, 0, 3, 11662, 18342, 5878, 5938,
                                                 12462, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19512, 0, 3, 11722, 18432, 5938, 5998,
                                                 12562, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19662, 0, 3, 11782, 18522, 5998, 6058,
                                                 12662, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19812, 0, 3, 11902, 18612, 6178, 6238,
                                                 12762, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19962, 0, 3, 11962, 18702, 6238, 6298,
                                                 12862, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20112, 0, 3, 12022, 18792, 6298, 6358,
                                                 12962, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20262, 0, 3, 12082, 18882, 6358, 6418,
                                                 13062, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20412, 0, 3, 12142, 18972, 6418, 6478,
                                                 13162, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20562, 0, 3, 18162, 18252, 12362, 19362,
                                                 6598, 6688, 13412, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20787, 0, 3, 18252, 18342, 12462, 19512,
                                                 6688, 6778, 13562, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21012, 0, 3, 18342, 18432, 12562, 19662,
                                                 6778, 6868, 13712, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21237, 0, 3, 18612, 18702, 12862, 20112,
                                                 7048, 7138, 14012, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21462, 0, 3, 18702, 18792, 12962, 20262,
                                                 7138, 7228, 14162, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21687, 0, 3, 18792, 18882, 13062, 20412,
                                                 7228, 7318, 14312, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 21912, 0, 3, 19062, 19212, 13262, 20562,
                                                 7498, 7624, 14462, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22227, 0, 3, 19212, 19362, 13412, 20787,
                                                 7624, 7750, 14672, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22542, 0, 3, 19362, 19512, 13562, 21012,
                                                 7750, 7876, 14882, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22857, 0, 3, 19812, 19962, 13862, 21237,
                                                 8128, 8254, 15092, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23172, 0, 3, 19962, 20112, 14012, 21462,
                                                 8254, 8380, 15302, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23487, 0, 3, 20112, 20262, 14162, 21687,
                                                 8380, 8506, 15512, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 23802, 0, 3, 20562, 20787, 14672, 22542,
                                                 8758, 8926, 16002, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 24222, 0, 3, 21237, 21462, 15302, 23487,
                                                 9262, 9430, 16562, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 24642, 0, 3, 21912, 22227, 15722, 23802,
                                                 9766, 9982, 16842, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 25182, 0, 3, 22857, 23172, 16282, 24222,
                                                 10414, 10630, 17202, ncols, alpha, beta, p);

            simdgeo::geom_i_x(buffer, 25722, 22857, 25182, 1, 15, ncols, alpha);

            simdgeo::geom_i_y(buffer, 26142, 22857, 25182, 1, 15, ncols, alpha);

            simdgeo::geom_i_z(buffer, 26562, 22857, 25182, 1, 15, ncols, alpha);

            simdgeo::geom_i_x(buffer, 26982, 21912, 24642, 1, 15, ncols, alpha);

            simdgeo::geom_i_y(buffer, 27402, 21912, 24642, 1, 15, ncols, alpha);

            simdgeo::geom_i_z(buffer, 27822, 21912, 24642, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 28242, 26982, 1260, ncols);

            simdfunc::contract_primitives(buffer, 29502, 25722, 1260, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 30762, 29502, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 30762, 9, nmax);

    simdtrf::transform_g_inner(buffer, 30762, 29922, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 117 * nvalues, nvalues, buffer, 30762, 9, nmax);

    simdtrf::transform_g_inner(buffer, 30762, 30342, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 234 * nvalues, nvalues, buffer, 30762, 9, nmax);

    simdtrf::transform_g_inner(buffer, 30762, 28242, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 351 * nvalues, nvalues, buffer, 30762, 9, nmax);

    simdtrf::transform_g_inner(buffer, 30762, 28662, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 468 * nvalues, nvalues, buffer, 30762, 9, nmax);

    simdtrf::transform_g_inner(buffer, 30762, 29082, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 585 * nvalues, nvalues, buffer, 30762, 9, nmax);
}

}  // namespace simdt2ceri
