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


#include "SimdElectronRepulsionGeom10RsRecLD.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryL1.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ld_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ld_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 17197, 15352, 1620, nvalues);

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

            compute_prim_ls_electron_repulsion_0(buffer, 1528, 0, 888, 916, 1240, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1573, 0, 916, 944, 1276, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1618, 0, 944, 972, 1312, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1663, 0, 1028, 1056, 1420, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1708, 0, 1056, 1084, 1456, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1753, 0, 1084, 1112, 1492, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1798, 0, 1168, 1204, 1528, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1853, 0, 1204, 1240, 1573, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1908, 0, 1240, 1276, 1618, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1963, 0, 1348, 1384, 1663, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2018, 0, 1384, 1420, 1708, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2073, 0, 1420, 1456, 1753, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 2128, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2131, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2134, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2137, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2140, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2143, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2146, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2149, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2152, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2155, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2158, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2161, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2164, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2167, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2170, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2173, 3, 29, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2176, 3, 9, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2185, 3, 10, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2194, 3, 11, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2203, 3, 12, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2212, 3, 13, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2221, 3, 14, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2230, 3, 15, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2239, 3, 16, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2248, 3, 21, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2257, 3, 22, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2266, 3, 23, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2275, 3, 24, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2284, 3, 25, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2293, 3, 26, 87, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2302, 3, 27, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2311, 3, 28, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2320, 0, 3, 36, 2176, 102, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2338, 0, 3, 39, 2185, 108, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2356, 0, 3, 42, 2194, 114, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2374, 0, 3, 45, 2203, 120, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2392, 0, 3, 48, 2212, 126, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2410, 0, 3, 51, 2221, 132, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2428, 0, 3, 54, 2230, 138, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2446, 0, 3, 57, 2239, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2464, 0, 3, 69, 2248, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2482, 0, 3, 72, 2257, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2500, 0, 3, 75, 2266, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2518, 0, 3, 78, 2275, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2536, 0, 3, 81, 2284, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2554, 0, 3, 84, 2293, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2572, 0, 3, 87, 2302, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2590, 0, 3, 90, 2311, 198, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2608, 0, 3, 102, 2338, 224, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2638, 0, 3, 108, 2356, 234, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2668, 0, 3, 114, 2374, 244, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2698, 0, 3, 120, 2392, 254, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2728, 0, 3, 126, 2410, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2758, 0, 3, 132, 2428, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2788, 0, 3, 138, 2446, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2818, 0, 3, 156, 2482, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2848, 0, 3, 162, 2500, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2878, 0, 3, 168, 2518, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2908, 0, 3, 174, 2536, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2938, 0, 3, 180, 2554, 354, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2968, 0, 3, 186, 2572, 364, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2998, 0, 3, 192, 2590, 374, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3028, 0, 3, 224, 2638, 399, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3073, 0, 3, 234, 2668, 414, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3118, 0, 3, 244, 2698, 429, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3163, 0, 3, 254, 2728, 444, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3208, 0, 3, 264, 2758, 459, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3253, 0, 3, 274, 2788, 474, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3298, 0, 3, 314, 2848, 504, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3343, 0, 3, 324, 2878, 519, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3388, 0, 3, 334, 2908, 534, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3433, 0, 3, 344, 2938, 549, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3478, 0, 3, 354, 2968, 564, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3523, 0, 3, 364, 2998, 579, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3568, 0, 3, 399, 3073, 636, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3631, 0, 3, 414, 3118, 657, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3694, 0, 3, 429, 3163, 678, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3757, 0, 3, 444, 3208, 699, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3820, 0, 3, 459, 3253, 720, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3883, 0, 3, 504, 3343, 783, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3946, 0, 3, 519, 3388, 804, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4009, 0, 3, 534, 3433, 825, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4072, 0, 3, 549, 3478, 846, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4135, 0, 3, 564, 3523, 867, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4198, 0, 3, 636, 3631, 916, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4282, 0, 3, 657, 3694, 944, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4366, 0, 3, 678, 3757, 972, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4450, 0, 3, 699, 3820, 1000, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4534, 0, 3, 783, 3946, 1056, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4618, 0, 3, 804, 4009, 1084, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4702, 0, 3, 825, 4072, 1112, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4786, 0, 3, 846, 4135, 1140, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4870, 0, 3, 916, 4282, 1240, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4978, 0, 3, 944, 4366, 1276, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5086, 0, 3, 972, 4450, 1312, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5194, 0, 3, 1056, 4618, 1420, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5302, 0, 3, 1084, 4702, 1456, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5410, 0, 3, 1112, 4786, 1492, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5518, 0, 3, 1240, 4978, 1573, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5653, 0, 3, 1276, 5086, 1618, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5788, 0, 3, 1420, 5302, 1708, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5923, 0, 3, 1456, 5410, 1753, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 6058, 0, 3, 1573, 5653, 1908, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 6223, 0, 3, 1708, 5923, 2073, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 6388, 3, 9, 10, 2131, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6394, 3, 10, 11, 2134, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6400, 3, 11, 12, 2137, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6406, 3, 12, 13, 2140, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6412, 3, 13, 14, 2143, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6418, 3, 14, 15, 2146, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6424, 3, 15, 16, 2149, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6430, 3, 21, 22, 2155, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6436, 3, 22, 23, 2158, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6442, 3, 23, 24, 2161, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6448, 3, 24, 25, 2164, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6454, 3, 25, 26, 2167, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6460, 3, 26, 27, 2170, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6466, 3, 27, 28, 2173, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 6472, 0, 3, 2128, 6388, 2185, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6490, 0, 3, 2131, 6394, 2194, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6508, 0, 3, 2134, 6400, 2203, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6526, 0, 3, 2137, 6406, 2212, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6544, 0, 3, 2140, 6412, 2221, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6562, 0, 3, 2143, 6418, 2230, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6580, 0, 3, 2146, 6424, 2239, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6598, 0, 3, 2152, 6430, 2257, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6616, 0, 3, 2155, 6436, 2266, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6634, 0, 3, 2158, 6442, 2275, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6652, 0, 3, 2161, 6448, 2284, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6670, 0, 3, 2164, 6454, 2293, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6688, 0, 3, 2167, 6460, 2302, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6706, 0, 3, 2170, 6466, 2311, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 6724, 0, 3, 2176, 6472, 96, 102, 2338,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6760, 0, 3, 2185, 6490, 102, 108, 2356,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6796, 0, 3, 2194, 6508, 108, 114, 2374,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6832, 0, 3, 2203, 6526, 114, 120, 2392,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6868, 0, 3, 2212, 6544, 120, 126, 2410,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6904, 0, 3, 2221, 6562, 126, 132, 2428,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6940, 0, 3, 2230, 6580, 132, 138, 2446,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6976, 0, 3, 2248, 6598, 150, 156, 2482,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7012, 0, 3, 2257, 6616, 156, 162, 2500,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7048, 0, 3, 2266, 6634, 162, 168, 2518,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7084, 0, 3, 2275, 6652, 168, 174, 2536,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7120, 0, 3, 2284, 6670, 174, 180, 2554,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7156, 0, 3, 2293, 6688, 180, 186, 2572,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7192, 0, 3, 2302, 6706, 186, 192, 2590,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7228, 0, 3, 2320, 6724, 204, 214, 2608,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7288, 0, 3, 2338, 6760, 214, 224, 2638,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7348, 0, 3, 2356, 6796, 224, 234, 2668,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7408, 0, 3, 2374, 6832, 234, 244, 2698,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7468, 0, 3, 2392, 6868, 244, 254, 2728,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7528, 0, 3, 2410, 6904, 254, 264, 2758,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7588, 0, 3, 2428, 6940, 264, 274, 2788,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7648, 0, 3, 2464, 6976, 294, 304, 2818,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7708, 0, 3, 2482, 7012, 304, 314, 2848,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7768, 0, 3, 2500, 7048, 314, 324, 2878,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7828, 0, 3, 2518, 7084, 324, 334, 2908,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7888, 0, 3, 2536, 7120, 334, 344, 2938,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7948, 0, 3, 2554, 7156, 344, 354, 2968,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8008, 0, 3, 2572, 7192, 354, 364, 2998,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8068, 0, 3, 6724, 6760, 2638, 7348, 384,
                                                 399, 3073, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8158, 0, 3, 6760, 6796, 2668, 7408, 399,
                                                 414, 3118, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8248, 0, 3, 6796, 6832, 2698, 7468, 414,
                                                 429, 3163, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8338, 0, 3, 6832, 6868, 2728, 7528, 429,
                                                 444, 3208, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8428, 0, 3, 6868, 6904, 2758, 7588, 444,
                                                 459, 3253, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8518, 0, 3, 6976, 7012, 2848, 7768, 489,
                                                 504, 3343, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8608, 0, 3, 7012, 7048, 2878, 7828, 504,
                                                 519, 3388, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8698, 0, 3, 7048, 7084, 2908, 7888, 519,
                                                 534, 3433, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8788, 0, 3, 7084, 7120, 2938, 7948, 534,
                                                 549, 3478, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8878, 0, 3, 7120, 7156, 2968, 8008, 549,
                                                 564, 3523, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8968, 0, 3, 7228, 7288, 3028, 8068, 594,
                                                 615, 3568, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9094, 0, 3, 7288, 7348, 3073, 8158, 615,
                                                 636, 3631, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9220, 0, 3, 7348, 7408, 3118, 8248, 636,
                                                 657, 3694, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9346, 0, 3, 7408, 7468, 3163, 8338, 657,
                                                 678, 3757, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9472, 0, 3, 7468, 7528, 3208, 8428, 678,
                                                 699, 3820, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9598, 0, 3, 7648, 7708, 3298, 8518, 741,
                                                 762, 3883, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9724, 0, 3, 7708, 7768, 3343, 8608, 762,
                                                 783, 3946, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9850, 0, 3, 7768, 7828, 3388, 8698, 783,
                                                 804, 4009, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 9976, 0, 3, 7828, 7888, 3433, 8788, 804,
                                                 825, 4072, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10102, 0, 3, 7888, 7948, 3478, 8878,
                                                 825, 846, 4135, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10228, 0, 3, 8068, 8158, 3631, 9220,
                                                 888, 916, 4282, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10396, 0, 3, 8158, 8248, 3694, 9346,
                                                 916, 944, 4366, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10564, 0, 3, 8248, 8338, 3757, 9472,
                                                 944, 972, 4450, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10732, 0, 3, 8518, 8608, 3946, 9850,
                                                 1028, 1056, 4618, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 10900, 0, 3, 8608, 8698, 4009, 9976,
                                                 1056, 1084, 4702, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11068, 0, 3, 8698, 8788, 4072, 10102,
                                                 1084, 1112, 4786, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11236, 0, 3, 8968, 9094, 4198, 10228,
                                                 1168, 1204, 4870, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11452, 0, 3, 9094, 9220, 4282, 10396,
                                                 1204, 1240, 4978, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11668, 0, 3, 9220, 9346, 4366, 10564,
                                                 1240, 1276, 5086, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11884, 0, 3, 9598, 9724, 4534, 10732,
                                                 1348, 1384, 5194, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12100, 0, 3, 9724, 9850, 4618, 10900,
                                                 1384, 1420, 5302, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 12316, 0, 3, 9850, 9976, 4702, 11068,
                                                 1420, 1456, 5410, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 12532, 0, 3, 10228, 10396, 4978, 11668,
                                                 1528, 1573, 5653, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 12802, 0, 3, 10732, 10900, 5302, 12316,
                                                 1663, 1708, 5923, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 13072, 0, 3, 11236, 11452, 5518, 12532,
                                                 1798, 1853, 6058, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 13402, 0, 3, 11884, 12100, 5788, 12802,
                                                 1963, 2018, 6223, ncols, alpha, beta, p);

            simdgeo::geom_l_x(buffer, 13732, 11884, 13402, 1, 6, ncols, alpha);

            simdgeo::geom_l_y(buffer, 14002, 11884, 13402, 1, 6, ncols, alpha);

            simdgeo::geom_l_z(buffer, 14272, 11884, 13402, 1, 6, ncols, alpha);

            simdgeo::geom_l_x(buffer, 14542, 11236, 13072, 1, 6, ncols, alpha);

            simdgeo::geom_l_y(buffer, 14812, 11236, 13072, 1, 6, ncols, alpha);

            simdgeo::geom_l_z(buffer, 15082, 11236, 13072, 1, 6, ncols, alpha);

            simdfunc::contract_primitives(buffer, 15352, 14542, 810, ncols);

            simdfunc::contract_primitives(buffer, 16162, 13732, 810, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 16972, 16162, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 16972, 5, nmax);

    simdtrf::transform_d_inner(buffer, 16972, 16432, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 85 * nvalues, nvalues, buffer, 16972, 5, nmax);

    simdtrf::transform_d_inner(buffer, 16972, 16702, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 170 * nvalues, nvalues, buffer, 16972, 5, nmax);

    simdtrf::transform_d_inner(buffer, 16972, 15352, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 255 * nvalues, nvalues, buffer, 16972, 5, nmax);

    simdtrf::transform_d_inner(buffer, 16972, 15622, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 340 * nvalues, nvalues, buffer, 16972, 5, nmax);

    simdtrf::transform_d_inner(buffer, 16972, 15892, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 425 * nvalues, nvalues, buffer, 16972, 5, nmax);
}

}  // namespace simdt2ceri
