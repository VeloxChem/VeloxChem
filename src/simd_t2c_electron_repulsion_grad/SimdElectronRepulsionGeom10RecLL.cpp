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


#include "SimdElectronRepulsionGeom10RecLL.hpp"

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
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGL.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHL.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIL.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKL.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLK.hpp"
#include "SimdElectronRepulsionVrrRecLL.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMH.hpp"
#include "SimdElectronRepulsionVrrRecMI.hpp"
#include "SimdElectronRepulsionVrrRecMK.hpp"
#include "SimdElectronRepulsionVrrRecML.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
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
#include "SimdGeometryL1.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_geom_10_ll_electron_repulsion(double               *values,
                                      const size_t          nvalues,
                                      const CBasisFunction &bra,
                                      const CBasisFunction &ket,
                                      const CSimdMatrix    &coordinates,
                                      CSimdMatrix          &buffer) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_10_ll_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 177970, 171130, 6075, nvalues);

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
                                            10, 11, 12, 13, 14, 15, 16, 17}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 24, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 27, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 30, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 33, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 36, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 39, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 23, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 75, 0, 7, 8, 30, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 81, 0, 8, 9, 33, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 87, 0, 9, 10, 36, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 93, 0, 10, 11, 39, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 99, 0, 11, 12, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 105, 0, 12, 13, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 111, 0, 13, 14, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 117, 0, 14, 15, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 123, 0, 15, 16, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 129, 0, 16, 17, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 135, 0, 17, 18, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 141, 0, 18, 19, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 147, 0, 19, 20, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 153, 0, 20, 21, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 159, 0, 21, 22, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 165, 0, 24, 27, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 175, 0, 27, 30, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 185, 0, 30, 33, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 195, 0, 33, 36, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 205, 0, 36, 39, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 215, 0, 39, 42, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 225, 0, 42, 45, 111, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 235, 0, 45, 48, 117, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 245, 0, 48, 51, 123, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 255, 0, 51, 54, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 265, 0, 54, 57, 135, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 275, 0, 57, 60, 141, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 285, 0, 60, 63, 147, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 295, 0, 63, 66, 153, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 305, 0, 66, 69, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 315, 0, 75, 81, 185, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 330, 0, 81, 87, 195, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 345, 0, 87, 93, 205, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 360, 0, 93, 99, 215, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 375, 0, 99, 105, 225, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 390, 0, 105, 111, 235, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 405, 0, 111, 117, 245, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 420, 0, 117, 123, 255, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 435, 0, 123, 129, 265, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 450, 0, 129, 135, 275, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 465, 0, 135, 141, 285, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 480, 0, 141, 147, 295, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 495, 0, 147, 153, 305, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 510, 0, 165, 175, 315, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 531, 0, 175, 185, 330, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 552, 0, 185, 195, 345, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 573, 0, 195, 205, 360, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 594, 0, 205, 215, 375, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 615, 0, 215, 225, 390, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 636, 0, 225, 235, 405, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 657, 0, 235, 245, 420, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 678, 0, 245, 255, 435, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 699, 0, 255, 265, 450, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 720, 0, 265, 275, 465, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 741, 0, 275, 285, 480, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 762, 0, 285, 295, 495, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 783, 0, 315, 330, 552, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 811, 0, 330, 345, 573, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 839, 0, 345, 360, 594, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 867, 0, 360, 375, 615, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 895, 0, 375, 390, 636, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 923, 0, 390, 405, 657, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 951, 0, 405, 420, 678, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 979, 0, 420, 435, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1007, 0, 435, 450, 720, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1035, 0, 450, 465, 741, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1063, 0, 465, 480, 762, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1091, 0, 510, 531, 783, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1127, 0, 531, 552, 811, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1163, 0, 552, 573, 839, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1199, 0, 573, 594, 867, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1235, 0, 594, 615, 895, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1271, 0, 615, 636, 923, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1307, 0, 636, 657, 951, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1343, 0, 657, 678, 979, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1379, 0, 678, 699, 1007, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1415, 0, 699, 720, 1035, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1451, 0, 720, 741, 1063, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1487, 0, 783, 811, 1163, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1532, 0, 811, 839, 1199, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1577, 0, 839, 867, 1235, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1622, 0, 867, 895, 1271, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1667, 0, 895, 923, 1307, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1712, 0, 923, 951, 1343, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1757, 0, 951, 979, 1379, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1802, 0, 979, 1007, 1415, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1847, 0, 1007, 1035, 1451, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1892, 0, 1091, 1127, 1487, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 1947, 0, 1127, 1163, 1532, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2002, 0, 1163, 1199, 1577, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2057, 0, 1199, 1235, 1622, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2112, 0, 1235, 1271, 1667, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2167, 0, 1271, 1307, 1712, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2222, 0, 1307, 1343, 1757, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2277, 0, 1343, 1379, 1802, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2332, 0, 1379, 1415, 1847, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 2387, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2390, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2393, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2396, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2399, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2402, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2405, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2408, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2411, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2414, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2417, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2420, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2423, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2426, 3, 23, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2429, 3, 9, 33, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2438, 3, 10, 36, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2447, 3, 11, 39, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2456, 3, 12, 42, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2465, 3, 13, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2474, 3, 14, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2483, 3, 15, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2492, 3, 16, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2501, 3, 17, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2510, 3, 18, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2519, 3, 19, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2528, 3, 20, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2537, 3, 21, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2546, 3, 22, 72, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2555, 0, 3, 30, 2429, 81, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2573, 0, 3, 33, 2438, 87, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2591, 0, 3, 36, 2447, 93, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2609, 0, 3, 39, 2456, 99, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2627, 0, 3, 42, 2465, 105, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2645, 0, 3, 45, 2474, 111, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2663, 0, 3, 48, 2483, 117, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2681, 0, 3, 51, 2492, 123, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2699, 0, 3, 54, 2501, 129, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2717, 0, 3, 57, 2510, 135, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2735, 0, 3, 60, 2519, 141, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2753, 0, 3, 63, 2528, 147, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2771, 0, 3, 66, 2537, 153, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2789, 0, 3, 69, 2546, 159, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2807, 0, 3, 81, 2573, 185, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2837, 0, 3, 87, 2591, 195, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2867, 0, 3, 93, 2609, 205, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2897, 0, 3, 99, 2627, 215, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2927, 0, 3, 105, 2645, 225, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2957, 0, 3, 111, 2663, 235, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2987, 0, 3, 117, 2681, 245, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3017, 0, 3, 123, 2699, 255, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3047, 0, 3, 129, 2717, 265, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3077, 0, 3, 135, 2735, 275, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3107, 0, 3, 141, 2753, 285, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3137, 0, 3, 147, 2771, 295, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3167, 0, 3, 153, 2789, 305, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3197, 0, 3, 185, 2837, 330, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3242, 0, 3, 195, 2867, 345, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3287, 0, 3, 205, 2897, 360, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3332, 0, 3, 215, 2927, 375, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3377, 0, 3, 225, 2957, 390, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3422, 0, 3, 235, 2987, 405, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3467, 0, 3, 245, 3017, 420, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3512, 0, 3, 255, 3047, 435, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3557, 0, 3, 265, 3077, 450, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3602, 0, 3, 275, 3107, 465, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3647, 0, 3, 285, 3137, 480, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3692, 0, 3, 295, 3167, 495, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3737, 0, 3, 330, 3242, 552, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3800, 0, 3, 345, 3287, 573, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3863, 0, 3, 360, 3332, 594, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3926, 0, 3, 375, 3377, 615, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3989, 0, 3, 390, 3422, 636, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4052, 0, 3, 405, 3467, 657, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4115, 0, 3, 420, 3512, 678, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4178, 0, 3, 435, 3557, 699, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4241, 0, 3, 450, 3602, 720, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4304, 0, 3, 465, 3647, 741, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4367, 0, 3, 480, 3692, 762, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4430, 0, 3, 552, 3800, 811, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4514, 0, 3, 573, 3863, 839, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4598, 0, 3, 594, 3926, 867, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4682, 0, 3, 615, 3989, 895, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4766, 0, 3, 636, 4052, 923, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4850, 0, 3, 657, 4115, 951, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4934, 0, 3, 678, 4178, 979, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5018, 0, 3, 699, 4241, 1007, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5102, 0, 3, 720, 4304, 1035, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5186, 0, 3, 741, 4367, 1063, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5270, 0, 3, 811, 4514, 1163, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5378, 0, 3, 839, 4598, 1199, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5486, 0, 3, 867, 4682, 1235, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5594, 0, 3, 895, 4766, 1271, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5702, 0, 3, 923, 4850, 1307, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5810, 0, 3, 951, 4934, 1343, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 5918, 0, 3, 979, 5018, 1379, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 6026, 0, 3, 1007, 5102, 1415, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6134, 0, 3, 1035, 5186, 1451, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6242, 0, 3, 1163, 5378, 1532, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6377, 0, 3, 1199, 5486, 1577, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6512, 0, 3, 1235, 5594, 1622, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6647, 0, 3, 1271, 5702, 1667, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6782, 0, 3, 1307, 5810, 1712, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 6917, 0, 3, 1343, 5918, 1757, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7052, 0, 3, 1379, 6026, 1802, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7187, 0, 3, 1415, 6134, 1847, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7322, 0, 3, 1532, 6377, 2002, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7487, 0, 3, 1577, 6512, 2057, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7652, 0, 3, 1622, 6647, 2112, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7817, 0, 3, 1667, 6782, 2167, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 7982, 0, 3, 1712, 6917, 2222, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 8147, 0, 3, 1757, 7052, 2277, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 8312, 0, 3, 1802, 7187, 2332, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 8477, 3, 9, 10, 2390, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8483, 3, 10, 11, 2393, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8489, 3, 11, 12, 2396, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8495, 3, 12, 13, 2399, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8501, 3, 13, 14, 2402, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8507, 3, 14, 15, 2405, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8513, 3, 15, 16, 2408, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8519, 3, 16, 17, 2411, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8525, 3, 17, 18, 2414, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8531, 3, 18, 19, 2417, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8537, 3, 19, 20, 2420, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8543, 3, 20, 21, 2423, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8549, 3, 21, 22, 2426, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 8555, 0, 3, 2387, 8477, 2438, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8573, 0, 3, 2390, 8483, 2447, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8591, 0, 3, 2393, 8489, 2456, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8609, 0, 3, 2396, 8495, 2465, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8627, 0, 3, 2399, 8501, 2474, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8645, 0, 3, 2402, 8507, 2483, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8663, 0, 3, 2405, 8513, 2492, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8681, 0, 3, 2408, 8519, 2501, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8699, 0, 3, 2411, 8525, 2510, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8717, 0, 3, 2414, 8531, 2519, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8735, 0, 3, 2417, 8537, 2528, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8753, 0, 3, 2420, 8543, 2537, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8771, 0, 3, 2423, 8549, 2546, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 8789, 0, 3, 2429, 8555, 75, 81, 2573,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8825, 0, 3, 2438, 8573, 81, 87, 2591,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8861, 0, 3, 2447, 8591, 87, 93, 2609,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8897, 0, 3, 2456, 8609, 93, 99, 2627,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8933, 0, 3, 2465, 8627, 99, 105, 2645,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 8969, 0, 3, 2474, 8645, 105, 111, 2663,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9005, 0, 3, 2483, 8663, 111, 117, 2681,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9041, 0, 3, 2492, 8681, 117, 123, 2699,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9077, 0, 3, 2501, 8699, 123, 129, 2717,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9113, 0, 3, 2510, 8717, 129, 135, 2735,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9149, 0, 3, 2519, 8735, 135, 141, 2753,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9185, 0, 3, 2528, 8753, 141, 147, 2771,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9221, 0, 3, 2537, 8771, 147, 153, 2789,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9257, 0, 3, 2555, 8789, 165, 175, 2807,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9317, 0, 3, 2573, 8825, 175, 185, 2837,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9377, 0, 3, 2591, 8861, 185, 195, 2867,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9437, 0, 3, 2609, 8897, 195, 205, 2897,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9497, 0, 3, 2627, 8933, 205, 215, 2927,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9557, 0, 3, 2645, 8969, 215, 225, 2957,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9617, 0, 3, 2663, 9005, 225, 235, 2987,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9677, 0, 3, 2681, 9041, 235, 245, 3017,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9737, 0, 3, 2699, 9077, 245, 255, 3047,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9797, 0, 3, 2717, 9113, 255, 265, 3077,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9857, 0, 3, 2735, 9149, 265, 275, 3107,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9917, 0, 3, 2753, 9185, 275, 285, 3137,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9977, 0, 3, 2771, 9221, 285, 295, 3167,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10037, 0, 3, 8789, 8825, 2837, 9377,
                                                 315, 330, 3242, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10127, 0, 3, 8825, 8861, 2867, 9437,
                                                 330, 345, 3287, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10217, 0, 3, 8861, 8897, 2897, 9497,
                                                 345, 360, 3332, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10307, 0, 3, 8897, 8933, 2927, 9557,
                                                 360, 375, 3377, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10397, 0, 3, 8933, 8969, 2957, 9617,
                                                 375, 390, 3422, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10487, 0, 3, 8969, 9005, 2987, 9677,
                                                 390, 405, 3467, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10577, 0, 3, 9005, 9041, 3017, 9737,
                                                 405, 420, 3512, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10667, 0, 3, 9041, 9077, 3047, 9797,
                                                 420, 435, 3557, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10757, 0, 3, 9077, 9113, 3077, 9857,
                                                 435, 450, 3602, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10847, 0, 3, 9113, 9149, 3107, 9917,
                                                 450, 465, 3647, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10937, 0, 3, 9149, 9185, 3137, 9977,
                                                 465, 480, 3692, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11027, 0, 3, 9257, 9317, 3197, 10037,
                                                 510, 531, 3737, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11153, 0, 3, 9317, 9377, 3242, 10127,
                                                 531, 552, 3800, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11279, 0, 3, 9377, 9437, 3287, 10217,
                                                 552, 573, 3863, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11405, 0, 3, 9437, 9497, 3332, 10307,
                                                 573, 594, 3926, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11531, 0, 3, 9497, 9557, 3377, 10397,
                                                 594, 615, 3989, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11657, 0, 3, 9557, 9617, 3422, 10487,
                                                 615, 636, 4052, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11783, 0, 3, 9617, 9677, 3467, 10577,
                                                 636, 657, 4115, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11909, 0, 3, 9677, 9737, 3512, 10667,
                                                 657, 678, 4178, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12035, 0, 3, 9737, 9797, 3557, 10757,
                                                 678, 699, 4241, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12161, 0, 3, 9797, 9857, 3602, 10847,
                                                 699, 720, 4304, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12287, 0, 3, 9857, 9917, 3647, 10937,
                                                 720, 741, 4367, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12413, 0, 3, 10037, 10127, 3800, 11279,
                                                 783, 811, 4514, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12581, 0, 3, 10127, 10217, 3863, 11405,
                                                 811, 839, 4598, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12749, 0, 3, 10217, 10307, 3926, 11531,
                                                 839, 867, 4682, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12917, 0, 3, 10307, 10397, 3989, 11657,
                                                 867, 895, 4766, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13085, 0, 3, 10397, 10487, 4052, 11783,
                                                 895, 923, 4850, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13253, 0, 3, 10487, 10577, 4115, 11909,
                                                 923, 951, 4934, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13421, 0, 3, 10577, 10667, 4178, 12035,
                                                 951, 979, 5018, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13589, 0, 3, 10667, 10757, 4241, 12161,
                                                 979, 1007, 5102, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13757, 0, 3, 10757, 10847, 4304, 12287,
                                                 1007, 1035, 5186, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13925, 0, 3, 11027, 11153, 4430, 12413,
                                                 1091, 1127, 5270, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14141, 0, 3, 11153, 11279, 4514, 12581,
                                                 1127, 1163, 5378, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14357, 0, 3, 11279, 11405, 4598, 12749,
                                                 1163, 1199, 5486, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14573, 0, 3, 11405, 11531, 4682, 12917,
                                                 1199, 1235, 5594, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14789, 0, 3, 11531, 11657, 4766, 13085,
                                                 1235, 1271, 5702, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15005, 0, 3, 11657, 11783, 4850, 13253,
                                                 1271, 1307, 5810, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15221, 0, 3, 11783, 11909, 4934, 13421,
                                                 1307, 1343, 5918, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15437, 0, 3, 11909, 12035, 5018, 13589,
                                                 1343, 1379, 6026, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15653, 0, 3, 12035, 12161, 5102, 13757,
                                                 1379, 1415, 6134, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 15869, 0, 3, 12413, 12581, 5378, 14357,
                                                 1487, 1532, 6377, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 16139, 0, 3, 12581, 12749, 5486, 14573,
                                                 1532, 1577, 6512, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 16409, 0, 3, 12749, 12917, 5594, 14789,
                                                 1577, 1622, 6647, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 16679, 0, 3, 12917, 13085, 5702, 15005,
                                                 1622, 1667, 6782, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 16949, 0, 3, 13085, 13253, 5810, 15221,
                                                 1667, 1712, 6917, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 17219, 0, 3, 13253, 13421, 5918, 15437,
                                                 1712, 1757, 7052, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 17489, 0, 3, 13421, 13589, 6026, 15653,
                                                 1757, 1802, 7187, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 17759, 0, 3, 13925, 14141, 6242, 15869,
                                                 1892, 1947, 7322, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 18089, 0, 3, 14141, 14357, 6377, 16139,
                                                 1947, 2002, 7487, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 18419, 0, 3, 14357, 14573, 6512, 16409,
                                                 2002, 2057, 7652, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 18749, 0, 3, 14573, 14789, 6647, 16679,
                                                 2057, 2112, 7817, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 19079, 0, 3, 14789, 15005, 6782, 16949,
                                                 2112, 2167, 7982, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 19409, 0, 3, 15005, 15221, 6917, 17219,
                                                 2167, 2222, 8147, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 19739, 0, 3, 15221, 15437, 7052, 17489,
                                                 2222, 2277, 8312, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20069, 3, 2387, 2390, 8483, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20079, 3, 2390, 2393, 8489, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20089, 3, 2393, 2396, 8495, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20099, 3, 2396, 2399, 8501, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20109, 3, 2399, 2402, 8507, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20119, 3, 2402, 2405, 8513, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20129, 3, 2405, 2408, 8519, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20139, 3, 2408, 2411, 8525, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20149, 3, 2411, 2414, 8531, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20159, 3, 2414, 2417, 8537, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20169, 3, 2417, 2420, 8543, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20179, 3, 2420, 2423, 8549, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 20189, 0, 3, 8477, 20069, 8573, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20219, 0, 3, 8483, 20079, 8591, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20249, 0, 3, 8489, 20089, 8609, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20279, 0, 3, 8495, 20099, 8627, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20309, 0, 3, 8501, 20109, 8645, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20339, 0, 3, 8507, 20119, 8663, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20369, 0, 3, 8513, 20129, 8681, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20399, 0, 3, 8519, 20139, 8699, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20429, 0, 3, 8525, 20149, 8717, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20459, 0, 3, 8531, 20159, 8735, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20489, 0, 3, 8537, 20169, 8753, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20519, 0, 3, 8543, 20179, 8771, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 20549, 0, 3, 8555, 20189, 2555, 2573,
                                                 8825, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20609, 0, 3, 8573, 20219, 2573, 2591,
                                                 8861, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20669, 0, 3, 8591, 20249, 2591, 2609,
                                                 8897, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20729, 0, 3, 8609, 20279, 2609, 2627,
                                                 8933, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20789, 0, 3, 8627, 20309, 2627, 2645,
                                                 8969, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20849, 0, 3, 8645, 20339, 2645, 2663,
                                                 9005, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20909, 0, 3, 8663, 20369, 2663, 2681,
                                                 9041, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20969, 0, 3, 8681, 20399, 2681, 2699,
                                                 9077, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21029, 0, 3, 8699, 20429, 2699, 2717,
                                                 9113, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21089, 0, 3, 8717, 20459, 2717, 2735,
                                                 9149, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21149, 0, 3, 8735, 20489, 2735, 2753,
                                                 9185, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21209, 0, 3, 8753, 20519, 2753, 2771,
                                                 9221, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21269, 0, 3, 8825, 20609, 2807, 2837,
                                                 9377, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21369, 0, 3, 8861, 20669, 2837, 2867,
                                                 9437, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21469, 0, 3, 8897, 20729, 2867, 2897,
                                                 9497, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21569, 0, 3, 8933, 20789, 2897, 2927,
                                                 9557, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21669, 0, 3, 8969, 20849, 2927, 2957,
                                                 9617, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21769, 0, 3, 9005, 20909, 2957, 2987,
                                                 9677, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21869, 0, 3, 9041, 20969, 2987, 3017,
                                                 9737, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21969, 0, 3, 9077, 21029, 3017, 3047,
                                                 9797, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 22069, 0, 3, 9113, 21089, 3047, 3077,
                                                 9857, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 22169, 0, 3, 9149, 21149, 3077, 3107,
                                                 9917, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 22269, 0, 3, 9185, 21209, 3107, 3137,
                                                 9977, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22369, 0, 3, 20549, 20609, 9377, 21369,
                                                 3197, 3242, 10127, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22519, 0, 3, 20609, 20669, 9437, 21469,
                                                 3242, 3287, 10217, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22669, 0, 3, 20669, 20729, 9497, 21569,
                                                 3287, 3332, 10307, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22819, 0, 3, 20729, 20789, 9557, 21669,
                                                 3332, 3377, 10397, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22969, 0, 3, 20789, 20849, 9617, 21769,
                                                 3377, 3422, 10487, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23119, 0, 3, 20849, 20909, 9677, 21869,
                                                 3422, 3467, 10577, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23269, 0, 3, 20909, 20969, 9737, 21969,
                                                 3467, 3512, 10667, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23419, 0, 3, 20969, 21029, 9797, 22069,
                                                 3512, 3557, 10757, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23569, 0, 3, 21029, 21089, 9857, 22169,
                                                 3557, 3602, 10847, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23719, 0, 3, 21089, 21149, 9917, 22269,
                                                 3602, 3647, 10937, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 23869, 0, 3, 21269, 21369, 10127, 22519,
                                                 3737, 3800, 11279, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24079, 0, 3, 21369, 21469, 10217, 22669,
                                                 3800, 3863, 11405, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24289, 0, 3, 21469, 21569, 10307, 22819,
                                                 3863, 3926, 11531, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24499, 0, 3, 21569, 21669, 10397, 22969,
                                                 3926, 3989, 11657, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24709, 0, 3, 21669, 21769, 10487, 23119,
                                                 3989, 4052, 11783, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24919, 0, 3, 21769, 21869, 10577, 23269,
                                                 4052, 4115, 11909, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25129, 0, 3, 21869, 21969, 10667, 23419,
                                                 4115, 4178, 12035, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25339, 0, 3, 21969, 22069, 10757, 23569,
                                                 4178, 4241, 12161, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25549, 0, 3, 22069, 22169, 10847, 23719,
                                                 4241, 4304, 12287, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 25759, 0, 3, 22369, 22519, 11279, 24079,
                                                 4430, 4514, 12581, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26039, 0, 3, 22519, 22669, 11405, 24289,
                                                 4514, 4598, 12749, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26319, 0, 3, 22669, 22819, 11531, 24499,
                                                 4598, 4682, 12917, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26599, 0, 3, 22819, 22969, 11657, 24709,
                                                 4682, 4766, 13085, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26879, 0, 3, 22969, 23119, 11783, 24919,
                                                 4766, 4850, 13253, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27159, 0, 3, 23119, 23269, 11909, 25129,
                                                 4850, 4934, 13421, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27439, 0, 3, 23269, 23419, 12035, 25339,
                                                 4934, 5018, 13589, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27719, 0, 3, 23419, 23569, 12161, 25549,
                                                 5018, 5102, 13757, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 27999, 0, 3, 23869, 24079, 12581, 26039,
                                                 5270, 5378, 14357, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 28359, 0, 3, 24079, 24289, 12749, 26319,
                                                 5378, 5486, 14573, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 28719, 0, 3, 24289, 24499, 12917, 26599,
                                                 5486, 5594, 14789, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29079, 0, 3, 24499, 24709, 13085, 26879,
                                                 5594, 5702, 15005, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29439, 0, 3, 24709, 24919, 13253, 27159,
                                                 5702, 5810, 15221, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29799, 0, 3, 24919, 25129, 13421, 27439,
                                                 5810, 5918, 15437, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 30159, 0, 3, 25129, 25339, 13589, 27719,
                                                 5918, 6026, 15653, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 30519, 0, 3, 25759, 26039, 14357, 28359,
                                                 6242, 6377, 16139, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 30969, 0, 3, 26039, 26319, 14573, 28719,
                                                 6377, 6512, 16409, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 31419, 0, 3, 26319, 26599, 14789, 29079,
                                                 6512, 6647, 16679, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 31869, 0, 3, 26599, 26879, 15005, 29439,
                                                 6647, 6782, 16949, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 32319, 0, 3, 26879, 27159, 15221, 29799,
                                                 6782, 6917, 17219, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 32769, 0, 3, 27159, 27439, 15437, 30159,
                                                 6917, 7052, 17489, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 33219, 0, 3, 27999, 28359, 16139, 30969,
                                                 7322, 7487, 18419, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 33769, 0, 3, 28359, 28719, 16409, 31419,
                                                 7487, 7652, 18749, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 34319, 0, 3, 28719, 29079, 16679, 31869,
                                                 7652, 7817, 19079, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 34869, 0, 3, 29079, 29439, 16949, 32319,
                                                 7817, 7982, 19409, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 35419, 0, 3, 29439, 29799, 17219, 32769,
                                                 7982, 8147, 19739, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 35969, 3, 8477, 8483, 20079, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 35984, 3, 8483, 8489, 20089, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 35999, 3, 8489, 8495, 20099, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36014, 3, 8495, 8501, 20109, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36029, 3, 8501, 8507, 20119, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36044, 3, 8507, 8513, 20129, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36059, 3, 8513, 8519, 20139, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36074, 3, 8519, 8525, 20149, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36089, 3, 8525, 8531, 20159, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36104, 3, 8531, 8537, 20169, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36119, 3, 8537, 8543, 20179, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 36134, 0, 3, 20069, 35969, 20219, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36179, 0, 3, 20079, 35984, 20249, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36224, 0, 3, 20089, 35999, 20279, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36269, 0, 3, 20099, 36014, 20309, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36314, 0, 3, 20109, 36029, 20339, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36359, 0, 3, 20119, 36044, 20369, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36404, 0, 3, 20129, 36059, 20399, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36449, 0, 3, 20139, 36074, 20429, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36494, 0, 3, 20149, 36089, 20459, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36539, 0, 3, 20159, 36104, 20489, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36584, 0, 3, 20169, 36119, 20519, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 36629, 0, 3, 20189, 36134, 8789, 8825,
                                                 20609, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 36719, 0, 3, 20219, 36179, 8825, 8861,
                                                 20669, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 36809, 0, 3, 20249, 36224, 8861, 8897,
                                                 20729, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 36899, 0, 3, 20279, 36269, 8897, 8933,
                                                 20789, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 36989, 0, 3, 20309, 36314, 8933, 8969,
                                                 20849, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37079, 0, 3, 20339, 36359, 8969, 9005,
                                                 20909, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37169, 0, 3, 20369, 36404, 9005, 9041,
                                                 20969, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37259, 0, 3, 20399, 36449, 9041, 9077,
                                                 21029, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37349, 0, 3, 20429, 36494, 9077, 9113,
                                                 21089, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37439, 0, 3, 20459, 36539, 9113, 9149,
                                                 21149, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37529, 0, 3, 20489, 36584, 9149, 9185,
                                                 21209, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37619, 0, 3, 20549, 36629, 9257, 9317,
                                                 21269, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37769, 0, 3, 20609, 36719, 9317, 9377,
                                                 21369, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37919, 0, 3, 20669, 36809, 9377, 9437,
                                                 21469, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38069, 0, 3, 20729, 36899, 9437, 9497,
                                                 21569, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38219, 0, 3, 20789, 36989, 9497, 9557,
                                                 21669, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38369, 0, 3, 20849, 37079, 9557, 9617,
                                                 21769, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38519, 0, 3, 20909, 37169, 9617, 9677,
                                                 21869, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38669, 0, 3, 20969, 37259, 9677, 9737,
                                                 21969, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38819, 0, 3, 21029, 37349, 9737, 9797,
                                                 22069, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38969, 0, 3, 21089, 37439, 9797, 9857,
                                                 22169, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39119, 0, 3, 21149, 37529, 9857, 9917,
                                                 22269, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39269, 0, 3, 36629, 36719, 21369, 37919,
                                                 10037, 10127, 22519, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39494, 0, 3, 36719, 36809, 21469, 38069,
                                                 10127, 10217, 22669, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39719, 0, 3, 36809, 36899, 21569, 38219,
                                                 10217, 10307, 22819, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39944, 0, 3, 36899, 36989, 21669, 38369,
                                                 10307, 10397, 22969, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 40169, 0, 3, 36989, 37079, 21769, 38519,
                                                 10397, 10487, 23119, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 40394, 0, 3, 37079, 37169, 21869, 38669,
                                                 10487, 10577, 23269, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 40619, 0, 3, 37169, 37259, 21969, 38819,
                                                 10577, 10667, 23419, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 40844, 0, 3, 37259, 37349, 22069, 38969,
                                                 10667, 10757, 23569, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 41069, 0, 3, 37349, 37439, 22169, 39119,
                                                 10757, 10847, 23719, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41294, 0, 3, 37619, 37769, 22369, 39269,
                                                 11027, 11153, 23869, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41609, 0, 3, 37769, 37919, 22519, 39494,
                                                 11153, 11279, 24079, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41924, 0, 3, 37919, 38069, 22669, 39719,
                                                 11279, 11405, 24289, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 42239, 0, 3, 38069, 38219, 22819, 39944,
                                                 11405, 11531, 24499, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 42554, 0, 3, 38219, 38369, 22969, 40169,
                                                 11531, 11657, 24709, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 42869, 0, 3, 38369, 38519, 23119, 40394,
                                                 11657, 11783, 24919, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 43184, 0, 3, 38519, 38669, 23269, 40619,
                                                 11783, 11909, 25129, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 43499, 0, 3, 38669, 38819, 23419, 40844,
                                                 11909, 12035, 25339, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 43814, 0, 3, 38819, 38969, 23569, 41069,
                                                 12035, 12161, 25549, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 44129, 0, 3, 39269, 39494, 24079, 41924,
                                                 12413, 12581, 26039, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 44549, 0, 3, 39494, 39719, 24289, 42239,
                                                 12581, 12749, 26319, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 44969, 0, 3, 39719, 39944, 24499, 42554,
                                                 12749, 12917, 26599, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 45389, 0, 3, 39944, 40169, 24709, 42869,
                                                 12917, 13085, 26879, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 45809, 0, 3, 40169, 40394, 24919, 43184,
                                                 13085, 13253, 27159, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 46229, 0, 3, 40394, 40619, 25129, 43499,
                                                 13253, 13421, 27439, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 46649, 0, 3, 40619, 40844, 25339, 43814,
                                                 13421, 13589, 27719, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 47069, 0, 3, 41294, 41609, 25759, 44129,
                                                 13925, 14141, 27999, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 47609, 0, 3, 41609, 41924, 26039, 44549,
                                                 14141, 14357, 28359, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 48149, 0, 3, 41924, 42239, 26319, 44969,
                                                 14357, 14573, 28719, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 48689, 0, 3, 42239, 42554, 26599, 45389,
                                                 14573, 14789, 29079, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 49229, 0, 3, 42554, 42869, 26879, 45809,
                                                 14789, 15005, 29439, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 49769, 0, 3, 42869, 43184, 27159, 46229,
                                                 15005, 15221, 29799, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 50309, 0, 3, 43184, 43499, 27439, 46649,
                                                 15221, 15437, 30159, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 50849, 0, 3, 44129, 44549, 28359, 48149,
                                                 15869, 16139, 30969, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 51524, 0, 3, 44549, 44969, 28719, 48689,
                                                 16139, 16409, 31419, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 52199, 0, 3, 44969, 45389, 29079, 49229,
                                                 16409, 16679, 31869, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 52874, 0, 3, 45389, 45809, 29439, 49769,
                                                 16679, 16949, 32319, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 53549, 0, 3, 45809, 46229, 29799, 50309,
                                                 16949, 17219, 32769, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 54224, 0, 3, 47069, 47609, 30519, 50849,
                                                 17759, 18089, 33219, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 55049, 0, 3, 47609, 48149, 30969, 51524,
                                                 18089, 18419, 33769, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 55874, 0, 3, 48149, 48689, 31419, 52199,
                                                 18419, 18749, 34319, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 56699, 0, 3, 48689, 49229, 31869, 52874,
                                                 18749, 19079, 34869, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 57524, 0, 3, 49229, 49769, 32319, 53549,
                                                 19079, 19409, 35419, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58349, 3, 20069, 20079, 35984, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58370, 3, 20079, 20089, 35999, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58391, 3, 20089, 20099, 36014, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58412, 3, 20099, 20109, 36029, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58433, 3, 20109, 20119, 36044, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58454, 3, 20119, 20129, 36059, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58475, 3, 20129, 20139, 36074, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58496, 3, 20139, 20149, 36089, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58517, 3, 20149, 20159, 36104, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 58538, 3, 20159, 20169, 36119, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 58559, 0, 3, 35969, 58349, 36179, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 58622, 0, 3, 35984, 58370, 36224, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 58685, 0, 3, 35999, 58391, 36269, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 58748, 0, 3, 36014, 58412, 36314, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 58811, 0, 3, 36029, 58433, 36359, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 58874, 0, 3, 36044, 58454, 36404, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 58937, 0, 3, 36059, 58475, 36449, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 59000, 0, 3, 36074, 58496, 36494, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 59063, 0, 3, 36089, 58517, 36539, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 59126, 0, 3, 36104, 58538, 36584, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 59189, 0, 3, 36134, 58559, 20549, 20609,
                                                 36719, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 59315, 0, 3, 36179, 58622, 20609, 20669,
                                                 36809, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 59441, 0, 3, 36224, 58685, 20669, 20729,
                                                 36899, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 59567, 0, 3, 36269, 58748, 20729, 20789,
                                                 36989, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 59693, 0, 3, 36314, 58811, 20789, 20849,
                                                 37079, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 59819, 0, 3, 36359, 58874, 20849, 20909,
                                                 37169, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 59945, 0, 3, 36404, 58937, 20909, 20969,
                                                 37259, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 60071, 0, 3, 36449, 59000, 20969, 21029,
                                                 37349, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 60197, 0, 3, 36494, 59063, 21029, 21089,
                                                 37439, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 60323, 0, 3, 36539, 59126, 21089, 21149,
                                                 37529, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 60449, 0, 3, 36719, 59315, 21269, 21369,
                                                 37919, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 60659, 0, 3, 36809, 59441, 21369, 21469,
                                                 38069, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 60869, 0, 3, 36899, 59567, 21469, 21569,
                                                 38219, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 61079, 0, 3, 36989, 59693, 21569, 21669,
                                                 38369, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 61289, 0, 3, 37079, 59819, 21669, 21769,
                                                 38519, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 61499, 0, 3, 37169, 59945, 21769, 21869,
                                                 38669, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 61709, 0, 3, 37259, 60071, 21869, 21969,
                                                 38819, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 61919, 0, 3, 37349, 60197, 21969, 22069,
                                                 38969, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 62129, 0, 3, 37439, 60323, 22069, 22169,
                                                 39119, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 62339, 0, 3, 59189, 59315, 37919, 60659,
                                                 22369, 22519, 39494, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 62654, 0, 3, 59315, 59441, 38069, 60869,
                                                 22519, 22669, 39719, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 62969, 0, 3, 59441, 59567, 38219, 61079,
                                                 22669, 22819, 39944, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 63284, 0, 3, 59567, 59693, 38369, 61289,
                                                 22819, 22969, 40169, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 63599, 0, 3, 59693, 59819, 38519, 61499,
                                                 22969, 23119, 40394, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 63914, 0, 3, 59819, 59945, 38669, 61709,
                                                 23119, 23269, 40619, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 64229, 0, 3, 59945, 60071, 38819, 61919,
                                                 23269, 23419, 40844, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 64544, 0, 3, 60071, 60197, 38969, 62129,
                                                 23419, 23569, 41069, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 64859, 0, 3, 60449, 60659, 39494, 62654,
                                                 23869, 24079, 41924, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 65300, 0, 3, 60659, 60869, 39719, 62969,
                                                 24079, 24289, 42239, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 65741, 0, 3, 60869, 61079, 39944, 63284,
                                                 24289, 24499, 42554, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 66182, 0, 3, 61079, 61289, 40169, 63599,
                                                 24499, 24709, 42869, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 66623, 0, 3, 61289, 61499, 40394, 63914,
                                                 24709, 24919, 43184, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 67064, 0, 3, 61499, 61709, 40619, 64229,
                                                 24919, 25129, 43499, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 67505, 0, 3, 61709, 61919, 40844, 64544,
                                                 25129, 25339, 43814, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 67946, 0, 3, 62339, 62654, 41924, 65300,
                                                 25759, 26039, 44549, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 68534, 0, 3, 62654, 62969, 42239, 65741,
                                                 26039, 26319, 44969, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 69122, 0, 3, 62969, 63284, 42554, 66182,
                                                 26319, 26599, 45389, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 69710, 0, 3, 63284, 63599, 42869, 66623,
                                                 26599, 26879, 45809, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 70298, 0, 3, 63599, 63914, 43184, 67064,
                                                 26879, 27159, 46229, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 70886, 0, 3, 63914, 64229, 43499, 67505,
                                                 27159, 27439, 46649, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 71474, 0, 3, 64859, 65300, 44549, 68534,
                                                 27999, 28359, 48149, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 72230, 0, 3, 65300, 65741, 44969, 69122,
                                                 28359, 28719, 48689, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 72986, 0, 3, 65741, 66182, 45389, 69710,
                                                 28719, 29079, 49229, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 73742, 0, 3, 66182, 66623, 45809, 70298,
                                                 29079, 29439, 49769, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 74498, 0, 3, 66623, 67064, 46229, 70886,
                                                 29439, 29799, 50309, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 75254, 0, 3, 67946, 68534, 48149, 72230,
                                                 30519, 30969, 51524, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 76199, 0, 3, 68534, 69122, 48689, 72986,
                                                 30969, 31419, 52199, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 77144, 0, 3, 69122, 69710, 49229, 73742,
                                                 31419, 31869, 52874, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 78089, 0, 3, 69710, 70298, 49769, 74498,
                                                 31869, 32319, 53549, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 79034, 0, 3, 71474, 72230, 51524, 76199,
                                                 33219, 33769, 55874, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 80189, 0, 3, 72230, 72986, 52199, 77144,
                                                 33769, 34319, 56699, ncols, alpha, beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 81344, 0, 3, 72986, 73742, 52874, 78089,
                                                 34319, 34869, 57524, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82499, 3, 35969, 35984, 58370, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82527, 3, 35984, 35999, 58391, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82555, 3, 35999, 36014, 58412, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82583, 3, 36014, 36029, 58433, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82611, 3, 36029, 36044, 58454, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82639, 3, 36044, 36059, 58475, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82667, 3, 36059, 36074, 58496, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82695, 3, 36074, 36089, 58517, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82723, 3, 36089, 36104, 58538, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 82751, 0, 3, 58349, 82499, 58622, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 82835, 0, 3, 58370, 82527, 58685, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 82919, 0, 3, 58391, 82555, 58748, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83003, 0, 3, 58412, 82583, 58811, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83087, 0, 3, 58433, 82611, 58874, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83171, 0, 3, 58454, 82639, 58937, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83255, 0, 3, 58475, 82667, 59000, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83339, 0, 3, 58496, 82695, 59063, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83423, 0, 3, 58517, 82723, 59126, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 83507, 0, 3, 58559, 82751, 36629, 36719,
                                                 59315, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 83675, 0, 3, 58622, 82835, 36719, 36809,
                                                 59441, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 83843, 0, 3, 58685, 82919, 36809, 36899,
                                                 59567, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84011, 0, 3, 58748, 83003, 36899, 36989,
                                                 59693, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84179, 0, 3, 58811, 83087, 36989, 37079,
                                                 59819, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84347, 0, 3, 58874, 83171, 37079, 37169,
                                                 59945, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84515, 0, 3, 58937, 83255, 37169, 37259,
                                                 60071, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84683, 0, 3, 59000, 83339, 37259, 37349,
                                                 60197, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84851, 0, 3, 59063, 83423, 37349, 37439,
                                                 60323, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 85019, 0, 3, 59189, 83507, 37619, 37769,
                                                 60449, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 85299, 0, 3, 59315, 83675, 37769, 37919,
                                                 60659, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 85579, 0, 3, 59441, 83843, 37919, 38069,
                                                 60869, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 85859, 0, 3, 59567, 84011, 38069, 38219,
                                                 61079, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 86139, 0, 3, 59693, 84179, 38219, 38369,
                                                 61289, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 86419, 0, 3, 59819, 84347, 38369, 38519,
                                                 61499, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 86699, 0, 3, 59945, 84515, 38519, 38669,
                                                 61709, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 86979, 0, 3, 60071, 84683, 38669, 38819,
                                                 61919, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 87259, 0, 3, 60197, 84851, 38819, 38969,
                                                 62129, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 87539, 0, 3, 83507, 83675, 60659, 85579,
                                                 39269, 39494, 62654, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 87959, 0, 3, 83675, 83843, 60869, 85859,
                                                 39494, 39719, 62969, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 88379, 0, 3, 83843, 84011, 61079, 86139,
                                                 39719, 39944, 63284, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 88799, 0, 3, 84011, 84179, 61289, 86419,
                                                 39944, 40169, 63599, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 89219, 0, 3, 84179, 84347, 61499, 86699,
                                                 40169, 40394, 63914, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 89639, 0, 3, 84347, 84515, 61709, 86979,
                                                 40394, 40619, 64229, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 90059, 0, 3, 84515, 84683, 61919, 87259,
                                                 40619, 40844, 64544, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 90479, 0, 3, 85019, 85299, 62339, 87539,
                                                 41294, 41609, 64859, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 91067, 0, 3, 85299, 85579, 62654, 87959,
                                                 41609, 41924, 65300, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 91655, 0, 3, 85579, 85859, 62969, 88379,
                                                 41924, 42239, 65741, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 92243, 0, 3, 85859, 86139, 63284, 88799,
                                                 42239, 42554, 66182, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 92831, 0, 3, 86139, 86419, 63599, 89219,
                                                 42554, 42869, 66623, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 93419, 0, 3, 86419, 86699, 63914, 89639,
                                                 42869, 43184, 67064, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 94007, 0, 3, 86699, 86979, 64229, 90059,
                                                 43184, 43499, 67505, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 94595, 0, 3, 87539, 87959, 65300, 91655,
                                                 44129, 44549, 68534, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 95379, 0, 3, 87959, 88379, 65741, 92243,
                                                 44549, 44969, 69122, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 96163, 0, 3, 88379, 88799, 66182, 92831,
                                                 44969, 45389, 69710, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 96947, 0, 3, 88799, 89219, 66623, 93419,
                                                 45389, 45809, 70298, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 97731, 0, 3, 89219, 89639, 67064, 94007,
                                                 45809, 46229, 70886, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 98515, 0, 3, 90479, 91067, 67946, 94595,
                                                 47069, 47609, 71474, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 99523, 0, 3, 91067, 91655, 68534, 95379,
                                                 47609, 48149, 72230, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 100531, 0, 3, 91655, 92243, 69122,
                                                 96163, 48149, 48689, 72986, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 101539, 0, 3, 92243, 92831, 69710,
                                                 96947, 48689, 49229, 73742, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 102547, 0, 3, 92831, 93419, 70298,
                                                 97731, 49229, 49769, 74498, ncols, alpha, beta,
                                                 p);

            compute_prim_li_electron_repulsion_0(buffer, 103555, 0, 3, 94595, 95379, 72230,
                                                 100531, 50849, 51524, 76199, ncols, alpha, beta,
                                                 p);

            compute_prim_li_electron_repulsion_0(buffer, 104815, 0, 3, 95379, 96163, 72986,
                                                 101539, 51524, 52199, 77144, ncols, alpha, beta,
                                                 p);

            compute_prim_li_electron_repulsion_0(buffer, 106075, 0, 3, 96163, 96947, 73742,
                                                 102547, 52199, 52874, 78089, ncols, alpha, beta,
                                                 p);

            compute_prim_mi_electron_repulsion_0(buffer, 107335, 0, 3, 98515, 99523, 75254,
                                                 103555, 54224, 55049, 79034, ncols, alpha, beta,
                                                 p);

            compute_prim_mi_electron_repulsion_0(buffer, 108875, 0, 3, 99523, 100531, 76199,
                                                 104815, 55049, 55874, 80189, ncols, alpha, beta,
                                                 p);

            compute_prim_mi_electron_repulsion_0(buffer, 110415, 0, 3, 100531, 101539, 77144,
                                                 106075, 55874, 56699, 81344, ncols, alpha, beta,
                                                 p);

            compute_prim_sk_electron_repulsion_0(buffer, 111955, 3, 58349, 58370, 82527, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 111991, 3, 58370, 58391, 82555, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112027, 3, 58391, 58412, 82583, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112063, 3, 58412, 58433, 82611, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112099, 3, 58433, 58454, 82639, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112135, 3, 58454, 58475, 82667, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112171, 3, 58475, 58496, 82695, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112207, 3, 58496, 58517, 82723, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112243, 0, 3, 82499, 111955, 82835,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112351, 0, 3, 82527, 111991, 82919,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112459, 0, 3, 82555, 112027, 83003,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112567, 0, 3, 82583, 112063, 83087,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112675, 0, 3, 82611, 112099, 83171,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112783, 0, 3, 82639, 112135, 83255,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112891, 0, 3, 82667, 112171, 83339,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112999, 0, 3, 82695, 112207, 83423,
                                                 ncols, p);

            compute_prim_dk_electron_repulsion_0(buffer, 113107, 0, 3, 82751, 112243, 59189,
                                                 59315, 83675, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 113323, 0, 3, 82835, 112351, 59315,
                                                 59441, 83843, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 113539, 0, 3, 82919, 112459, 59441,
                                                 59567, 84011, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 113755, 0, 3, 83003, 112567, 59567,
                                                 59693, 84179, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 113971, 0, 3, 83087, 112675, 59693,
                                                 59819, 84347, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 114187, 0, 3, 83171, 112783, 59819,
                                                 59945, 84515, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 114403, 0, 3, 83255, 112891, 59945,
                                                 60071, 84683, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 114619, 0, 3, 83339, 112999, 60071,
                                                 60197, 84851, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 114835, 0, 3, 83675, 113323, 60449,
                                                 60659, 85579, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 115195, 0, 3, 83843, 113539, 60659,
                                                 60869, 85859, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 115555, 0, 3, 84011, 113755, 60869,
                                                 61079, 86139, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 115915, 0, 3, 84179, 113971, 61079,
                                                 61289, 86419, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 116275, 0, 3, 84347, 114187, 61289,
                                                 61499, 86699, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 116635, 0, 3, 84515, 114403, 61499,
                                                 61709, 86979, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 116995, 0, 3, 84683, 114619, 61709,
                                                 61919, 87259, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 117355, 0, 3, 113107, 113323, 85579,
                                                 115195, 62339, 62654, 87959, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 117895, 0, 3, 113323, 113539, 85859,
                                                 115555, 62654, 62969, 88379, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 118435, 0, 3, 113539, 113755, 86139,
                                                 115915, 62969, 63284, 88799, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 118975, 0, 3, 113755, 113971, 86419,
                                                 116275, 63284, 63599, 89219, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 119515, 0, 3, 113971, 114187, 86699,
                                                 116635, 63599, 63914, 89639, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 120055, 0, 3, 114187, 114403, 86979,
                                                 116995, 63914, 64229, 90059, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 120595, 0, 3, 114835, 115195, 87959,
                                                 117895, 64859, 65300, 91655, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 121351, 0, 3, 115195, 115555, 88379,
                                                 118435, 65300, 65741, 92243, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 122107, 0, 3, 115555, 115915, 88799,
                                                 118975, 65741, 66182, 92831, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 122863, 0, 3, 115915, 116275, 89219,
                                                 119515, 66182, 66623, 93419, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 123619, 0, 3, 116275, 116635, 89639,
                                                 120055, 66623, 67064, 94007, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 124375, 0, 3, 117355, 117895, 91655,
                                                 121351, 67946, 68534, 95379, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 125383, 0, 3, 117895, 118435, 92243,
                                                 122107, 68534, 69122, 96163, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 126391, 0, 3, 118435, 118975, 92831,
                                                 122863, 69122, 69710, 96947, ncols, alpha, beta,
                                                 p);

            compute_prim_ik_electron_repulsion_0(buffer, 127399, 0, 3, 118975, 119515, 93419,
                                                 123619, 69710, 70298, 97731, ncols, alpha, beta,
                                                 p);

            compute_prim_kk_electron_repulsion_0(buffer, 128407, 0, 3, 120595, 121351, 95379,
                                                 125383, 71474, 72230, 100531, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 129703, 0, 3, 121351, 122107, 96163,
                                                 126391, 72230, 72986, 101539, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 130999, 0, 3, 122107, 122863, 96947,
                                                 127399, 72986, 73742, 102547, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 132295, 0, 3, 124375, 125383, 100531,
                                                 129703, 75254, 76199, 104815, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 133915, 0, 3, 125383, 126391, 101539,
                                                 130999, 76199, 77144, 106075, ncols, alpha,
                                                 beta, p);

            compute_prim_mk_electron_repulsion_0(buffer, 135535, 0, 3, 128407, 129703, 104815,
                                                 133915, 79034, 80189, 110415, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 137515, 3, 82499, 82527, 111991, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 137560, 3, 82527, 82555, 112027, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 137605, 3, 82555, 82583, 112063, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 137650, 3, 82583, 82611, 112099, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 137695, 3, 82611, 82639, 112135, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 137740, 3, 82639, 82667, 112171, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 137785, 3, 82667, 82695, 112207, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 137830, 0, 3, 111955, 137515, 112351,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 137965, 0, 3, 111991, 137560, 112459,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 138100, 0, 3, 112027, 137605, 112567,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 138235, 0, 3, 112063, 137650, 112675,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 138370, 0, 3, 112099, 137695, 112783,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 138505, 0, 3, 112135, 137740, 112891,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 138640, 0, 3, 112171, 137785, 112999,
                                                 ncols, p);

            compute_prim_dl_electron_repulsion_0(buffer, 138775, 0, 3, 112243, 137830, 83507,
                                                 83675, 113323, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 139045, 0, 3, 112351, 137965, 83675,
                                                 83843, 113539, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 139315, 0, 3, 112459, 138100, 83843,
                                                 84011, 113755, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 139585, 0, 3, 112567, 138235, 84011,
                                                 84179, 113971, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 139855, 0, 3, 112675, 138370, 84179,
                                                 84347, 114187, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 140125, 0, 3, 112783, 138505, 84347,
                                                 84515, 114403, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 140395, 0, 3, 112891, 138640, 84515,
                                                 84683, 114619, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 140665, 0, 3, 113107, 138775, 85019,
                                                 85299, 114835, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 141115, 0, 3, 113323, 139045, 85299,
                                                 85579, 115195, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 141565, 0, 3, 113539, 139315, 85579,
                                                 85859, 115555, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 142015, 0, 3, 113755, 139585, 85859,
                                                 86139, 115915, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 142465, 0, 3, 113971, 139855, 86139,
                                                 86419, 116275, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 142915, 0, 3, 114187, 140125, 86419,
                                                 86699, 116635, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 143365, 0, 3, 114403, 140395, 86699,
                                                 86979, 116995, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 143815, 0, 3, 138775, 139045, 115195,
                                                 141565, 87539, 87959, 117895, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 144490, 0, 3, 139045, 139315, 115555,
                                                 142015, 87959, 88379, 118435, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 145165, 0, 3, 139315, 139585, 115915,
                                                 142465, 88379, 88799, 118975, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 145840, 0, 3, 139585, 139855, 116275,
                                                 142915, 88799, 89219, 119515, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 146515, 0, 3, 139855, 140125, 116635,
                                                 143365, 89219, 89639, 120055, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 147190, 0, 3, 140665, 141115, 117355,
                                                 143815, 90479, 91067, 120595, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 148135, 0, 3, 141115, 141565, 117895,
                                                 144490, 91067, 91655, 121351, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 149080, 0, 3, 141565, 142015, 118435,
                                                 145165, 91655, 92243, 122107, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 150025, 0, 3, 142015, 142465, 118975,
                                                 145840, 92243, 92831, 122863, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 150970, 0, 3, 142465, 142915, 119515,
                                                 146515, 92831, 93419, 123619, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 151915, 0, 3, 143815, 144490, 121351,
                                                 149080, 94595, 95379, 125383, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 153175, 0, 3, 144490, 145165, 122107,
                                                 150025, 95379, 96163, 126391, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 154435, 0, 3, 145165, 145840, 122863,
                                                 150970, 96163, 96947, 127399, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 155695, 0, 3, 147190, 148135, 124375,
                                                 151915, 98515, 99523, 128407, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 157315, 0, 3, 148135, 149080, 125383,
                                                 153175, 99523, 100531, 129703, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 158935, 0, 3, 149080, 150025, 126391,
                                                 154435, 100531, 101539, 130999, ncols, alpha,
                                                 beta, p);

            compute_prim_ll_electron_repulsion_0(buffer, 160555, 0, 3, 151915, 153175, 129703,
                                                 158935, 103555, 104815, 133915, ncols, alpha,
                                                 beta, p);

            compute_prim_ml_electron_repulsion_0(buffer, 162580, 0, 3, 155695, 157315, 132295,
                                                 160555, 107335, 108875, 135535, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_l_x(buffer, 165055, 155695, 162580, 1, 45, ncols, alpha);

            simdgeo::geom_l_y(buffer, 167080, 155695, 162580, 1, 45, ncols, alpha);

            simdgeo::geom_l_z(buffer, 169105, 155695, 162580, 1, 45, ncols, alpha);

            simdfunc::contract_primitives(buffer, 171130, 165055, 6075, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 177205, 171130, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 177205, 17, nmax);

    simdtrf::transform_l_inner(buffer, 177205, 173155, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 289 * nvalues, nvalues, buffer, 177205, 17, nmax);

    simdtrf::transform_l_inner(buffer, 177205, 175180, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 578 * nvalues, nvalues, buffer, 177205, 17, nmax);
}

}  // namespace simdt2ceri
