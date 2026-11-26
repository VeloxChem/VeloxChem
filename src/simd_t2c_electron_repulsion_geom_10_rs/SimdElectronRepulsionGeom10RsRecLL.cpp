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


#include "SimdElectronRepulsionGeom10RsRecLL.hpp"

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
compute_rs_geom_10_ll_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ll_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 355169, 342254, 12150, nvalues);

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
                                                9, 10, 11, 12, 13, 14, 15, 16, 17}, ncols, fj,
                                                mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 24, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13, 14, 15, 16, 17}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 84, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 87, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 90, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 93, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 96, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 99, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 102, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 105, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 108, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 111, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 114, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 117, 0, 33, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 120, 0, 34, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 123, 0, 35, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 126, 0, 36, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 129, 0, 37, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 132, 0, 38, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 135, 0, 39, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 138, 0, 40, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 141, 0, 41, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 7, 8, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 8, 9, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 9, 10, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 10, 11, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 11, 12, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 12, 13, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 13, 14, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 14, 15, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 192, 0, 15, 16, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 198, 0, 16, 17, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 204, 0, 17, 18, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 210, 0, 18, 19, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 216, 0, 19, 20, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 222, 0, 20, 21, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 228, 0, 21, 22, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 234, 0, 25, 26, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 240, 0, 26, 27, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 246, 0, 27, 28, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 252, 0, 28, 29, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 258, 0, 29, 30, 111, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 264, 0, 30, 31, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 270, 0, 31, 32, 117, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 276, 0, 32, 33, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 282, 0, 33, 34, 123, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 288, 0, 34, 35, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 294, 0, 35, 36, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 300, 0, 36, 37, 132, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 306, 0, 37, 38, 135, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 312, 0, 38, 39, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 318, 0, 39, 40, 141, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 324, 0, 42, 45, 144, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 334, 0, 45, 48, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 344, 0, 48, 51, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 354, 0, 51, 54, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 364, 0, 54, 57, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 374, 0, 57, 60, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 384, 0, 60, 63, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 394, 0, 63, 66, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 404, 0, 66, 69, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 414, 0, 69, 72, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 424, 0, 72, 75, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 434, 0, 75, 78, 210, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 444, 0, 78, 81, 216, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 454, 0, 81, 84, 222, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 464, 0, 84, 87, 228, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 474, 0, 93, 96, 234, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 484, 0, 96, 99, 240, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 494, 0, 99, 102, 246, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 504, 0, 102, 105, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 514, 0, 105, 108, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 524, 0, 108, 111, 264, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 534, 0, 111, 114, 270, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 544, 0, 114, 117, 276, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 554, 0, 117, 120, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 564, 0, 120, 123, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 574, 0, 123, 126, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 584, 0, 126, 129, 300, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 594, 0, 129, 132, 306, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 604, 0, 132, 135, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 614, 0, 135, 138, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 624, 0, 144, 150, 344, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 639, 0, 150, 156, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 654, 0, 156, 162, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 669, 0, 162, 168, 374, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 684, 0, 168, 174, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 699, 0, 174, 180, 394, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 714, 0, 180, 186, 404, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 729, 0, 186, 192, 414, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 744, 0, 192, 198, 424, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 759, 0, 198, 204, 434, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 774, 0, 204, 210, 444, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 789, 0, 210, 216, 454, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 804, 0, 216, 222, 464, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 819, 0, 234, 240, 494, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 834, 0, 240, 246, 504, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 849, 0, 246, 252, 514, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 864, 0, 252, 258, 524, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 879, 0, 258, 264, 534, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 894, 0, 264, 270, 544, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 909, 0, 270, 276, 554, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 924, 0, 276, 282, 564, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 939, 0, 282, 288, 574, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 954, 0, 288, 294, 584, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 969, 0, 294, 300, 594, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 984, 0, 300, 306, 604, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 999, 0, 306, 312, 614, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1014, 0, 324, 334, 624, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1035, 0, 334, 344, 639, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1056, 0, 344, 354, 654, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1077, 0, 354, 364, 669, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1098, 0, 364, 374, 684, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1119, 0, 374, 384, 699, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1140, 0, 384, 394, 714, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1161, 0, 394, 404, 729, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1182, 0, 404, 414, 744, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1203, 0, 414, 424, 759, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1224, 0, 424, 434, 774, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1245, 0, 434, 444, 789, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1266, 0, 444, 454, 804, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1287, 0, 474, 484, 819, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1308, 0, 484, 494, 834, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1329, 0, 494, 504, 849, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1350, 0, 504, 514, 864, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1371, 0, 514, 524, 879, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1392, 0, 524, 534, 894, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1413, 0, 534, 544, 909, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1434, 0, 544, 554, 924, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1455, 0, 554, 564, 939, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1476, 0, 564, 574, 954, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1497, 0, 574, 584, 969, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1518, 0, 584, 594, 984, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1539, 0, 594, 604, 999, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1560, 0, 624, 639, 1056, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1588, 0, 639, 654, 1077, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1616, 0, 654, 669, 1098, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1644, 0, 669, 684, 1119, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1672, 0, 684, 699, 1140, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1700, 0, 699, 714, 1161, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1728, 0, 714, 729, 1182, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1756, 0, 729, 744, 1203, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1784, 0, 744, 759, 1224, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1812, 0, 759, 774, 1245, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1840, 0, 774, 789, 1266, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1868, 0, 819, 834, 1329, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1896, 0, 834, 849, 1350, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1924, 0, 849, 864, 1371, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1952, 0, 864, 879, 1392, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1980, 0, 879, 894, 1413, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 2008, 0, 894, 909, 1434, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 2036, 0, 909, 924, 1455, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 2064, 0, 924, 939, 1476, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 2092, 0, 939, 954, 1497, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 2120, 0, 954, 969, 1518, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 2148, 0, 969, 984, 1539, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2176, 0, 1014, 1035, 1560, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2212, 0, 1035, 1056, 1588, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2248, 0, 1056, 1077, 1616, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2284, 0, 1077, 1098, 1644, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2320, 0, 1098, 1119, 1672, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2356, 0, 1119, 1140, 1700, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2392, 0, 1140, 1161, 1728, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2428, 0, 1161, 1182, 1756, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2464, 0, 1182, 1203, 1784, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2500, 0, 1203, 1224, 1812, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2536, 0, 1224, 1245, 1840, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2572, 0, 1287, 1308, 1868, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2608, 0, 1308, 1329, 1896, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2644, 0, 1329, 1350, 1924, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2680, 0, 1350, 1371, 1952, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2716, 0, 1371, 1392, 1980, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2752, 0, 1392, 1413, 2008, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2788, 0, 1413, 1434, 2036, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2824, 0, 1434, 1455, 2064, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2860, 0, 1455, 1476, 2092, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2896, 0, 1476, 1497, 2120, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2932, 0, 1497, 1518, 2148, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2968, 0, 1560, 1588, 2248, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3013, 0, 1588, 1616, 2284, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3058, 0, 1616, 1644, 2320, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3103, 0, 1644, 1672, 2356, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3148, 0, 1672, 1700, 2392, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3193, 0, 1700, 1728, 2428, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3238, 0, 1728, 1756, 2464, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3283, 0, 1756, 1784, 2500, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3328, 0, 1784, 1812, 2536, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3373, 0, 1868, 1896, 2644, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3418, 0, 1896, 1924, 2680, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3463, 0, 1924, 1952, 2716, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3508, 0, 1952, 1980, 2752, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3553, 0, 1980, 2008, 2788, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3598, 0, 2008, 2036, 2824, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3643, 0, 2036, 2064, 2860, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3688, 0, 2064, 2092, 2896, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3733, 0, 2092, 2120, 2932, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3778, 0, 2176, 2212, 2968, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3833, 0, 2212, 2248, 3013, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3888, 0, 2248, 2284, 3058, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3943, 0, 2284, 2320, 3103, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3998, 0, 2320, 2356, 3148, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4053, 0, 2356, 2392, 3193, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4108, 0, 2392, 2428, 3238, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4163, 0, 2428, 2464, 3283, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4218, 0, 2464, 2500, 3328, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4273, 0, 2572, 2608, 3373, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4328, 0, 2608, 2644, 3418, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4383, 0, 2644, 2680, 3463, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4438, 0, 2680, 2716, 3508, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4493, 0, 2716, 2752, 3553, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4548, 0, 2752, 2788, 3598, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4603, 0, 2788, 2824, 3643, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4658, 0, 2824, 2860, 3688, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4713, 0, 2860, 2896, 3733, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 4768, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4771, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4774, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4777, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4780, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4783, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4786, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4789, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4792, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4795, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4798, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4801, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4804, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4807, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4810, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4813, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4816, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4819, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4822, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4825, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4828, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4831, 3, 35, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4834, 3, 36, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4837, 3, 37, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4840, 3, 38, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4843, 3, 39, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4846, 3, 40, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4849, 3, 41, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 4852, 3, 9, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4861, 3, 10, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4870, 3, 11, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4879, 3, 12, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4888, 3, 13, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4897, 3, 14, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4906, 3, 15, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4915, 3, 16, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4924, 3, 17, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4933, 3, 18, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4942, 3, 19, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4951, 3, 20, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4960, 3, 21, 87, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4969, 3, 22, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4978, 3, 27, 102, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4987, 3, 28, 105, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4996, 3, 29, 108, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5005, 3, 30, 111, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5014, 3, 31, 114, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5023, 3, 32, 117, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5032, 3, 33, 120, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5041, 3, 34, 123, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5050, 3, 35, 126, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5059, 3, 36, 129, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5068, 3, 37, 132, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5077, 3, 38, 135, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5086, 3, 39, 138, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 5095, 3, 40, 141, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5104, 0, 3, 48, 4852, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5122, 0, 3, 51, 4861, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5140, 0, 3, 54, 4870, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5158, 0, 3, 57, 4879, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5176, 0, 3, 60, 4888, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5194, 0, 3, 63, 4897, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5212, 0, 3, 66, 4906, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5230, 0, 3, 69, 4915, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5248, 0, 3, 72, 4924, 198, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5266, 0, 3, 75, 4933, 204, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5284, 0, 3, 78, 4942, 210, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5302, 0, 3, 81, 4951, 216, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5320, 0, 3, 84, 4960, 222, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5338, 0, 3, 87, 4969, 228, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5356, 0, 3, 99, 4978, 240, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5374, 0, 3, 102, 4987, 246, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5392, 0, 3, 105, 4996, 252, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5410, 0, 3, 108, 5005, 258, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5428, 0, 3, 111, 5014, 264, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5446, 0, 3, 114, 5023, 270, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5464, 0, 3, 117, 5032, 276, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5482, 0, 3, 120, 5041, 282, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5500, 0, 3, 123, 5050, 288, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5518, 0, 3, 126, 5059, 294, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5536, 0, 3, 129, 5068, 300, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5554, 0, 3, 132, 5077, 306, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5572, 0, 3, 135, 5086, 312, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 5590, 0, 3, 138, 5095, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5608, 0, 3, 150, 5122, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5638, 0, 3, 156, 5140, 354, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5668, 0, 3, 162, 5158, 364, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5698, 0, 3, 168, 5176, 374, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5728, 0, 3, 174, 5194, 384, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5758, 0, 3, 180, 5212, 394, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5788, 0, 3, 186, 5230, 404, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5818, 0, 3, 192, 5248, 414, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5848, 0, 3, 198, 5266, 424, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5878, 0, 3, 204, 5284, 434, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5908, 0, 3, 210, 5302, 444, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5938, 0, 3, 216, 5320, 454, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5968, 0, 3, 222, 5338, 464, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5998, 0, 3, 240, 5374, 494, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6028, 0, 3, 246, 5392, 504, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6058, 0, 3, 252, 5410, 514, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6088, 0, 3, 258, 5428, 524, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6118, 0, 3, 264, 5446, 534, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6148, 0, 3, 270, 5464, 544, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6178, 0, 3, 276, 5482, 554, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6208, 0, 3, 282, 5500, 564, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6238, 0, 3, 288, 5518, 574, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6268, 0, 3, 294, 5536, 584, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6298, 0, 3, 300, 5554, 594, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6328, 0, 3, 306, 5572, 604, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 6358, 0, 3, 312, 5590, 614, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6388, 0, 3, 344, 5638, 639, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6433, 0, 3, 354, 5668, 654, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6478, 0, 3, 364, 5698, 669, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6523, 0, 3, 374, 5728, 684, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6568, 0, 3, 384, 5758, 699, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6613, 0, 3, 394, 5788, 714, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6658, 0, 3, 404, 5818, 729, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6703, 0, 3, 414, 5848, 744, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6748, 0, 3, 424, 5878, 759, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6793, 0, 3, 434, 5908, 774, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6838, 0, 3, 444, 5938, 789, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6883, 0, 3, 454, 5968, 804, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6928, 0, 3, 494, 6028, 834, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6973, 0, 3, 504, 6058, 849, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7018, 0, 3, 514, 6088, 864, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7063, 0, 3, 524, 6118, 879, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7108, 0, 3, 534, 6148, 894, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7153, 0, 3, 544, 6178, 909, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7198, 0, 3, 554, 6208, 924, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7243, 0, 3, 564, 6238, 939, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7288, 0, 3, 574, 6268, 954, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7333, 0, 3, 584, 6298, 969, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7378, 0, 3, 594, 6328, 984, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 7423, 0, 3, 604, 6358, 999, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7468, 0, 3, 639, 6433, 1056, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7531, 0, 3, 654, 6478, 1077, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7594, 0, 3, 669, 6523, 1098, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7657, 0, 3, 684, 6568, 1119, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7720, 0, 3, 699, 6613, 1140, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7783, 0, 3, 714, 6658, 1161, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7846, 0, 3, 729, 6703, 1182, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7909, 0, 3, 744, 6748, 1203, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7972, 0, 3, 759, 6793, 1224, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8035, 0, 3, 774, 6838, 1245, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8098, 0, 3, 789, 6883, 1266, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8161, 0, 3, 834, 6973, 1329, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8224, 0, 3, 849, 7018, 1350, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8287, 0, 3, 864, 7063, 1371, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8350, 0, 3, 879, 7108, 1392, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8413, 0, 3, 894, 7153, 1413, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8476, 0, 3, 909, 7198, 1434, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8539, 0, 3, 924, 7243, 1455, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8602, 0, 3, 939, 7288, 1476, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8665, 0, 3, 954, 7333, 1497, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8728, 0, 3, 969, 7378, 1518, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 8791, 0, 3, 984, 7423, 1539, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 8854, 0, 3, 1056, 7531, 1588, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8938, 0, 3, 1077, 7594, 1616, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9022, 0, 3, 1098, 7657, 1644, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9106, 0, 3, 1119, 7720, 1672, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9190, 0, 3, 1140, 7783, 1700, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9274, 0, 3, 1161, 7846, 1728, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9358, 0, 3, 1182, 7909, 1756, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9442, 0, 3, 1203, 7972, 1784, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9526, 0, 3, 1224, 8035, 1812, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9610, 0, 3, 1245, 8098, 1840, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9694, 0, 3, 1329, 8224, 1896, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9778, 0, 3, 1350, 8287, 1924, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9862, 0, 3, 1371, 8350, 1952, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9946, 0, 3, 1392, 8413, 1980, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 10030, 0, 3, 1413, 8476, 2008, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 10114, 0, 3, 1434, 8539, 2036, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 10198, 0, 3, 1455, 8602, 2064, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 10282, 0, 3, 1476, 8665, 2092, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 10366, 0, 3, 1497, 8728, 2120, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 10450, 0, 3, 1518, 8791, 2148, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10534, 0, 3, 1588, 8938, 2248, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10642, 0, 3, 1616, 9022, 2284, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10750, 0, 3, 1644, 9106, 2320, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10858, 0, 3, 1672, 9190, 2356, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10966, 0, 3, 1700, 9274, 2392, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11074, 0, 3, 1728, 9358, 2428, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11182, 0, 3, 1756, 9442, 2464, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11290, 0, 3, 1784, 9526, 2500, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11398, 0, 3, 1812, 9610, 2536, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11506, 0, 3, 1896, 9778, 2644, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11614, 0, 3, 1924, 9862, 2680, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11722, 0, 3, 1952, 9946, 2716, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11830, 0, 3, 1980, 10030, 2752, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11938, 0, 3, 2008, 10114, 2788, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 12046, 0, 3, 2036, 10198, 2824, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 12154, 0, 3, 2064, 10282, 2860, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 12262, 0, 3, 2092, 10366, 2896, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 12370, 0, 3, 2120, 10450, 2932, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12478, 0, 3, 2248, 10642, 3013, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12613, 0, 3, 2284, 10750, 3058, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12748, 0, 3, 2320, 10858, 3103, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12883, 0, 3, 2356, 10966, 3148, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13018, 0, 3, 2392, 11074, 3193, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13153, 0, 3, 2428, 11182, 3238, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13288, 0, 3, 2464, 11290, 3283, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13423, 0, 3, 2500, 11398, 3328, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13558, 0, 3, 2644, 11614, 3418, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13693, 0, 3, 2680, 11722, 3463, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13828, 0, 3, 2716, 11830, 3508, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13963, 0, 3, 2752, 11938, 3553, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 14098, 0, 3, 2788, 12046, 3598, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 14233, 0, 3, 2824, 12154, 3643, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 14368, 0, 3, 2860, 12262, 3688, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 14503, 0, 3, 2896, 12370, 3733, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14638, 0, 3, 3013, 12613, 3888, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14803, 0, 3, 3058, 12748, 3943, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14968, 0, 3, 3103, 12883, 3998, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15133, 0, 3, 3148, 13018, 4053, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15298, 0, 3, 3193, 13153, 4108, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15463, 0, 3, 3238, 13288, 4163, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15628, 0, 3, 3283, 13423, 4218, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15793, 0, 3, 3418, 13693, 4383, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15958, 0, 3, 3463, 13828, 4438, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 16123, 0, 3, 3508, 13963, 4493, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 16288, 0, 3, 3553, 14098, 4548, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 16453, 0, 3, 3598, 14233, 4603, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 16618, 0, 3, 3643, 14368, 4658, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 16783, 0, 3, 3688, 14503, 4713, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 16948, 3, 9, 10, 4771, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 16954, 3, 10, 11, 4774, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 16960, 3, 11, 12, 4777, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 16966, 3, 12, 13, 4780, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 16972, 3, 13, 14, 4783, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 16978, 3, 14, 15, 4786, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 16984, 3, 15, 16, 4789, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 16990, 3, 16, 17, 4792, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 16996, 3, 17, 18, 4795, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17002, 3, 18, 19, 4798, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17008, 3, 19, 20, 4801, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17014, 3, 20, 21, 4804, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17020, 3, 21, 22, 4807, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17026, 3, 27, 28, 4813, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17032, 3, 28, 29, 4816, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17038, 3, 29, 30, 4819, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17044, 3, 30, 31, 4822, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17050, 3, 31, 32, 4825, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17056, 3, 32, 33, 4828, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17062, 3, 33, 34, 4831, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17068, 3, 34, 35, 4834, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17074, 3, 35, 36, 4837, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17080, 3, 36, 37, 4840, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17086, 3, 37, 38, 4843, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17092, 3, 38, 39, 4846, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 17098, 3, 39, 40, 4849, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 17104, 0, 3, 4768, 16948, 4861, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17122, 0, 3, 4771, 16954, 4870, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17140, 0, 3, 4774, 16960, 4879, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17158, 0, 3, 4777, 16966, 4888, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17176, 0, 3, 4780, 16972, 4897, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17194, 0, 3, 4783, 16978, 4906, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17212, 0, 3, 4786, 16984, 4915, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17230, 0, 3, 4789, 16990, 4924, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17248, 0, 3, 4792, 16996, 4933, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17266, 0, 3, 4795, 17002, 4942, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17284, 0, 3, 4798, 17008, 4951, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17302, 0, 3, 4801, 17014, 4960, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17320, 0, 3, 4804, 17020, 4969, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17338, 0, 3, 4810, 17026, 4987, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17356, 0, 3, 4813, 17032, 4996, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17374, 0, 3, 4816, 17038, 5005, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17392, 0, 3, 4819, 17044, 5014, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17410, 0, 3, 4822, 17050, 5023, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17428, 0, 3, 4825, 17056, 5032, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17446, 0, 3, 4828, 17062, 5041, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17464, 0, 3, 4831, 17068, 5050, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17482, 0, 3, 4834, 17074, 5059, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17500, 0, 3, 4837, 17080, 5068, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17518, 0, 3, 4840, 17086, 5077, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17536, 0, 3, 4843, 17092, 5086, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 17554, 0, 3, 4846, 17098, 5095, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 17572, 0, 3, 4852, 17104, 144, 150,
                                                 5122, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17608, 0, 3, 4861, 17122, 150, 156,
                                                 5140, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17644, 0, 3, 4870, 17140, 156, 162,
                                                 5158, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17680, 0, 3, 4879, 17158, 162, 168,
                                                 5176, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17716, 0, 3, 4888, 17176, 168, 174,
                                                 5194, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17752, 0, 3, 4897, 17194, 174, 180,
                                                 5212, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17788, 0, 3, 4906, 17212, 180, 186,
                                                 5230, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17824, 0, 3, 4915, 17230, 186, 192,
                                                 5248, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17860, 0, 3, 4924, 17248, 192, 198,
                                                 5266, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17896, 0, 3, 4933, 17266, 198, 204,
                                                 5284, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17932, 0, 3, 4942, 17284, 204, 210,
                                                 5302, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17968, 0, 3, 4951, 17302, 210, 216,
                                                 5320, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18004, 0, 3, 4960, 17320, 216, 222,
                                                 5338, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18040, 0, 3, 4978, 17338, 234, 240,
                                                 5374, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18076, 0, 3, 4987, 17356, 240, 246,
                                                 5392, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18112, 0, 3, 4996, 17374, 246, 252,
                                                 5410, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18148, 0, 3, 5005, 17392, 252, 258,
                                                 5428, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18184, 0, 3, 5014, 17410, 258, 264,
                                                 5446, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18220, 0, 3, 5023, 17428, 264, 270,
                                                 5464, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18256, 0, 3, 5032, 17446, 270, 276,
                                                 5482, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18292, 0, 3, 5041, 17464, 276, 282,
                                                 5500, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18328, 0, 3, 5050, 17482, 282, 288,
                                                 5518, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18364, 0, 3, 5059, 17500, 288, 294,
                                                 5536, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18400, 0, 3, 5068, 17518, 294, 300,
                                                 5554, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18436, 0, 3, 5077, 17536, 300, 306,
                                                 5572, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 18472, 0, 3, 5086, 17554, 306, 312,
                                                 5590, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18508, 0, 3, 5104, 17572, 324, 334,
                                                 5608, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18568, 0, 3, 5122, 17608, 334, 344,
                                                 5638, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18628, 0, 3, 5140, 17644, 344, 354,
                                                 5668, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18688, 0, 3, 5158, 17680, 354, 364,
                                                 5698, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18748, 0, 3, 5176, 17716, 364, 374,
                                                 5728, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18808, 0, 3, 5194, 17752, 374, 384,
                                                 5758, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18868, 0, 3, 5212, 17788, 384, 394,
                                                 5788, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18928, 0, 3, 5230, 17824, 394, 404,
                                                 5818, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18988, 0, 3, 5248, 17860, 404, 414,
                                                 5848, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19048, 0, 3, 5266, 17896, 414, 424,
                                                 5878, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19108, 0, 3, 5284, 17932, 424, 434,
                                                 5908, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19168, 0, 3, 5302, 17968, 434, 444,
                                                 5938, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19228, 0, 3, 5320, 18004, 444, 454,
                                                 5968, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19288, 0, 3, 5356, 18040, 474, 484,
                                                 5998, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19348, 0, 3, 5374, 18076, 484, 494,
                                                 6028, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19408, 0, 3, 5392, 18112, 494, 504,
                                                 6058, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19468, 0, 3, 5410, 18148, 504, 514,
                                                 6088, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19528, 0, 3, 5428, 18184, 514, 524,
                                                 6118, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19588, 0, 3, 5446, 18220, 524, 534,
                                                 6148, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19648, 0, 3, 5464, 18256, 534, 544,
                                                 6178, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19708, 0, 3, 5482, 18292, 544, 554,
                                                 6208, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19768, 0, 3, 5500, 18328, 554, 564,
                                                 6238, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19828, 0, 3, 5518, 18364, 564, 574,
                                                 6268, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19888, 0, 3, 5536, 18400, 574, 584,
                                                 6298, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 19948, 0, 3, 5554, 18436, 584, 594,
                                                 6328, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 20008, 0, 3, 5572, 18472, 594, 604,
                                                 6358, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20068, 0, 3, 17572, 17608, 5638, 18628,
                                                 624, 639, 6433, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20158, 0, 3, 17608, 17644, 5668, 18688,
                                                 639, 654, 6478, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20248, 0, 3, 17644, 17680, 5698, 18748,
                                                 654, 669, 6523, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20338, 0, 3, 17680, 17716, 5728, 18808,
                                                 669, 684, 6568, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20428, 0, 3, 17716, 17752, 5758, 18868,
                                                 684, 699, 6613, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20518, 0, 3, 17752, 17788, 5788, 18928,
                                                 699, 714, 6658, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20608, 0, 3, 17788, 17824, 5818, 18988,
                                                 714, 729, 6703, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20698, 0, 3, 17824, 17860, 5848, 19048,
                                                 729, 744, 6748, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20788, 0, 3, 17860, 17896, 5878, 19108,
                                                 744, 759, 6793, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20878, 0, 3, 17896, 17932, 5908, 19168,
                                                 759, 774, 6838, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20968, 0, 3, 17932, 17968, 5938, 19228,
                                                 774, 789, 6883, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21058, 0, 3, 18040, 18076, 6028, 19408,
                                                 819, 834, 6973, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21148, 0, 3, 18076, 18112, 6058, 19468,
                                                 834, 849, 7018, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21238, 0, 3, 18112, 18148, 6088, 19528,
                                                 849, 864, 7063, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21328, 0, 3, 18148, 18184, 6118, 19588,
                                                 864, 879, 7108, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21418, 0, 3, 18184, 18220, 6148, 19648,
                                                 879, 894, 7153, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21508, 0, 3, 18220, 18256, 6178, 19708,
                                                 894, 909, 7198, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21598, 0, 3, 18256, 18292, 6208, 19768,
                                                 909, 924, 7243, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21688, 0, 3, 18292, 18328, 6238, 19828,
                                                 924, 939, 7288, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21778, 0, 3, 18328, 18364, 6268, 19888,
                                                 939, 954, 7333, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21868, 0, 3, 18364, 18400, 6298, 19948,
                                                 954, 969, 7378, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 21958, 0, 3, 18400, 18436, 6328, 20008,
                                                 969, 984, 7423, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22048, 0, 3, 18508, 18568, 6388, 20068,
                                                 1014, 1035, 7468, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22174, 0, 3, 18568, 18628, 6433, 20158,
                                                 1035, 1056, 7531, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22300, 0, 3, 18628, 18688, 6478, 20248,
                                                 1056, 1077, 7594, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22426, 0, 3, 18688, 18748, 6523, 20338,
                                                 1077, 1098, 7657, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22552, 0, 3, 18748, 18808, 6568, 20428,
                                                 1098, 1119, 7720, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22678, 0, 3, 18808, 18868, 6613, 20518,
                                                 1119, 1140, 7783, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22804, 0, 3, 18868, 18928, 6658, 20608,
                                                 1140, 1161, 7846, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22930, 0, 3, 18928, 18988, 6703, 20698,
                                                 1161, 1182, 7909, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 23056, 0, 3, 18988, 19048, 6748, 20788,
                                                 1182, 1203, 7972, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 23182, 0, 3, 19048, 19108, 6793, 20878,
                                                 1203, 1224, 8035, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 23308, 0, 3, 19108, 19168, 6838, 20968,
                                                 1224, 1245, 8098, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 23434, 0, 3, 19288, 19348, 6928, 21058,
                                                 1287, 1308, 8161, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 23560, 0, 3, 19348, 19408, 6973, 21148,
                                                 1308, 1329, 8224, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 23686, 0, 3, 19408, 19468, 7018, 21238,
                                                 1329, 1350, 8287, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 23812, 0, 3, 19468, 19528, 7063, 21328,
                                                 1350, 1371, 8350, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 23938, 0, 3, 19528, 19588, 7108, 21418,
                                                 1371, 1392, 8413, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 24064, 0, 3, 19588, 19648, 7153, 21508,
                                                 1392, 1413, 8476, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 24190, 0, 3, 19648, 19708, 7198, 21598,
                                                 1413, 1434, 8539, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 24316, 0, 3, 19708, 19768, 7243, 21688,
                                                 1434, 1455, 8602, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 24442, 0, 3, 19768, 19828, 7288, 21778,
                                                 1455, 1476, 8665, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 24568, 0, 3, 19828, 19888, 7333, 21868,
                                                 1476, 1497, 8728, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 24694, 0, 3, 19888, 19948, 7378, 21958,
                                                 1497, 1518, 8791, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 24820, 0, 3, 20068, 20158, 7531, 22300,
                                                 1560, 1588, 8938, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 24988, 0, 3, 20158, 20248, 7594, 22426,
                                                 1588, 1616, 9022, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 25156, 0, 3, 20248, 20338, 7657, 22552,
                                                 1616, 1644, 9106, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 25324, 0, 3, 20338, 20428, 7720, 22678,
                                                 1644, 1672, 9190, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 25492, 0, 3, 20428, 20518, 7783, 22804,
                                                 1672, 1700, 9274, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 25660, 0, 3, 20518, 20608, 7846, 22930,
                                                 1700, 1728, 9358, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 25828, 0, 3, 20608, 20698, 7909, 23056,
                                                 1728, 1756, 9442, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 25996, 0, 3, 20698, 20788, 7972, 23182,
                                                 1756, 1784, 9526, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 26164, 0, 3, 20788, 20878, 8035, 23308,
                                                 1784, 1812, 9610, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 26332, 0, 3, 21058, 21148, 8224, 23686,
                                                 1868, 1896, 9778, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 26500, 0, 3, 21148, 21238, 8287, 23812,
                                                 1896, 1924, 9862, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 26668, 0, 3, 21238, 21328, 8350, 23938,
                                                 1924, 1952, 9946, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 26836, 0, 3, 21328, 21418, 8413, 24064,
                                                 1952, 1980, 10030, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 27004, 0, 3, 21418, 21508, 8476, 24190,
                                                 1980, 2008, 10114, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 27172, 0, 3, 21508, 21598, 8539, 24316,
                                                 2008, 2036, 10198, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 27340, 0, 3, 21598, 21688, 8602, 24442,
                                                 2036, 2064, 10282, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 27508, 0, 3, 21688, 21778, 8665, 24568,
                                                 2064, 2092, 10366, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 27676, 0, 3, 21778, 21868, 8728, 24694,
                                                 2092, 2120, 10450, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 27844, 0, 3, 22048, 22174, 8854, 24820,
                                                 2176, 2212, 10534, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 28060, 0, 3, 22174, 22300, 8938, 24988,
                                                 2212, 2248, 10642, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 28276, 0, 3, 22300, 22426, 9022, 25156,
                                                 2248, 2284, 10750, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 28492, 0, 3, 22426, 22552, 9106, 25324,
                                                 2284, 2320, 10858, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 28708, 0, 3, 22552, 22678, 9190, 25492,
                                                 2320, 2356, 10966, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 28924, 0, 3, 22678, 22804, 9274, 25660,
                                                 2356, 2392, 11074, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 29140, 0, 3, 22804, 22930, 9358, 25828,
                                                 2392, 2428, 11182, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 29356, 0, 3, 22930, 23056, 9442, 25996,
                                                 2428, 2464, 11290, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 29572, 0, 3, 23056, 23182, 9526, 26164,
                                                 2464, 2500, 11398, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 29788, 0, 3, 23434, 23560, 9694, 26332,
                                                 2572, 2608, 11506, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 30004, 0, 3, 23560, 23686, 9778, 26500,
                                                 2608, 2644, 11614, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 30220, 0, 3, 23686, 23812, 9862, 26668,
                                                 2644, 2680, 11722, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 30436, 0, 3, 23812, 23938, 9946, 26836,
                                                 2680, 2716, 11830, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 30652, 0, 3, 23938, 24064, 10030, 27004,
                                                 2716, 2752, 11938, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 30868, 0, 3, 24064, 24190, 10114, 27172,
                                                 2752, 2788, 12046, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 31084, 0, 3, 24190, 24316, 10198, 27340,
                                                 2788, 2824, 12154, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 31300, 0, 3, 24316, 24442, 10282, 27508,
                                                 2824, 2860, 12262, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 31516, 0, 3, 24442, 24568, 10366, 27676,
                                                 2860, 2896, 12370, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 31732, 0, 3, 24820, 24988, 10642, 28276,
                                                 2968, 3013, 12613, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 32002, 0, 3, 24988, 25156, 10750, 28492,
                                                 3013, 3058, 12748, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 32272, 0, 3, 25156, 25324, 10858, 28708,
                                                 3058, 3103, 12883, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 32542, 0, 3, 25324, 25492, 10966, 28924,
                                                 3103, 3148, 13018, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 32812, 0, 3, 25492, 25660, 11074, 29140,
                                                 3148, 3193, 13153, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 33082, 0, 3, 25660, 25828, 11182, 29356,
                                                 3193, 3238, 13288, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 33352, 0, 3, 25828, 25996, 11290, 29572,
                                                 3238, 3283, 13423, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 33622, 0, 3, 26332, 26500, 11614, 30220,
                                                 3373, 3418, 13693, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 33892, 0, 3, 26500, 26668, 11722, 30436,
                                                 3418, 3463, 13828, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 34162, 0, 3, 26668, 26836, 11830, 30652,
                                                 3463, 3508, 13963, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 34432, 0, 3, 26836, 27004, 11938, 30868,
                                                 3508, 3553, 14098, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 34702, 0, 3, 27004, 27172, 12046, 31084,
                                                 3553, 3598, 14233, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 34972, 0, 3, 27172, 27340, 12154, 31300,
                                                 3598, 3643, 14368, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 35242, 0, 3, 27340, 27508, 12262, 31516,
                                                 3643, 3688, 14503, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 35512, 0, 3, 27844, 28060, 12478, 31732,
                                                 3778, 3833, 14638, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 35842, 0, 3, 28060, 28276, 12613, 32002,
                                                 3833, 3888, 14803, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 36172, 0, 3, 28276, 28492, 12748, 32272,
                                                 3888, 3943, 14968, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 36502, 0, 3, 28492, 28708, 12883, 32542,
                                                 3943, 3998, 15133, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 36832, 0, 3, 28708, 28924, 13018, 32812,
                                                 3998, 4053, 15298, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 37162, 0, 3, 28924, 29140, 13153, 33082,
                                                 4053, 4108, 15463, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 37492, 0, 3, 29140, 29356, 13288, 33352,
                                                 4108, 4163, 15628, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 37822, 0, 3, 29788, 30004, 13558, 33622,
                                                 4273, 4328, 15793, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 38152, 0, 3, 30004, 30220, 13693, 33892,
                                                 4328, 4383, 15958, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 38482, 0, 3, 30220, 30436, 13828, 34162,
                                                 4383, 4438, 16123, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 38812, 0, 3, 30436, 30652, 13963, 34432,
                                                 4438, 4493, 16288, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 39142, 0, 3, 30652, 30868, 14098, 34702,
                                                 4493, 4548, 16453, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 39472, 0, 3, 30868, 31084, 14233, 34972,
                                                 4548, 4603, 16618, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 39802, 0, 3, 31084, 31300, 14368, 35242,
                                                 4603, 4658, 16783, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40132, 3, 4768, 4771, 16954, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40142, 3, 4771, 4774, 16960, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40152, 3, 4774, 4777, 16966, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40162, 3, 4777, 4780, 16972, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40172, 3, 4780, 4783, 16978, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40182, 3, 4783, 4786, 16984, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40192, 3, 4786, 4789, 16990, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40202, 3, 4789, 4792, 16996, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40212, 3, 4792, 4795, 17002, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40222, 3, 4795, 4798, 17008, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40232, 3, 4798, 4801, 17014, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40242, 3, 4801, 4804, 17020, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40252, 3, 4810, 4813, 17032, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40262, 3, 4813, 4816, 17038, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40272, 3, 4816, 4819, 17044, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40282, 3, 4819, 4822, 17050, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40292, 3, 4822, 4825, 17056, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40302, 3, 4825, 4828, 17062, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40312, 3, 4828, 4831, 17068, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40322, 3, 4831, 4834, 17074, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40332, 3, 4834, 4837, 17080, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40342, 3, 4837, 4840, 17086, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40352, 3, 4840, 4843, 17092, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 40362, 3, 4843, 4846, 17098, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 40372, 0, 3, 16948, 40132, 17122, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40402, 0, 3, 16954, 40142, 17140, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40432, 0, 3, 16960, 40152, 17158, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40462, 0, 3, 16966, 40162, 17176, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40492, 0, 3, 16972, 40172, 17194, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40522, 0, 3, 16978, 40182, 17212, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40552, 0, 3, 16984, 40192, 17230, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40582, 0, 3, 16990, 40202, 17248, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40612, 0, 3, 16996, 40212, 17266, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40642, 0, 3, 17002, 40222, 17284, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40672, 0, 3, 17008, 40232, 17302, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40702, 0, 3, 17014, 40242, 17320, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40732, 0, 3, 17026, 40252, 17356, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40762, 0, 3, 17032, 40262, 17374, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40792, 0, 3, 17038, 40272, 17392, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40822, 0, 3, 17044, 40282, 17410, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40852, 0, 3, 17050, 40292, 17428, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40882, 0, 3, 17056, 40302, 17446, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40912, 0, 3, 17062, 40312, 17464, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40942, 0, 3, 17068, 40322, 17482, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 40972, 0, 3, 17074, 40332, 17500, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 41002, 0, 3, 17080, 40342, 17518, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 41032, 0, 3, 17086, 40352, 17536, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 41062, 0, 3, 17092, 40362, 17554, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 41092, 0, 3, 17104, 40372, 5104, 5122,
                                                 17608, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41152, 0, 3, 17122, 40402, 5122, 5140,
                                                 17644, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41212, 0, 3, 17140, 40432, 5140, 5158,
                                                 17680, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41272, 0, 3, 17158, 40462, 5158, 5176,
                                                 17716, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41332, 0, 3, 17176, 40492, 5176, 5194,
                                                 17752, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41392, 0, 3, 17194, 40522, 5194, 5212,
                                                 17788, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41452, 0, 3, 17212, 40552, 5212, 5230,
                                                 17824, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41512, 0, 3, 17230, 40582, 5230, 5248,
                                                 17860, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41572, 0, 3, 17248, 40612, 5248, 5266,
                                                 17896, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41632, 0, 3, 17266, 40642, 5266, 5284,
                                                 17932, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41692, 0, 3, 17284, 40672, 5284, 5302,
                                                 17968, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41752, 0, 3, 17302, 40702, 5302, 5320,
                                                 18004, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41812, 0, 3, 17338, 40732, 5356, 5374,
                                                 18076, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41872, 0, 3, 17356, 40762, 5374, 5392,
                                                 18112, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41932, 0, 3, 17374, 40792, 5392, 5410,
                                                 18148, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 41992, 0, 3, 17392, 40822, 5410, 5428,
                                                 18184, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 42052, 0, 3, 17410, 40852, 5428, 5446,
                                                 18220, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 42112, 0, 3, 17428, 40882, 5446, 5464,
                                                 18256, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 42172, 0, 3, 17446, 40912, 5464, 5482,
                                                 18292, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 42232, 0, 3, 17464, 40942, 5482, 5500,
                                                 18328, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 42292, 0, 3, 17482, 40972, 5500, 5518,
                                                 18364, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 42352, 0, 3, 17500, 41002, 5518, 5536,
                                                 18400, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 42412, 0, 3, 17518, 41032, 5536, 5554,
                                                 18436, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 42472, 0, 3, 17536, 41062, 5554, 5572,
                                                 18472, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 42532, 0, 3, 17608, 41152, 5608, 5638,
                                                 18628, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 42632, 0, 3, 17644, 41212, 5638, 5668,
                                                 18688, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 42732, 0, 3, 17680, 41272, 5668, 5698,
                                                 18748, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 42832, 0, 3, 17716, 41332, 5698, 5728,
                                                 18808, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 42932, 0, 3, 17752, 41392, 5728, 5758,
                                                 18868, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43032, 0, 3, 17788, 41452, 5758, 5788,
                                                 18928, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43132, 0, 3, 17824, 41512, 5788, 5818,
                                                 18988, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43232, 0, 3, 17860, 41572, 5818, 5848,
                                                 19048, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43332, 0, 3, 17896, 41632, 5848, 5878,
                                                 19108, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43432, 0, 3, 17932, 41692, 5878, 5908,
                                                 19168, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43532, 0, 3, 17968, 41752, 5908, 5938,
                                                 19228, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43632, 0, 3, 18076, 41872, 5998, 6028,
                                                 19408, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43732, 0, 3, 18112, 41932, 6028, 6058,
                                                 19468, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43832, 0, 3, 18148, 41992, 6058, 6088,
                                                 19528, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 43932, 0, 3, 18184, 42052, 6088, 6118,
                                                 19588, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 44032, 0, 3, 18220, 42112, 6118, 6148,
                                                 19648, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 44132, 0, 3, 18256, 42172, 6148, 6178,
                                                 19708, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 44232, 0, 3, 18292, 42232, 6178, 6208,
                                                 19768, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 44332, 0, 3, 18328, 42292, 6208, 6238,
                                                 19828, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 44432, 0, 3, 18364, 42352, 6238, 6268,
                                                 19888, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 44532, 0, 3, 18400, 42412, 6268, 6298,
                                                 19948, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 44632, 0, 3, 18436, 42472, 6298, 6328,
                                                 20008, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 44732, 0, 3, 41092, 41152, 18628, 42632,
                                                 6388, 6433, 20158, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 44882, 0, 3, 41152, 41212, 18688, 42732,
                                                 6433, 6478, 20248, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 45032, 0, 3, 41212, 41272, 18748, 42832,
                                                 6478, 6523, 20338, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 45182, 0, 3, 41272, 41332, 18808, 42932,
                                                 6523, 6568, 20428, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 45332, 0, 3, 41332, 41392, 18868, 43032,
                                                 6568, 6613, 20518, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 45482, 0, 3, 41392, 41452, 18928, 43132,
                                                 6613, 6658, 20608, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 45632, 0, 3, 41452, 41512, 18988, 43232,
                                                 6658, 6703, 20698, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 45782, 0, 3, 41512, 41572, 19048, 43332,
                                                 6703, 6748, 20788, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 45932, 0, 3, 41572, 41632, 19108, 43432,
                                                 6748, 6793, 20878, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 46082, 0, 3, 41632, 41692, 19168, 43532,
                                                 6793, 6838, 20968, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 46232, 0, 3, 41812, 41872, 19408, 43732,
                                                 6928, 6973, 21148, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 46382, 0, 3, 41872, 41932, 19468, 43832,
                                                 6973, 7018, 21238, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 46532, 0, 3, 41932, 41992, 19528, 43932,
                                                 7018, 7063, 21328, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 46682, 0, 3, 41992, 42052, 19588, 44032,
                                                 7063, 7108, 21418, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 46832, 0, 3, 42052, 42112, 19648, 44132,
                                                 7108, 7153, 21508, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 46982, 0, 3, 42112, 42172, 19708, 44232,
                                                 7153, 7198, 21598, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 47132, 0, 3, 42172, 42232, 19768, 44332,
                                                 7198, 7243, 21688, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 47282, 0, 3, 42232, 42292, 19828, 44432,
                                                 7243, 7288, 21778, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 47432, 0, 3, 42292, 42352, 19888, 44532,
                                                 7288, 7333, 21868, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 47582, 0, 3, 42352, 42412, 19948, 44632,
                                                 7333, 7378, 21958, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 47732, 0, 3, 42532, 42632, 20158, 44882,
                                                 7468, 7531, 22300, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 47942, 0, 3, 42632, 42732, 20248, 45032,
                                                 7531, 7594, 22426, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 48152, 0, 3, 42732, 42832, 20338, 45182,
                                                 7594, 7657, 22552, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 48362, 0, 3, 42832, 42932, 20428, 45332,
                                                 7657, 7720, 22678, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 48572, 0, 3, 42932, 43032, 20518, 45482,
                                                 7720, 7783, 22804, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 48782, 0, 3, 43032, 43132, 20608, 45632,
                                                 7783, 7846, 22930, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 48992, 0, 3, 43132, 43232, 20698, 45782,
                                                 7846, 7909, 23056, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 49202, 0, 3, 43232, 43332, 20788, 45932,
                                                 7909, 7972, 23182, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 49412, 0, 3, 43332, 43432, 20878, 46082,
                                                 7972, 8035, 23308, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 49622, 0, 3, 43632, 43732, 21148, 46382,
                                                 8161, 8224, 23686, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 49832, 0, 3, 43732, 43832, 21238, 46532,
                                                 8224, 8287, 23812, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 50042, 0, 3, 43832, 43932, 21328, 46682,
                                                 8287, 8350, 23938, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 50252, 0, 3, 43932, 44032, 21418, 46832,
                                                 8350, 8413, 24064, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 50462, 0, 3, 44032, 44132, 21508, 46982,
                                                 8413, 8476, 24190, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 50672, 0, 3, 44132, 44232, 21598, 47132,
                                                 8476, 8539, 24316, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 50882, 0, 3, 44232, 44332, 21688, 47282,
                                                 8539, 8602, 24442, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 51092, 0, 3, 44332, 44432, 21778, 47432,
                                                 8602, 8665, 24568, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 51302, 0, 3, 44432, 44532, 21868, 47582,
                                                 8665, 8728, 24694, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 51512, 0, 3, 44732, 44882, 22300, 47942,
                                                 8854, 8938, 24988, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 51792, 0, 3, 44882, 45032, 22426, 48152,
                                                 8938, 9022, 25156, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 52072, 0, 3, 45032, 45182, 22552, 48362,
                                                 9022, 9106, 25324, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 52352, 0, 3, 45182, 45332, 22678, 48572,
                                                 9106, 9190, 25492, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 52632, 0, 3, 45332, 45482, 22804, 48782,
                                                 9190, 9274, 25660, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 52912, 0, 3, 45482, 45632, 22930, 48992,
                                                 9274, 9358, 25828, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 53192, 0, 3, 45632, 45782, 23056, 49202,
                                                 9358, 9442, 25996, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 53472, 0, 3, 45782, 45932, 23182, 49412,
                                                 9442, 9526, 26164, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 53752, 0, 3, 46232, 46382, 23686, 49832,
                                                 9694, 9778, 26500, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 54032, 0, 3, 46382, 46532, 23812, 50042,
                                                 9778, 9862, 26668, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 54312, 0, 3, 46532, 46682, 23938, 50252,
                                                 9862, 9946, 26836, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 54592, 0, 3, 46682, 46832, 24064, 50462,
                                                 9946, 10030, 27004, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 54872, 0, 3, 46832, 46982, 24190, 50672,
                                                 10030, 10114, 27172, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 55152, 0, 3, 46982, 47132, 24316, 50882,
                                                 10114, 10198, 27340, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 55432, 0, 3, 47132, 47282, 24442, 51092,
                                                 10198, 10282, 27508, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 55712, 0, 3, 47282, 47432, 24568, 51302,
                                                 10282, 10366, 27676, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 55992, 0, 3, 47732, 47942, 24988, 51792,
                                                 10534, 10642, 28276, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 56352, 0, 3, 47942, 48152, 25156, 52072,
                                                 10642, 10750, 28492, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 56712, 0, 3, 48152, 48362, 25324, 52352,
                                                 10750, 10858, 28708, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 57072, 0, 3, 48362, 48572, 25492, 52632,
                                                 10858, 10966, 28924, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 57432, 0, 3, 48572, 48782, 25660, 52912,
                                                 10966, 11074, 29140, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 57792, 0, 3, 48782, 48992, 25828, 53192,
                                                 11074, 11182, 29356, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 58152, 0, 3, 48992, 49202, 25996, 53472,
                                                 11182, 11290, 29572, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 58512, 0, 3, 49622, 49832, 26500, 54032,
                                                 11506, 11614, 30220, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 58872, 0, 3, 49832, 50042, 26668, 54312,
                                                 11614, 11722, 30436, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 59232, 0, 3, 50042, 50252, 26836, 54592,
                                                 11722, 11830, 30652, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 59592, 0, 3, 50252, 50462, 27004, 54872,
                                                 11830, 11938, 30868, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 59952, 0, 3, 50462, 50672, 27172, 55152,
                                                 11938, 12046, 31084, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 60312, 0, 3, 50672, 50882, 27340, 55432,
                                                 12046, 12154, 31300, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 60672, 0, 3, 50882, 51092, 27508, 55712,
                                                 12154, 12262, 31516, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 61032, 0, 3, 51512, 51792, 28276, 56352,
                                                 12478, 12613, 32002, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 61482, 0, 3, 51792, 52072, 28492, 56712,
                                                 12613, 12748, 32272, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 61932, 0, 3, 52072, 52352, 28708, 57072,
                                                 12748, 12883, 32542, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 62382, 0, 3, 52352, 52632, 28924, 57432,
                                                 12883, 13018, 32812, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 62832, 0, 3, 52632, 52912, 29140, 57792,
                                                 13018, 13153, 33082, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 63282, 0, 3, 52912, 53192, 29356, 58152,
                                                 13153, 13288, 33352, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 63732, 0, 3, 53752, 54032, 30220, 58872,
                                                 13558, 13693, 33892, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 64182, 0, 3, 54032, 54312, 30436, 59232,
                                                 13693, 13828, 34162, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 64632, 0, 3, 54312, 54592, 30652, 59592,
                                                 13828, 13963, 34432, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 65082, 0, 3, 54592, 54872, 30868, 59952,
                                                 13963, 14098, 34702, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 65532, 0, 3, 54872, 55152, 31084, 60312,
                                                 14098, 14233, 34972, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 65982, 0, 3, 55152, 55432, 31300, 60672,
                                                 14233, 14368, 35242, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 66432, 0, 3, 55992, 56352, 32002, 61482,
                                                 14638, 14803, 36172, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 66982, 0, 3, 56352, 56712, 32272, 61932,
                                                 14803, 14968, 36502, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 67532, 0, 3, 56712, 57072, 32542, 62382,
                                                 14968, 15133, 36832, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 68082, 0, 3, 57072, 57432, 32812, 62832,
                                                 15133, 15298, 37162, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 68632, 0, 3, 57432, 57792, 33082, 63282,
                                                 15298, 15463, 37492, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 69182, 0, 3, 58512, 58872, 33892, 64182,
                                                 15793, 15958, 38482, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 69732, 0, 3, 58872, 59232, 34162, 64632,
                                                 15958, 16123, 38812, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 70282, 0, 3, 59232, 59592, 34432, 65082,
                                                 16123, 16288, 39142, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 70832, 0, 3, 59592, 59952, 34702, 65532,
                                                 16288, 16453, 39472, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 71382, 0, 3, 59952, 60312, 34972, 65982,
                                                 16453, 16618, 39802, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 71932, 3, 16948, 16954, 40142, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 71947, 3, 16954, 16960, 40152, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 71962, 3, 16960, 16966, 40162, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 71977, 3, 16966, 16972, 40172, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 71992, 3, 16972, 16978, 40182, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72007, 3, 16978, 16984, 40192, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72022, 3, 16984, 16990, 40202, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72037, 3, 16990, 16996, 40212, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72052, 3, 16996, 17002, 40222, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72067, 3, 17002, 17008, 40232, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72082, 3, 17008, 17014, 40242, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72097, 3, 17026, 17032, 40262, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72112, 3, 17032, 17038, 40272, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72127, 3, 17038, 17044, 40282, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72142, 3, 17044, 17050, 40292, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72157, 3, 17050, 17056, 40302, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72172, 3, 17056, 17062, 40312, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72187, 3, 17062, 17068, 40322, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72202, 3, 17068, 17074, 40332, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72217, 3, 17074, 17080, 40342, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72232, 3, 17080, 17086, 40352, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 72247, 3, 17086, 17092, 40362, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 72262, 0, 3, 40132, 71932, 40402, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72307, 0, 3, 40142, 71947, 40432, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72352, 0, 3, 40152, 71962, 40462, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72397, 0, 3, 40162, 71977, 40492, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72442, 0, 3, 40172, 71992, 40522, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72487, 0, 3, 40182, 72007, 40552, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72532, 0, 3, 40192, 72022, 40582, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72577, 0, 3, 40202, 72037, 40612, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72622, 0, 3, 40212, 72052, 40642, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72667, 0, 3, 40222, 72067, 40672, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72712, 0, 3, 40232, 72082, 40702, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72757, 0, 3, 40252, 72097, 40762, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72802, 0, 3, 40262, 72112, 40792, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72847, 0, 3, 40272, 72127, 40822, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72892, 0, 3, 40282, 72142, 40852, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72937, 0, 3, 40292, 72157, 40882, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 72982, 0, 3, 40302, 72172, 40912, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 73027, 0, 3, 40312, 72187, 40942, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 73072, 0, 3, 40322, 72202, 40972, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 73117, 0, 3, 40332, 72217, 41002, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 73162, 0, 3, 40342, 72232, 41032, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 73207, 0, 3, 40352, 72247, 41062, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 73252, 0, 3, 40372, 72262, 17572, 17608,
                                                 41152, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 73342, 0, 3, 40402, 72307, 17608, 17644,
                                                 41212, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 73432, 0, 3, 40432, 72352, 17644, 17680,
                                                 41272, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 73522, 0, 3, 40462, 72397, 17680, 17716,
                                                 41332, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 73612, 0, 3, 40492, 72442, 17716, 17752,
                                                 41392, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 73702, 0, 3, 40522, 72487, 17752, 17788,
                                                 41452, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 73792, 0, 3, 40552, 72532, 17788, 17824,
                                                 41512, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 73882, 0, 3, 40582, 72577, 17824, 17860,
                                                 41572, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 73972, 0, 3, 40612, 72622, 17860, 17896,
                                                 41632, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74062, 0, 3, 40642, 72667, 17896, 17932,
                                                 41692, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74152, 0, 3, 40672, 72712, 17932, 17968,
                                                 41752, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74242, 0, 3, 40732, 72757, 18040, 18076,
                                                 41872, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74332, 0, 3, 40762, 72802, 18076, 18112,
                                                 41932, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74422, 0, 3, 40792, 72847, 18112, 18148,
                                                 41992, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74512, 0, 3, 40822, 72892, 18148, 18184,
                                                 42052, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74602, 0, 3, 40852, 72937, 18184, 18220,
                                                 42112, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74692, 0, 3, 40882, 72982, 18220, 18256,
                                                 42172, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74782, 0, 3, 40912, 73027, 18256, 18292,
                                                 42232, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74872, 0, 3, 40942, 73072, 18292, 18328,
                                                 42292, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 74962, 0, 3, 40972, 73117, 18328, 18364,
                                                 42352, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 75052, 0, 3, 41002, 73162, 18364, 18400,
                                                 42412, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 75142, 0, 3, 41032, 73207, 18400, 18436,
                                                 42472, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 75232, 0, 3, 41092, 73252, 18508, 18568,
                                                 42532, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 75382, 0, 3, 41152, 73342, 18568, 18628,
                                                 42632, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 75532, 0, 3, 41212, 73432, 18628, 18688,
                                                 42732, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 75682, 0, 3, 41272, 73522, 18688, 18748,
                                                 42832, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 75832, 0, 3, 41332, 73612, 18748, 18808,
                                                 42932, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 75982, 0, 3, 41392, 73702, 18808, 18868,
                                                 43032, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 76132, 0, 3, 41452, 73792, 18868, 18928,
                                                 43132, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 76282, 0, 3, 41512, 73882, 18928, 18988,
                                                 43232, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 76432, 0, 3, 41572, 73972, 18988, 19048,
                                                 43332, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 76582, 0, 3, 41632, 74062, 19048, 19108,
                                                 43432, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 76732, 0, 3, 41692, 74152, 19108, 19168,
                                                 43532, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 76882, 0, 3, 41812, 74242, 19288, 19348,
                                                 43632, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 77032, 0, 3, 41872, 74332, 19348, 19408,
                                                 43732, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 77182, 0, 3, 41932, 74422, 19408, 19468,
                                                 43832, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 77332, 0, 3, 41992, 74512, 19468, 19528,
                                                 43932, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 77482, 0, 3, 42052, 74602, 19528, 19588,
                                                 44032, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 77632, 0, 3, 42112, 74692, 19588, 19648,
                                                 44132, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 77782, 0, 3, 42172, 74782, 19648, 19708,
                                                 44232, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 77932, 0, 3, 42232, 74872, 19708, 19768,
                                                 44332, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 78082, 0, 3, 42292, 74962, 19768, 19828,
                                                 44432, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 78232, 0, 3, 42352, 75052, 19828, 19888,
                                                 44532, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 78382, 0, 3, 42412, 75142, 19888, 19948,
                                                 44632, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 78532, 0, 3, 73252, 73342, 42632, 75532,
                                                 20068, 20158, 44882, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 78757, 0, 3, 73342, 73432, 42732, 75682,
                                                 20158, 20248, 45032, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 78982, 0, 3, 73432, 73522, 42832, 75832,
                                                 20248, 20338, 45182, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 79207, 0, 3, 73522, 73612, 42932, 75982,
                                                 20338, 20428, 45332, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 79432, 0, 3, 73612, 73702, 43032, 76132,
                                                 20428, 20518, 45482, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 79657, 0, 3, 73702, 73792, 43132, 76282,
                                                 20518, 20608, 45632, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 79882, 0, 3, 73792, 73882, 43232, 76432,
                                                 20608, 20698, 45782, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 80107, 0, 3, 73882, 73972, 43332, 76582,
                                                 20698, 20788, 45932, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 80332, 0, 3, 73972, 74062, 43432, 76732,
                                                 20788, 20878, 46082, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 80557, 0, 3, 74242, 74332, 43732, 77182,
                                                 21058, 21148, 46382, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 80782, 0, 3, 74332, 74422, 43832, 77332,
                                                 21148, 21238, 46532, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 81007, 0, 3, 74422, 74512, 43932, 77482,
                                                 21238, 21328, 46682, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 81232, 0, 3, 74512, 74602, 44032, 77632,
                                                 21328, 21418, 46832, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 81457, 0, 3, 74602, 74692, 44132, 77782,
                                                 21418, 21508, 46982, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 81682, 0, 3, 74692, 74782, 44232, 77932,
                                                 21508, 21598, 47132, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 81907, 0, 3, 74782, 74872, 44332, 78082,
                                                 21598, 21688, 47282, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 82132, 0, 3, 74872, 74962, 44432, 78232,
                                                 21688, 21778, 47432, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 82357, 0, 3, 74962, 75052, 44532, 78382,
                                                 21778, 21868, 47582, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 82582, 0, 3, 75232, 75382, 44732, 78532,
                                                 22048, 22174, 47732, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 82897, 0, 3, 75382, 75532, 44882, 78757,
                                                 22174, 22300, 47942, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 83212, 0, 3, 75532, 75682, 45032, 78982,
                                                 22300, 22426, 48152, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 83527, 0, 3, 75682, 75832, 45182, 79207,
                                                 22426, 22552, 48362, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 83842, 0, 3, 75832, 75982, 45332, 79432,
                                                 22552, 22678, 48572, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 84157, 0, 3, 75982, 76132, 45482, 79657,
                                                 22678, 22804, 48782, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 84472, 0, 3, 76132, 76282, 45632, 79882,
                                                 22804, 22930, 48992, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 84787, 0, 3, 76282, 76432, 45782, 80107,
                                                 22930, 23056, 49202, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 85102, 0, 3, 76432, 76582, 45932, 80332,
                                                 23056, 23182, 49412, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 85417, 0, 3, 76882, 77032, 46232, 80557,
                                                 23434, 23560, 49622, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 85732, 0, 3, 77032, 77182, 46382, 80782,
                                                 23560, 23686, 49832, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 86047, 0, 3, 77182, 77332, 46532, 81007,
                                                 23686, 23812, 50042, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 86362, 0, 3, 77332, 77482, 46682, 81232,
                                                 23812, 23938, 50252, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 86677, 0, 3, 77482, 77632, 46832, 81457,
                                                 23938, 24064, 50462, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 86992, 0, 3, 77632, 77782, 46982, 81682,
                                                 24064, 24190, 50672, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 87307, 0, 3, 77782, 77932, 47132, 81907,
                                                 24190, 24316, 50882, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 87622, 0, 3, 77932, 78082, 47282, 82132,
                                                 24316, 24442, 51092, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 87937, 0, 3, 78082, 78232, 47432, 82357,
                                                 24442, 24568, 51302, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 88252, 0, 3, 78532, 78757, 47942, 83212,
                                                 24820, 24988, 51792, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 88672, 0, 3, 78757, 78982, 48152, 83527,
                                                 24988, 25156, 52072, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 89092, 0, 3, 78982, 79207, 48362, 83842,
                                                 25156, 25324, 52352, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 89512, 0, 3, 79207, 79432, 48572, 84157,
                                                 25324, 25492, 52632, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 89932, 0, 3, 79432, 79657, 48782, 84472,
                                                 25492, 25660, 52912, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 90352, 0, 3, 79657, 79882, 48992, 84787,
                                                 25660, 25828, 53192, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 90772, 0, 3, 79882, 80107, 49202, 85102,
                                                 25828, 25996, 53472, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 91192, 0, 3, 80557, 80782, 49832, 86047,
                                                 26332, 26500, 54032, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 91612, 0, 3, 80782, 81007, 50042, 86362,
                                                 26500, 26668, 54312, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 92032, 0, 3, 81007, 81232, 50252, 86677,
                                                 26668, 26836, 54592, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 92452, 0, 3, 81232, 81457, 50462, 86992,
                                                 26836, 27004, 54872, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 92872, 0, 3, 81457, 81682, 50672, 87307,
                                                 27004, 27172, 55152, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 93292, 0, 3, 81682, 81907, 50882, 87622,
                                                 27172, 27340, 55432, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 93712, 0, 3, 81907, 82132, 51092, 87937,
                                                 27340, 27508, 55712, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 94132, 0, 3, 82582, 82897, 51512, 88252,
                                                 27844, 28060, 55992, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 94672, 0, 3, 82897, 83212, 51792, 88672,
                                                 28060, 28276, 56352, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 95212, 0, 3, 83212, 83527, 52072, 89092,
                                                 28276, 28492, 56712, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 95752, 0, 3, 83527, 83842, 52352, 89512,
                                                 28492, 28708, 57072, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 96292, 0, 3, 83842, 84157, 52632, 89932,
                                                 28708, 28924, 57432, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 96832, 0, 3, 84157, 84472, 52912, 90352,
                                                 28924, 29140, 57792, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 97372, 0, 3, 84472, 84787, 53192, 90772,
                                                 29140, 29356, 58152, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 97912, 0, 3, 85417, 85732, 53752, 91192,
                                                 29788, 30004, 58512, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 98452, 0, 3, 85732, 86047, 54032, 91612,
                                                 30004, 30220, 58872, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 98992, 0, 3, 86047, 86362, 54312, 92032,
                                                 30220, 30436, 59232, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 99532, 0, 3, 86362, 86677, 54592, 92452,
                                                 30436, 30652, 59592, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 100072, 0, 3, 86677, 86992, 54872,
                                                 92872, 30652, 30868, 59952, ncols, alpha, beta,
                                                 p);

            compute_prim_kg_electron_repulsion_0(buffer, 100612, 0, 3, 86992, 87307, 55152,
                                                 93292, 30868, 31084, 60312, ncols, alpha, beta,
                                                 p);

            compute_prim_kg_electron_repulsion_0(buffer, 101152, 0, 3, 87307, 87622, 55432,
                                                 93712, 31084, 31300, 60672, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 101692, 0, 3, 88252, 88672, 56352,
                                                 95212, 31732, 32002, 61482, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 102367, 0, 3, 88672, 89092, 56712,
                                                 95752, 32002, 32272, 61932, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 103042, 0, 3, 89092, 89512, 57072,
                                                 96292, 32272, 32542, 62382, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 103717, 0, 3, 89512, 89932, 57432,
                                                 96832, 32542, 32812, 62832, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 104392, 0, 3, 89932, 90352, 57792,
                                                 97372, 32812, 33082, 63282, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 105067, 0, 3, 91192, 91612, 58872,
                                                 98992, 33622, 33892, 64182, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 105742, 0, 3, 91612, 92032, 59232,
                                                 99532, 33892, 34162, 64632, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 106417, 0, 3, 92032, 92452, 59592,
                                                 100072, 34162, 34432, 65082, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 107092, 0, 3, 92452, 92872, 59952,
                                                 100612, 34432, 34702, 65532, ncols, alpha, beta,
                                                 p);

            compute_prim_lg_electron_repulsion_0(buffer, 107767, 0, 3, 92872, 93292, 60312,
                                                 101152, 34702, 34972, 65982, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 108442, 0, 3, 94132, 94672, 61032,
                                                 101692, 35512, 35842, 66432, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 109267, 0, 3, 94672, 95212, 61482,
                                                 102367, 35842, 36172, 66982, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 110092, 0, 3, 95212, 95752, 61932,
                                                 103042, 36172, 36502, 67532, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 110917, 0, 3, 95752, 96292, 62382,
                                                 103717, 36502, 36832, 68082, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 111742, 0, 3, 96292, 96832, 62832,
                                                 104392, 36832, 37162, 68632, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 112567, 0, 3, 97912, 98452, 63732,
                                                 105067, 37822, 38152, 69182, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 113392, 0, 3, 98452, 98992, 64182,
                                                 105742, 38152, 38482, 69732, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 114217, 0, 3, 98992, 99532, 64632,
                                                 106417, 38482, 38812, 70282, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 115042, 0, 3, 99532, 100072, 65082,
                                                 107092, 38812, 39142, 70832, ncols, alpha, beta,
                                                 p);

            compute_prim_mg_electron_repulsion_0(buffer, 115867, 0, 3, 100072, 100612, 65532,
                                                 107767, 39142, 39472, 71382, ncols, alpha, beta,
                                                 p);

            compute_prim_sh_electron_repulsion_0(buffer, 116692, 3, 40132, 40142, 71947, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116713, 3, 40142, 40152, 71962, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116734, 3, 40152, 40162, 71977, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116755, 3, 40162, 40172, 71992, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116776, 3, 40172, 40182, 72007, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116797, 3, 40182, 40192, 72022, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116818, 3, 40192, 40202, 72037, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116839, 3, 40202, 40212, 72052, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116860, 3, 40212, 40222, 72067, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116881, 3, 40222, 40232, 72082, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116902, 3, 40252, 40262, 72112, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116923, 3, 40262, 40272, 72127, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116944, 3, 40272, 40282, 72142, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116965, 3, 40282, 40292, 72157, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 116986, 3, 40292, 40302, 72172, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 117007, 3, 40302, 40312, 72187, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 117028, 3, 40312, 40322, 72202, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 117049, 3, 40322, 40332, 72217, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 117070, 3, 40332, 40342, 72232, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 117091, 3, 40342, 40352, 72247, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117112, 0, 3, 71932, 116692, 72307,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117175, 0, 3, 71947, 116713, 72352,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117238, 0, 3, 71962, 116734, 72397,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117301, 0, 3, 71977, 116755, 72442,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117364, 0, 3, 71992, 116776, 72487,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117427, 0, 3, 72007, 116797, 72532,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117490, 0, 3, 72022, 116818, 72577,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117553, 0, 3, 72037, 116839, 72622,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117616, 0, 3, 72052, 116860, 72667,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117679, 0, 3, 72067, 116881, 72712,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117742, 0, 3, 72097, 116902, 72802,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117805, 0, 3, 72112, 116923, 72847,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117868, 0, 3, 72127, 116944, 72892,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117931, 0, 3, 72142, 116965, 72937,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 117994, 0, 3, 72157, 116986, 72982,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 118057, 0, 3, 72172, 117007, 73027,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 118120, 0, 3, 72187, 117028, 73072,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 118183, 0, 3, 72202, 117049, 73117,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 118246, 0, 3, 72217, 117070, 73162,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 118309, 0, 3, 72232, 117091, 73207,
                                                 ncols, p);

            compute_prim_dh_electron_repulsion_0(buffer, 118372, 0, 3, 72262, 117112, 41092,
                                                 41152, 73342, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 118498, 0, 3, 72307, 117175, 41152,
                                                 41212, 73432, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 118624, 0, 3, 72352, 117238, 41212,
                                                 41272, 73522, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 118750, 0, 3, 72397, 117301, 41272,
                                                 41332, 73612, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 118876, 0, 3, 72442, 117364, 41332,
                                                 41392, 73702, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 119002, 0, 3, 72487, 117427, 41392,
                                                 41452, 73792, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 119128, 0, 3, 72532, 117490, 41452,
                                                 41512, 73882, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 119254, 0, 3, 72577, 117553, 41512,
                                                 41572, 73972, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 119380, 0, 3, 72622, 117616, 41572,
                                                 41632, 74062, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 119506, 0, 3, 72667, 117679, 41632,
                                                 41692, 74152, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 119632, 0, 3, 72757, 117742, 41812,
                                                 41872, 74332, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 119758, 0, 3, 72802, 117805, 41872,
                                                 41932, 74422, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 119884, 0, 3, 72847, 117868, 41932,
                                                 41992, 74512, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 120010, 0, 3, 72892, 117931, 41992,
                                                 42052, 74602, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 120136, 0, 3, 72937, 117994, 42052,
                                                 42112, 74692, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 120262, 0, 3, 72982, 118057, 42112,
                                                 42172, 74782, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 120388, 0, 3, 73027, 118120, 42172,
                                                 42232, 74872, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 120514, 0, 3, 73072, 118183, 42232,
                                                 42292, 74962, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 120640, 0, 3, 73117, 118246, 42292,
                                                 42352, 75052, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 120766, 0, 3, 73162, 118309, 42352,
                                                 42412, 75142, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 120892, 0, 3, 73342, 118498, 42532,
                                                 42632, 75532, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 121102, 0, 3, 73432, 118624, 42632,
                                                 42732, 75682, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 121312, 0, 3, 73522, 118750, 42732,
                                                 42832, 75832, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 121522, 0, 3, 73612, 118876, 42832,
                                                 42932, 75982, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 121732, 0, 3, 73702, 119002, 42932,
                                                 43032, 76132, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 121942, 0, 3, 73792, 119128, 43032,
                                                 43132, 76282, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 122152, 0, 3, 73882, 119254, 43132,
                                                 43232, 76432, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 122362, 0, 3, 73972, 119380, 43232,
                                                 43332, 76582, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 122572, 0, 3, 74062, 119506, 43332,
                                                 43432, 76732, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 122782, 0, 3, 74332, 119758, 43632,
                                                 43732, 77182, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 122992, 0, 3, 74422, 119884, 43732,
                                                 43832, 77332, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 123202, 0, 3, 74512, 120010, 43832,
                                                 43932, 77482, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 123412, 0, 3, 74602, 120136, 43932,
                                                 44032, 77632, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 123622, 0, 3, 74692, 120262, 44032,
                                                 44132, 77782, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 123832, 0, 3, 74782, 120388, 44132,
                                                 44232, 77932, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 124042, 0, 3, 74872, 120514, 44232,
                                                 44332, 78082, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 124252, 0, 3, 74962, 120640, 44332,
                                                 44432, 78232, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 124462, 0, 3, 75052, 120766, 44432,
                                                 44532, 78382, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 124672, 0, 3, 118372, 118498, 75532,
                                                 121102, 44732, 44882, 78757, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 124987, 0, 3, 118498, 118624, 75682,
                                                 121312, 44882, 45032, 78982, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 125302, 0, 3, 118624, 118750, 75832,
                                                 121522, 45032, 45182, 79207, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 125617, 0, 3, 118750, 118876, 75982,
                                                 121732, 45182, 45332, 79432, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 125932, 0, 3, 118876, 119002, 76132,
                                                 121942, 45332, 45482, 79657, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 126247, 0, 3, 119002, 119128, 76282,
                                                 122152, 45482, 45632, 79882, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 126562, 0, 3, 119128, 119254, 76432,
                                                 122362, 45632, 45782, 80107, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 126877, 0, 3, 119254, 119380, 76582,
                                                 122572, 45782, 45932, 80332, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 127192, 0, 3, 119632, 119758, 77182,
                                                 122992, 46232, 46382, 80782, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 127507, 0, 3, 119758, 119884, 77332,
                                                 123202, 46382, 46532, 81007, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 127822, 0, 3, 119884, 120010, 77482,
                                                 123412, 46532, 46682, 81232, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 128137, 0, 3, 120010, 120136, 77632,
                                                 123622, 46682, 46832, 81457, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 128452, 0, 3, 120136, 120262, 77782,
                                                 123832, 46832, 46982, 81682, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 128767, 0, 3, 120262, 120388, 77932,
                                                 124042, 46982, 47132, 81907, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 129082, 0, 3, 120388, 120514, 78082,
                                                 124252, 47132, 47282, 82132, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 129397, 0, 3, 120514, 120640, 78232,
                                                 124462, 47282, 47432, 82357, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 129712, 0, 3, 120892, 121102, 78757,
                                                 124987, 47732, 47942, 83212, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 130153, 0, 3, 121102, 121312, 78982,
                                                 125302, 47942, 48152, 83527, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 130594, 0, 3, 121312, 121522, 79207,
                                                 125617, 48152, 48362, 83842, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 131035, 0, 3, 121522, 121732, 79432,
                                                 125932, 48362, 48572, 84157, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 131476, 0, 3, 121732, 121942, 79657,
                                                 126247, 48572, 48782, 84472, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 131917, 0, 3, 121942, 122152, 79882,
                                                 126562, 48782, 48992, 84787, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 132358, 0, 3, 122152, 122362, 80107,
                                                 126877, 48992, 49202, 85102, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 132799, 0, 3, 122782, 122992, 80782,
                                                 127507, 49622, 49832, 86047, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 133240, 0, 3, 122992, 123202, 81007,
                                                 127822, 49832, 50042, 86362, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 133681, 0, 3, 123202, 123412, 81232,
                                                 128137, 50042, 50252, 86677, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 134122, 0, 3, 123412, 123622, 81457,
                                                 128452, 50252, 50462, 86992, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 134563, 0, 3, 123622, 123832, 81682,
                                                 128767, 50462, 50672, 87307, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 135004, 0, 3, 123832, 124042, 81907,
                                                 129082, 50672, 50882, 87622, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 135445, 0, 3, 124042, 124252, 82132,
                                                 129397, 50882, 51092, 87937, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 135886, 0, 3, 124672, 124987, 83212,
                                                 130153, 51512, 51792, 88672, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 136474, 0, 3, 124987, 125302, 83527,
                                                 130594, 51792, 52072, 89092, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 137062, 0, 3, 125302, 125617, 83842,
                                                 131035, 52072, 52352, 89512, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 137650, 0, 3, 125617, 125932, 84157,
                                                 131476, 52352, 52632, 89932, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 138238, 0, 3, 125932, 126247, 84472,
                                                 131917, 52632, 52912, 90352, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 138826, 0, 3, 126247, 126562, 84787,
                                                 132358, 52912, 53192, 90772, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 139414, 0, 3, 127192, 127507, 86047,
                                                 133240, 53752, 54032, 91612, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 140002, 0, 3, 127507, 127822, 86362,
                                                 133681, 54032, 54312, 92032, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 140590, 0, 3, 127822, 128137, 86677,
                                                 134122, 54312, 54592, 92452, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 141178, 0, 3, 128137, 128452, 86992,
                                                 134563, 54592, 54872, 92872, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 141766, 0, 3, 128452, 128767, 87307,
                                                 135004, 54872, 55152, 93292, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 142354, 0, 3, 128767, 129082, 87622,
                                                 135445, 55152, 55432, 93712, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 142942, 0, 3, 129712, 130153, 88672,
                                                 136474, 55992, 56352, 95212, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 143698, 0, 3, 130153, 130594, 89092,
                                                 137062, 56352, 56712, 95752, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 144454, 0, 3, 130594, 131035, 89512,
                                                 137650, 56712, 57072, 96292, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 145210, 0, 3, 131035, 131476, 89932,
                                                 138238, 57072, 57432, 96832, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 145966, 0, 3, 131476, 131917, 90352,
                                                 138826, 57432, 57792, 97372, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 146722, 0, 3, 132799, 133240, 91612,
                                                 140002, 58512, 58872, 98992, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 147478, 0, 3, 133240, 133681, 92032,
                                                 140590, 58872, 59232, 99532, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 148234, 0, 3, 133681, 134122, 92452,
                                                 141178, 59232, 59592, 100072, ncols, alpha,
                                                 beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 148990, 0, 3, 134122, 134563, 92872,
                                                 141766, 59592, 59952, 100612, ncols, alpha,
                                                 beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 149746, 0, 3, 134563, 135004, 93292,
                                                 142354, 59952, 60312, 101152, ncols, alpha,
                                                 beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 150502, 0, 3, 135886, 136474, 95212,
                                                 143698, 61032, 61482, 102367, ncols, alpha,
                                                 beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 151447, 0, 3, 136474, 137062, 95752,
                                                 144454, 61482, 61932, 103042, ncols, alpha,
                                                 beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 152392, 0, 3, 137062, 137650, 96292,
                                                 145210, 61932, 62382, 103717, ncols, alpha,
                                                 beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 153337, 0, 3, 137650, 138238, 96832,
                                                 145966, 62382, 62832, 104392, ncols, alpha,
                                                 beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 154282, 0, 3, 139414, 140002, 98992,
                                                 147478, 63732, 64182, 105742, ncols, alpha,
                                                 beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 155227, 0, 3, 140002, 140590, 99532,
                                                 148234, 64182, 64632, 106417, ncols, alpha,
                                                 beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 156172, 0, 3, 140590, 141178, 100072,
                                                 148990, 64632, 65082, 107092, ncols, alpha,
                                                 beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 157117, 0, 3, 141178, 141766, 100612,
                                                 149746, 65082, 65532, 107767, ncols, alpha,
                                                 beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 158062, 0, 3, 142942, 143698, 102367,
                                                 151447, 66432, 66982, 110092, ncols, alpha,
                                                 beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 159217, 0, 3, 143698, 144454, 103042,
                                                 152392, 66982, 67532, 110917, ncols, alpha,
                                                 beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 160372, 0, 3, 144454, 145210, 103717,
                                                 153337, 67532, 68082, 111742, ncols, alpha,
                                                 beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 161527, 0, 3, 146722, 147478, 105742,
                                                 155227, 69182, 69732, 114217, ncols, alpha,
                                                 beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 162682, 0, 3, 147478, 148234, 106417,
                                                 156172, 69732, 70282, 115042, ncols, alpha,
                                                 beta, p);

            compute_prim_mh_electron_repulsion_0(buffer, 163837, 0, 3, 148234, 148990, 107092,
                                                 157117, 70282, 70832, 115867, ncols, alpha,
                                                 beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 164992, 3, 71932, 71947, 116713, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165020, 3, 71947, 71962, 116734, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165048, 3, 71962, 71977, 116755, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165076, 3, 71977, 71992, 116776, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165104, 3, 71992, 72007, 116797, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165132, 3, 72007, 72022, 116818, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165160, 3, 72022, 72037, 116839, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165188, 3, 72037, 72052, 116860, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165216, 3, 72052, 72067, 116881, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165244, 3, 72097, 72112, 116923, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165272, 3, 72112, 72127, 116944, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165300, 3, 72127, 72142, 116965, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165328, 3, 72142, 72157, 116986, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165356, 3, 72157, 72172, 117007, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165384, 3, 72172, 72187, 117028, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165412, 3, 72187, 72202, 117049, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165440, 3, 72202, 72217, 117070, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 165468, 3, 72217, 72232, 117091, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 165496, 0, 3, 116692, 164992, 117175,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 165580, 0, 3, 116713, 165020, 117238,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 165664, 0, 3, 116734, 165048, 117301,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 165748, 0, 3, 116755, 165076, 117364,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 165832, 0, 3, 116776, 165104, 117427,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 165916, 0, 3, 116797, 165132, 117490,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166000, 0, 3, 116818, 165160, 117553,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166084, 0, 3, 116839, 165188, 117616,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166168, 0, 3, 116860, 165216, 117679,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166252, 0, 3, 116902, 165244, 117805,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166336, 0, 3, 116923, 165272, 117868,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166420, 0, 3, 116944, 165300, 117931,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166504, 0, 3, 116965, 165328, 117994,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166588, 0, 3, 116986, 165356, 118057,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166672, 0, 3, 117007, 165384, 118120,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166756, 0, 3, 117028, 165412, 118183,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166840, 0, 3, 117049, 165440, 118246,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 166924, 0, 3, 117070, 165468, 118309,
                                                 ncols, p);

            compute_prim_di_electron_repulsion_0(buffer, 167008, 0, 3, 117112, 165496, 73252,
                                                 73342, 118498, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 167176, 0, 3, 117175, 165580, 73342,
                                                 73432, 118624, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 167344, 0, 3, 117238, 165664, 73432,
                                                 73522, 118750, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 167512, 0, 3, 117301, 165748, 73522,
                                                 73612, 118876, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 167680, 0, 3, 117364, 165832, 73612,
                                                 73702, 119002, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 167848, 0, 3, 117427, 165916, 73702,
                                                 73792, 119128, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 168016, 0, 3, 117490, 166000, 73792,
                                                 73882, 119254, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 168184, 0, 3, 117553, 166084, 73882,
                                                 73972, 119380, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 168352, 0, 3, 117616, 166168, 73972,
                                                 74062, 119506, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 168520, 0, 3, 117742, 166252, 74242,
                                                 74332, 119758, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 168688, 0, 3, 117805, 166336, 74332,
                                                 74422, 119884, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 168856, 0, 3, 117868, 166420, 74422,
                                                 74512, 120010, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 169024, 0, 3, 117931, 166504, 74512,
                                                 74602, 120136, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 169192, 0, 3, 117994, 166588, 74602,
                                                 74692, 120262, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 169360, 0, 3, 118057, 166672, 74692,
                                                 74782, 120388, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 169528, 0, 3, 118120, 166756, 74782,
                                                 74872, 120514, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 169696, 0, 3, 118183, 166840, 74872,
                                                 74962, 120640, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 169864, 0, 3, 118246, 166924, 74962,
                                                 75052, 120766, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 170032, 0, 3, 118372, 167008, 75232,
                                                 75382, 120892, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 170312, 0, 3, 118498, 167176, 75382,
                                                 75532, 121102, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 170592, 0, 3, 118624, 167344, 75532,
                                                 75682, 121312, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 170872, 0, 3, 118750, 167512, 75682,
                                                 75832, 121522, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 171152, 0, 3, 118876, 167680, 75832,
                                                 75982, 121732, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 171432, 0, 3, 119002, 167848, 75982,
                                                 76132, 121942, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 171712, 0, 3, 119128, 168016, 76132,
                                                 76282, 122152, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 171992, 0, 3, 119254, 168184, 76282,
                                                 76432, 122362, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 172272, 0, 3, 119380, 168352, 76432,
                                                 76582, 122572, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 172552, 0, 3, 119632, 168520, 76882,
                                                 77032, 122782, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 172832, 0, 3, 119758, 168688, 77032,
                                                 77182, 122992, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 173112, 0, 3, 119884, 168856, 77182,
                                                 77332, 123202, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 173392, 0, 3, 120010, 169024, 77332,
                                                 77482, 123412, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 173672, 0, 3, 120136, 169192, 77482,
                                                 77632, 123622, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 173952, 0, 3, 120262, 169360, 77632,
                                                 77782, 123832, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 174232, 0, 3, 120388, 169528, 77782,
                                                 77932, 124042, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 174512, 0, 3, 120514, 169696, 77932,
                                                 78082, 124252, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 174792, 0, 3, 120640, 169864, 78082,
                                                 78232, 124462, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 175072, 0, 3, 167008, 167176, 121102,
                                                 170592, 78532, 78757, 124987, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 175492, 0, 3, 167176, 167344, 121312,
                                                 170872, 78757, 78982, 125302, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 175912, 0, 3, 167344, 167512, 121522,
                                                 171152, 78982, 79207, 125617, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 176332, 0, 3, 167512, 167680, 121732,
                                                 171432, 79207, 79432, 125932, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 176752, 0, 3, 167680, 167848, 121942,
                                                 171712, 79432, 79657, 126247, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 177172, 0, 3, 167848, 168016, 122152,
                                                 171992, 79657, 79882, 126562, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 177592, 0, 3, 168016, 168184, 122362,
                                                 172272, 79882, 80107, 126877, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 178012, 0, 3, 168520, 168688, 122992,
                                                 173112, 80557, 80782, 127507, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 178432, 0, 3, 168688, 168856, 123202,
                                                 173392, 80782, 81007, 127822, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 178852, 0, 3, 168856, 169024, 123412,
                                                 173672, 81007, 81232, 128137, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 179272, 0, 3, 169024, 169192, 123622,
                                                 173952, 81232, 81457, 128452, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 179692, 0, 3, 169192, 169360, 123832,
                                                 174232, 81457, 81682, 128767, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 180112, 0, 3, 169360, 169528, 124042,
                                                 174512, 81682, 81907, 129082, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 180532, 0, 3, 169528, 169696, 124252,
                                                 174792, 81907, 82132, 129397, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 180952, 0, 3, 170032, 170312, 124672,
                                                 175072, 82582, 82897, 129712, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 181540, 0, 3, 170312, 170592, 124987,
                                                 175492, 82897, 83212, 130153, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 182128, 0, 3, 170592, 170872, 125302,
                                                 175912, 83212, 83527, 130594, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 182716, 0, 3, 170872, 171152, 125617,
                                                 176332, 83527, 83842, 131035, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 183304, 0, 3, 171152, 171432, 125932,
                                                 176752, 83842, 84157, 131476, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 183892, 0, 3, 171432, 171712, 126247,
                                                 177172, 84157, 84472, 131917, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 184480, 0, 3, 171712, 171992, 126562,
                                                 177592, 84472, 84787, 132358, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 185068, 0, 3, 172552, 172832, 127192,
                                                 178012, 85417, 85732, 132799, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 185656, 0, 3, 172832, 173112, 127507,
                                                 178432, 85732, 86047, 133240, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 186244, 0, 3, 173112, 173392, 127822,
                                                 178852, 86047, 86362, 133681, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 186832, 0, 3, 173392, 173672, 128137,
                                                 179272, 86362, 86677, 134122, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 187420, 0, 3, 173672, 173952, 128452,
                                                 179692, 86677, 86992, 134563, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 188008, 0, 3, 173952, 174232, 128767,
                                                 180112, 86992, 87307, 135004, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 188596, 0, 3, 174232, 174512, 129082,
                                                 180532, 87307, 87622, 135445, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 189184, 0, 3, 175072, 175492, 130153,
                                                 182128, 88252, 88672, 136474, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 189968, 0, 3, 175492, 175912, 130594,
                                                 182716, 88672, 89092, 137062, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 190752, 0, 3, 175912, 176332, 131035,
                                                 183304, 89092, 89512, 137650, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 191536, 0, 3, 176332, 176752, 131476,
                                                 183892, 89512, 89932, 138238, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 192320, 0, 3, 176752, 177172, 131917,
                                                 184480, 89932, 90352, 138826, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 193104, 0, 3, 178012, 178432, 133240,
                                                 186244, 91192, 91612, 140002, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 193888, 0, 3, 178432, 178852, 133681,
                                                 186832, 91612, 92032, 140590, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 194672, 0, 3, 178852, 179272, 134122,
                                                 187420, 92032, 92452, 141178, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 195456, 0, 3, 179272, 179692, 134563,
                                                 188008, 92452, 92872, 141766, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 196240, 0, 3, 179692, 180112, 135004,
                                                 188596, 92872, 93292, 142354, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 197024, 0, 3, 180952, 181540, 135886,
                                                 189184, 94132, 94672, 142942, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 198032, 0, 3, 181540, 182128, 136474,
                                                 189968, 94672, 95212, 143698, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 199040, 0, 3, 182128, 182716, 137062,
                                                 190752, 95212, 95752, 144454, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 200048, 0, 3, 182716, 183304, 137650,
                                                 191536, 95752, 96292, 145210, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 201056, 0, 3, 183304, 183892, 138238,
                                                 192320, 96292, 96832, 145966, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 202064, 0, 3, 185068, 185656, 139414,
                                                 193104, 97912, 98452, 146722, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 203072, 0, 3, 185656, 186244, 140002,
                                                 193888, 98452, 98992, 147478, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 204080, 0, 3, 186244, 186832, 140590,
                                                 194672, 98992, 99532, 148234, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 205088, 0, 3, 186832, 187420, 141178,
                                                 195456, 99532, 100072, 148990, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 206096, 0, 3, 187420, 188008, 141766,
                                                 196240, 100072, 100612, 149746, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 207104, 0, 3, 189184, 189968, 143698,
                                                 199040, 101692, 102367, 151447, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 208364, 0, 3, 189968, 190752, 144454,
                                                 200048, 102367, 103042, 152392, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 209624, 0, 3, 190752, 191536, 145210,
                                                 201056, 103042, 103717, 153337, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 210884, 0, 3, 193104, 193888, 147478,
                                                 204080, 105067, 105742, 155227, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 212144, 0, 3, 193888, 194672, 148234,
                                                 205088, 105742, 106417, 156172, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 213404, 0, 3, 194672, 195456, 148990,
                                                 206096, 106417, 107092, 157117, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 214664, 0, 3, 197024, 198032, 150502,
                                                 207104, 108442, 109267, 158062, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 216204, 0, 3, 198032, 199040, 151447,
                                                 208364, 109267, 110092, 159217, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 217744, 0, 3, 199040, 200048, 152392,
                                                 209624, 110092, 110917, 160372, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 219284, 0, 3, 202064, 203072, 154282,
                                                 210884, 112567, 113392, 161527, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 220824, 0, 3, 203072, 204080, 155227,
                                                 212144, 113392, 114217, 162682, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 222364, 0, 3, 204080, 205088, 156172,
                                                 213404, 114217, 115042, 163837, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 223904, 3, 116692, 116713, 165020,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 223940, 3, 116713, 116734, 165048,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 223976, 3, 116734, 116755, 165076,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224012, 3, 116755, 116776, 165104,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224048, 3, 116776, 116797, 165132,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224084, 3, 116797, 116818, 165160,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224120, 3, 116818, 116839, 165188,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224156, 3, 116839, 116860, 165216,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224192, 3, 116902, 116923, 165272,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224228, 3, 116923, 116944, 165300,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224264, 3, 116944, 116965, 165328,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224300, 3, 116965, 116986, 165356,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224336, 3, 116986, 117007, 165384,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224372, 3, 117007, 117028, 165412,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224408, 3, 117028, 117049, 165440,
                                                 ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 224444, 3, 117049, 117070, 165468,
                                                 ncols, alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 224480, 0, 3, 164992, 223904, 165580,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 224588, 0, 3, 165020, 223940, 165664,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 224696, 0, 3, 165048, 223976, 165748,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 224804, 0, 3, 165076, 224012, 165832,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 224912, 0, 3, 165104, 224048, 165916,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225020, 0, 3, 165132, 224084, 166000,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225128, 0, 3, 165160, 224120, 166084,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225236, 0, 3, 165188, 224156, 166168,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225344, 0, 3, 165244, 224192, 166336,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225452, 0, 3, 165272, 224228, 166420,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225560, 0, 3, 165300, 224264, 166504,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225668, 0, 3, 165328, 224300, 166588,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225776, 0, 3, 165356, 224336, 166672,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225884, 0, 3, 165384, 224372, 166756,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 225992, 0, 3, 165412, 224408, 166840,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 226100, 0, 3, 165440, 224444, 166924,
                                                 ncols, p);

            compute_prim_dk_electron_repulsion_0(buffer, 226208, 0, 3, 165496, 224480, 118372,
                                                 118498, 167176, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 226424, 0, 3, 165580, 224588, 118498,
                                                 118624, 167344, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 226640, 0, 3, 165664, 224696, 118624,
                                                 118750, 167512, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 226856, 0, 3, 165748, 224804, 118750,
                                                 118876, 167680, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 227072, 0, 3, 165832, 224912, 118876,
                                                 119002, 167848, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 227288, 0, 3, 165916, 225020, 119002,
                                                 119128, 168016, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 227504, 0, 3, 166000, 225128, 119128,
                                                 119254, 168184, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 227720, 0, 3, 166084, 225236, 119254,
                                                 119380, 168352, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 227936, 0, 3, 166252, 225344, 119632,
                                                 119758, 168688, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 228152, 0, 3, 166336, 225452, 119758,
                                                 119884, 168856, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 228368, 0, 3, 166420, 225560, 119884,
                                                 120010, 169024, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 228584, 0, 3, 166504, 225668, 120010,
                                                 120136, 169192, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 228800, 0, 3, 166588, 225776, 120136,
                                                 120262, 169360, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 229016, 0, 3, 166672, 225884, 120262,
                                                 120388, 169528, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 229232, 0, 3, 166756, 225992, 120388,
                                                 120514, 169696, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 229448, 0, 3, 166840, 226100, 120514,
                                                 120640, 169864, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 229664, 0, 3, 167176, 226424, 120892,
                                                 121102, 170592, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 230024, 0, 3, 167344, 226640, 121102,
                                                 121312, 170872, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 230384, 0, 3, 167512, 226856, 121312,
                                                 121522, 171152, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 230744, 0, 3, 167680, 227072, 121522,
                                                 121732, 171432, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 231104, 0, 3, 167848, 227288, 121732,
                                                 121942, 171712, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 231464, 0, 3, 168016, 227504, 121942,
                                                 122152, 171992, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 231824, 0, 3, 168184, 227720, 122152,
                                                 122362, 172272, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 232184, 0, 3, 168688, 228152, 122782,
                                                 122992, 173112, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 232544, 0, 3, 168856, 228368, 122992,
                                                 123202, 173392, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 232904, 0, 3, 169024, 228584, 123202,
                                                 123412, 173672, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 233264, 0, 3, 169192, 228800, 123412,
                                                 123622, 173952, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 233624, 0, 3, 169360, 229016, 123622,
                                                 123832, 174232, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 233984, 0, 3, 169528, 229232, 123832,
                                                 124042, 174512, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 234344, 0, 3, 169696, 229448, 124042,
                                                 124252, 174792, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 234704, 0, 3, 226208, 226424, 170592,
                                                 230024, 124672, 124987, 175492, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 235244, 0, 3, 226424, 226640, 170872,
                                                 230384, 124987, 125302, 175912, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 235784, 0, 3, 226640, 226856, 171152,
                                                 230744, 125302, 125617, 176332, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 236324, 0, 3, 226856, 227072, 171432,
                                                 231104, 125617, 125932, 176752, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 236864, 0, 3, 227072, 227288, 171712,
                                                 231464, 125932, 126247, 177172, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 237404, 0, 3, 227288, 227504, 171992,
                                                 231824, 126247, 126562, 177592, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 237944, 0, 3, 227936, 228152, 173112,
                                                 232544, 127192, 127507, 178432, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 238484, 0, 3, 228152, 228368, 173392,
                                                 232904, 127507, 127822, 178852, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 239024, 0, 3, 228368, 228584, 173672,
                                                 233264, 127822, 128137, 179272, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 239564, 0, 3, 228584, 228800, 173952,
                                                 233624, 128137, 128452, 179692, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 240104, 0, 3, 228800, 229016, 174232,
                                                 233984, 128452, 128767, 180112, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 240644, 0, 3, 229016, 229232, 174512,
                                                 234344, 128767, 129082, 180532, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 241184, 0, 3, 229664, 230024, 175492,
                                                 235244, 129712, 130153, 182128, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 241940, 0, 3, 230024, 230384, 175912,
                                                 235784, 130153, 130594, 182716, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 242696, 0, 3, 230384, 230744, 176332,
                                                 236324, 130594, 131035, 183304, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 243452, 0, 3, 230744, 231104, 176752,
                                                 236864, 131035, 131476, 183892, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 244208, 0, 3, 231104, 231464, 177172,
                                                 237404, 131476, 131917, 184480, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 244964, 0, 3, 232184, 232544, 178432,
                                                 238484, 132799, 133240, 186244, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 245720, 0, 3, 232544, 232904, 178852,
                                                 239024, 133240, 133681, 186832, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 246476, 0, 3, 232904, 233264, 179272,
                                                 239564, 133681, 134122, 187420, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 247232, 0, 3, 233264, 233624, 179692,
                                                 240104, 134122, 134563, 188008, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 247988, 0, 3, 233624, 233984, 180112,
                                                 240644, 134563, 135004, 188596, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 248744, 0, 3, 234704, 235244, 182128,
                                                 241940, 135886, 136474, 189968, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 249752, 0, 3, 235244, 235784, 182716,
                                                 242696, 136474, 137062, 190752, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 250760, 0, 3, 235784, 236324, 183304,
                                                 243452, 137062, 137650, 191536, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 251768, 0, 3, 236324, 236864, 183892,
                                                 244208, 137650, 138238, 192320, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 252776, 0, 3, 237944, 238484, 186244,
                                                 245720, 139414, 140002, 193888, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 253784, 0, 3, 238484, 239024, 186832,
                                                 246476, 140002, 140590, 194672, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 254792, 0, 3, 239024, 239564, 187420,
                                                 247232, 140590, 141178, 195456, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 255800, 0, 3, 239564, 240104, 188008,
                                                 247988, 141178, 141766, 196240, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 256808, 0, 3, 241184, 241940, 189968,
                                                 249752, 142942, 143698, 199040, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 258104, 0, 3, 241940, 242696, 190752,
                                                 250760, 143698, 144454, 200048, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 259400, 0, 3, 242696, 243452, 191536,
                                                 251768, 144454, 145210, 201056, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 260696, 0, 3, 244964, 245720, 193888,
                                                 253784, 146722, 147478, 204080, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 261992, 0, 3, 245720, 246476, 194672,
                                                 254792, 147478, 148234, 205088, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 263288, 0, 3, 246476, 247232, 195456,
                                                 255800, 148234, 148990, 206096, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 264584, 0, 3, 248744, 249752, 199040,
                                                 258104, 150502, 151447, 208364, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 266204, 0, 3, 249752, 250760, 200048,
                                                 259400, 151447, 152392, 209624, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 267824, 0, 3, 252776, 253784, 204080,
                                                 261992, 154282, 155227, 212144, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 269444, 0, 3, 253784, 254792, 205088,
                                                 263288, 155227, 156172, 213404, ncols, alpha,
                                                 beta, p);

            compute_prim_mk_electron_repulsion_0(buffer, 271064, 0, 3, 256808, 258104, 208364,
                                                 266204, 158062, 159217, 217744, ncols, alpha,
                                                 beta, p);

            compute_prim_mk_electron_repulsion_0(buffer, 273044, 0, 3, 260696, 261992, 212144,
                                                 269444, 161527, 162682, 222364, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275024, 3, 164992, 165020, 223940,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275069, 3, 165020, 165048, 223976,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275114, 3, 165048, 165076, 224012,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275159, 3, 165076, 165104, 224048,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275204, 3, 165104, 165132, 224084,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275249, 3, 165132, 165160, 224120,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275294, 3, 165160, 165188, 224156,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275339, 3, 165244, 165272, 224228,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275384, 3, 165272, 165300, 224264,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275429, 3, 165300, 165328, 224300,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275474, 3, 165328, 165356, 224336,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275519, 3, 165356, 165384, 224372,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275564, 3, 165384, 165412, 224408,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 275609, 3, 165412, 165440, 224444,
                                                 ncols, alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 275654, 0, 3, 223904, 275024, 224588,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 275789, 0, 3, 223940, 275069, 224696,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 275924, 0, 3, 223976, 275114, 224804,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 276059, 0, 3, 224012, 275159, 224912,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 276194, 0, 3, 224048, 275204, 225020,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 276329, 0, 3, 224084, 275249, 225128,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 276464, 0, 3, 224120, 275294, 225236,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 276599, 0, 3, 224192, 275339, 225452,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 276734, 0, 3, 224228, 275384, 225560,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 276869, 0, 3, 224264, 275429, 225668,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 277004, 0, 3, 224300, 275474, 225776,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 277139, 0, 3, 224336, 275519, 225884,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 277274, 0, 3, 224372, 275564, 225992,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 277409, 0, 3, 224408, 275609, 226100,
                                                 ncols, p);

            compute_prim_dl_electron_repulsion_0(buffer, 277544, 0, 3, 224480, 275654, 167008,
                                                 167176, 226424, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 277814, 0, 3, 224588, 275789, 167176,
                                                 167344, 226640, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 278084, 0, 3, 224696, 275924, 167344,
                                                 167512, 226856, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 278354, 0, 3, 224804, 276059, 167512,
                                                 167680, 227072, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 278624, 0, 3, 224912, 276194, 167680,
                                                 167848, 227288, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 278894, 0, 3, 225020, 276329, 167848,
                                                 168016, 227504, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 279164, 0, 3, 225128, 276464, 168016,
                                                 168184, 227720, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 279434, 0, 3, 225344, 276599, 168520,
                                                 168688, 228152, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 279704, 0, 3, 225452, 276734, 168688,
                                                 168856, 228368, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 279974, 0, 3, 225560, 276869, 168856,
                                                 169024, 228584, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 280244, 0, 3, 225668, 277004, 169024,
                                                 169192, 228800, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 280514, 0, 3, 225776, 277139, 169192,
                                                 169360, 229016, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 280784, 0, 3, 225884, 277274, 169360,
                                                 169528, 229232, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 281054, 0, 3, 225992, 277409, 169528,
                                                 169696, 229448, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 281324, 0, 3, 226208, 277544, 170032,
                                                 170312, 229664, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 281774, 0, 3, 226424, 277814, 170312,
                                                 170592, 230024, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 282224, 0, 3, 226640, 278084, 170592,
                                                 170872, 230384, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 282674, 0, 3, 226856, 278354, 170872,
                                                 171152, 230744, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 283124, 0, 3, 227072, 278624, 171152,
                                                 171432, 231104, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 283574, 0, 3, 227288, 278894, 171432,
                                                 171712, 231464, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 284024, 0, 3, 227504, 279164, 171712,
                                                 171992, 231824, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 284474, 0, 3, 227936, 279434, 172552,
                                                 172832, 232184, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 284924, 0, 3, 228152, 279704, 172832,
                                                 173112, 232544, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 285374, 0, 3, 228368, 279974, 173112,
                                                 173392, 232904, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 285824, 0, 3, 228584, 280244, 173392,
                                                 173672, 233264, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 286274, 0, 3, 228800, 280514, 173672,
                                                 173952, 233624, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 286724, 0, 3, 229016, 280784, 173952,
                                                 174232, 233984, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 287174, 0, 3, 229232, 281054, 174232,
                                                 174512, 234344, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 287624, 0, 3, 277544, 277814, 230024,
                                                 282224, 175072, 175492, 235244, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 288299, 0, 3, 277814, 278084, 230384,
                                                 282674, 175492, 175912, 235784, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 288974, 0, 3, 278084, 278354, 230744,
                                                 283124, 175912, 176332, 236324, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 289649, 0, 3, 278354, 278624, 231104,
                                                 283574, 176332, 176752, 236864, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 290324, 0, 3, 278624, 278894, 231464,
                                                 284024, 176752, 177172, 237404, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 290999, 0, 3, 279434, 279704, 232544,
                                                 285374, 178012, 178432, 238484, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 291674, 0, 3, 279704, 279974, 232904,
                                                 285824, 178432, 178852, 239024, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 292349, 0, 3, 279974, 280244, 233264,
                                                 286274, 178852, 179272, 239564, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 293024, 0, 3, 280244, 280514, 233624,
                                                 286724, 179272, 179692, 240104, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 293699, 0, 3, 280514, 280784, 233984,
                                                 287174, 179692, 180112, 240644, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 294374, 0, 3, 281324, 281774, 234704,
                                                 287624, 180952, 181540, 241184, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 295319, 0, 3, 281774, 282224, 235244,
                                                 288299, 181540, 182128, 241940, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 296264, 0, 3, 282224, 282674, 235784,
                                                 288974, 182128, 182716, 242696, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 297209, 0, 3, 282674, 283124, 236324,
                                                 289649, 182716, 183304, 243452, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 298154, 0, 3, 283124, 283574, 236864,
                                                 290324, 183304, 183892, 244208, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 299099, 0, 3, 284474, 284924, 237944,
                                                 290999, 185068, 185656, 244964, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 300044, 0, 3, 284924, 285374, 238484,
                                                 291674, 185656, 186244, 245720, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 300989, 0, 3, 285374, 285824, 239024,
                                                 292349, 186244, 186832, 246476, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 301934, 0, 3, 285824, 286274, 239564,
                                                 293024, 186832, 187420, 247232, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 302879, 0, 3, 286274, 286724, 240104,
                                                 293699, 187420, 188008, 247988, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 303824, 0, 3, 287624, 288299, 241940,
                                                 296264, 189184, 189968, 249752, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 305084, 0, 3, 288299, 288974, 242696,
                                                 297209, 189968, 190752, 250760, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 306344, 0, 3, 288974, 289649, 243452,
                                                 298154, 190752, 191536, 251768, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 307604, 0, 3, 290999, 291674, 245720,
                                                 300989, 193104, 193888, 253784, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 308864, 0, 3, 291674, 292349, 246476,
                                                 301934, 193888, 194672, 254792, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 310124, 0, 3, 292349, 293024, 247232,
                                                 302879, 194672, 195456, 255800, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 311384, 0, 3, 294374, 295319, 248744,
                                                 303824, 197024, 198032, 256808, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 313004, 0, 3, 295319, 296264, 249752,
                                                 305084, 198032, 199040, 258104, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 314624, 0, 3, 296264, 297209, 250760,
                                                 306344, 199040, 200048, 259400, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 316244, 0, 3, 299099, 300044, 252776,
                                                 307604, 202064, 203072, 260696, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 317864, 0, 3, 300044, 300989, 253784,
                                                 308864, 203072, 204080, 261992, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 319484, 0, 3, 300989, 301934, 254792,
                                                 310124, 204080, 205088, 263288, ncols, alpha,
                                                 beta, p);

            compute_prim_ll_electron_repulsion_0(buffer, 321104, 0, 3, 303824, 305084, 258104,
                                                 314624, 207104, 208364, 266204, ncols, alpha,
                                                 beta, p);

            compute_prim_ll_electron_repulsion_0(buffer, 323129, 0, 3, 307604, 308864, 261992,
                                                 319484, 210884, 212144, 269444, ncols, alpha,
                                                 beta, p);

            compute_prim_ml_electron_repulsion_0(buffer, 325154, 0, 3, 311384, 313004, 264584,
                                                 321104, 214664, 216204, 271064, ncols, alpha,
                                                 beta, p);

            compute_prim_ml_electron_repulsion_0(buffer, 327629, 0, 3, 316244, 317864, 267824,
                                                 323129, 219284, 220824, 273044, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_l_x(buffer, 330104, 316244, 327629, 1, 45, ncols, alpha);

            simdgeo::geom_l_y(buffer, 332129, 316244, 327629, 1, 45, ncols, alpha);

            simdgeo::geom_l_z(buffer, 334154, 316244, 327629, 1, 45, ncols, alpha);

            simdgeo::geom_l_x(buffer, 336179, 311384, 325154, 1, 45, ncols, alpha);

            simdgeo::geom_l_y(buffer, 338204, 311384, 325154, 1, 45, ncols, alpha);

            simdgeo::geom_l_z(buffer, 340229, 311384, 325154, 1, 45, ncols, alpha);

            simdfunc::contract_primitives(buffer, 342254, 336179, 6075, ncols);

            simdfunc::contract_primitives(buffer, 348329, 330104, 6075, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 354404, 348329, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 354404, 17, nmax);

    simdtrf::transform_l_inner(buffer, 354404, 350354, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 289 * nvalues, nvalues, buffer, 354404, 17, nmax);

    simdtrf::transform_l_inner(buffer, 354404, 352379, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 578 * nvalues, nvalues, buffer, 354404, 17, nmax);

    simdtrf::transform_l_inner(buffer, 354404, 342254, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 867 * nvalues, nvalues, buffer, 354404, 17, nmax);

    simdtrf::transform_l_inner(buffer, 354404, 344279, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 1156 * nvalues, nvalues, buffer, 354404, 17, nmax);

    simdtrf::transform_l_inner(buffer, 354404, 346304, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 1445 * nvalues, nvalues, buffer, 354404, 17, nmax);
}

}  // namespace simdt2ceri
