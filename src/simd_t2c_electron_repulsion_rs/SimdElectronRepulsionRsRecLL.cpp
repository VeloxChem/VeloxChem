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


#include "SimdElectronRepulsionRsRecLL.hpp"

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
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_ll_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ll_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 240157, 235342, 4050, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 16, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 24, 16, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 42, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 45, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 48, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 51, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 54, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 57, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 60, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 63, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 66, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 69, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 72, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 75, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 78, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 81, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 84, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 87, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 90, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 93, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 96, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 99, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 102, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 105, 0, 33, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 108, 0, 34, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 111, 0, 35, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 114, 0, 36, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 117, 0, 37, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 120, 0, 38, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 123, 0, 39, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 126, 0, 40, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 129, 0, 41, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 132, 0, 7, 8, 42, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 138, 0, 8, 9, 45, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 144, 0, 9, 10, 48, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 150, 0, 10, 11, 51, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 156, 0, 11, 12, 54, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 162, 0, 12, 13, 57, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 168, 0, 13, 14, 60, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 174, 0, 14, 15, 63, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 180, 0, 15, 16, 66, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 186, 0, 16, 17, 69, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 192, 0, 17, 18, 72, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 198, 0, 18, 19, 75, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 204, 0, 19, 20, 78, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 210, 0, 20, 21, 81, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 216, 0, 21, 22, 84, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 222, 0, 25, 26, 87, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 228, 0, 26, 27, 90, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 234, 0, 27, 28, 93, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 240, 0, 28, 29, 96, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 246, 0, 29, 30, 99, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 252, 0, 30, 31, 102, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 258, 0, 31, 32, 105, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 264, 0, 32, 33, 108, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 270, 0, 33, 34, 111, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 276, 0, 34, 35, 114, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 282, 0, 35, 36, 117, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 288, 0, 36, 37, 120, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 294, 0, 37, 38, 123, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 300, 0, 38, 39, 126, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 306, 0, 39, 40, 129, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 312, 0, 42, 45, 144, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 322, 0, 45, 48, 150, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 332, 0, 48, 51, 156, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 342, 0, 51, 54, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 352, 0, 54, 57, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 362, 0, 57, 60, 174, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 372, 0, 60, 63, 180, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 382, 0, 63, 66, 186, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 392, 0, 66, 69, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 402, 0, 69, 72, 198, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 412, 0, 72, 75, 204, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 422, 0, 75, 78, 210, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 432, 0, 78, 81, 216, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 442, 0, 87, 90, 234, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 452, 0, 90, 93, 240, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 462, 0, 93, 96, 246, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 472, 0, 96, 99, 252, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 482, 0, 99, 102, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 492, 0, 102, 105, 264, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 502, 0, 105, 108, 270, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 512, 0, 108, 111, 276, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 522, 0, 111, 114, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 532, 0, 114, 117, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 542, 0, 117, 120, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 552, 0, 120, 123, 300, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 562, 0, 123, 126, 306, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 572, 0, 132, 138, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 587, 0, 138, 144, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 602, 0, 144, 150, 332, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 617, 0, 150, 156, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 632, 0, 156, 162, 352, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 647, 0, 162, 168, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 662, 0, 168, 174, 372, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 677, 0, 174, 180, 382, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 692, 0, 180, 186, 392, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 707, 0, 186, 192, 402, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 722, 0, 192, 198, 412, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 737, 0, 198, 204, 422, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 752, 0, 204, 210, 432, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 767, 0, 222, 228, 442, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 782, 0, 228, 234, 452, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 797, 0, 234, 240, 462, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 812, 0, 240, 246, 472, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 827, 0, 246, 252, 482, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 842, 0, 252, 258, 492, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 857, 0, 258, 264, 502, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 872, 0, 264, 270, 512, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 887, 0, 270, 276, 522, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 902, 0, 276, 282, 532, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 917, 0, 282, 288, 542, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 932, 0, 288, 294, 552, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 947, 0, 294, 300, 562, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 962, 0, 312, 322, 602, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 983, 0, 322, 332, 617, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1004, 0, 332, 342, 632, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1025, 0, 342, 352, 647, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1046, 0, 352, 362, 662, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1067, 0, 362, 372, 677, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1088, 0, 372, 382, 692, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1109, 0, 382, 392, 707, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1130, 0, 392, 402, 722, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1151, 0, 402, 412, 737, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1172, 0, 412, 422, 752, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1193, 0, 442, 452, 797, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1214, 0, 452, 462, 812, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1235, 0, 462, 472, 827, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1256, 0, 472, 482, 842, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1277, 0, 482, 492, 857, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1298, 0, 492, 502, 872, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1319, 0, 502, 512, 887, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1340, 0, 512, 522, 902, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1361, 0, 522, 532, 917, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1382, 0, 532, 542, 932, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1403, 0, 542, 552, 947, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1424, 0, 572, 587, 962, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1452, 0, 587, 602, 983, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1480, 0, 602, 617, 1004, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1508, 0, 617, 632, 1025, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1536, 0, 632, 647, 1046, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1564, 0, 647, 662, 1067, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1592, 0, 662, 677, 1088, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1620, 0, 677, 692, 1109, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1648, 0, 692, 707, 1130, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1676, 0, 707, 722, 1151, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1704, 0, 722, 737, 1172, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1732, 0, 767, 782, 1193, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1760, 0, 782, 797, 1214, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1788, 0, 797, 812, 1235, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1816, 0, 812, 827, 1256, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1844, 0, 827, 842, 1277, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1872, 0, 842, 857, 1298, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1900, 0, 857, 872, 1319, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1928, 0, 872, 887, 1340, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1956, 0, 887, 902, 1361, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1984, 0, 902, 917, 1382, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 2012, 0, 917, 932, 1403, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2040, 0, 962, 983, 1480, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2076, 0, 983, 1004, 1508, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2112, 0, 1004, 1025, 1536, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2148, 0, 1025, 1046, 1564, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2184, 0, 1046, 1067, 1592, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2220, 0, 1067, 1088, 1620, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2256, 0, 1088, 1109, 1648, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2292, 0, 1109, 1130, 1676, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2328, 0, 1130, 1151, 1704, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2364, 0, 1193, 1214, 1788, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2400, 0, 1214, 1235, 1816, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2436, 0, 1235, 1256, 1844, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2472, 0, 1256, 1277, 1872, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2508, 0, 1277, 1298, 1900, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2544, 0, 1298, 1319, 1928, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2580, 0, 1319, 1340, 1956, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2616, 0, 1340, 1361, 1984, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2652, 0, 1361, 1382, 2012, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2688, 0, 1424, 1452, 2040, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2733, 0, 1452, 1480, 2076, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2778, 0, 1480, 1508, 2112, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2823, 0, 1508, 1536, 2148, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2868, 0, 1536, 1564, 2184, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2913, 0, 1564, 1592, 2220, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2958, 0, 1592, 1620, 2256, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3003, 0, 1620, 1648, 2292, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3048, 0, 1648, 1676, 2328, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3093, 0, 1732, 1760, 2364, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3138, 0, 1760, 1788, 2400, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3183, 0, 1788, 1816, 2436, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3228, 0, 1816, 1844, 2472, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3273, 0, 1844, 1872, 2508, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3318, 0, 1872, 1900, 2544, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3363, 0, 1900, 1928, 2580, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3408, 0, 1928, 1956, 2616, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3453, 0, 1956, 1984, 2652, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 3498, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3501, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3504, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3507, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3510, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3513, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3516, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3519, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3522, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3525, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3528, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3531, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3534, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3537, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3540, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3543, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3546, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3549, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3552, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3555, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3558, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3561, 3, 35, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3564, 3, 36, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3567, 3, 37, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3570, 3, 38, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3573, 3, 39, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3576, 3, 40, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3579, 3, 41, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 3582, 3, 9, 45, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3591, 3, 10, 48, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3600, 3, 11, 51, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3609, 3, 12, 54, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3618, 3, 13, 57, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3627, 3, 14, 60, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3636, 3, 15, 63, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3645, 3, 16, 66, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3654, 3, 17, 69, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3663, 3, 18, 72, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3672, 3, 19, 75, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3681, 3, 20, 78, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3690, 3, 21, 81, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3699, 3, 22, 84, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3708, 3, 27, 90, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3717, 3, 28, 93, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3726, 3, 29, 96, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3735, 3, 30, 99, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3744, 3, 31, 102, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3753, 3, 32, 105, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3762, 3, 33, 108, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3771, 3, 34, 111, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3780, 3, 35, 114, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3789, 3, 36, 117, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3798, 3, 37, 120, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3807, 3, 38, 123, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3816, 3, 39, 126, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3825, 3, 40, 129, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3834, 0, 3, 45, 3591, 144, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3852, 0, 3, 48, 3600, 150, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3870, 0, 3, 51, 3609, 156, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3888, 0, 3, 54, 3618, 162, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3906, 0, 3, 57, 3627, 168, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3924, 0, 3, 60, 3636, 174, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3942, 0, 3, 63, 3645, 180, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3960, 0, 3, 66, 3654, 186, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3978, 0, 3, 69, 3663, 192, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3996, 0, 3, 72, 3672, 198, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4014, 0, 3, 75, 3681, 204, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4032, 0, 3, 78, 3690, 210, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4050, 0, 3, 81, 3699, 216, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4068, 0, 3, 90, 3717, 234, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4086, 0, 3, 93, 3726, 240, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4104, 0, 3, 96, 3735, 246, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4122, 0, 3, 99, 3744, 252, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4140, 0, 3, 102, 3753, 258, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4158, 0, 3, 105, 3762, 264, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4176, 0, 3, 108, 3771, 270, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4194, 0, 3, 111, 3780, 276, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4212, 0, 3, 114, 3789, 282, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4230, 0, 3, 117, 3798, 288, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4248, 0, 3, 120, 3807, 294, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4266, 0, 3, 123, 3816, 300, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4284, 0, 3, 126, 3825, 306, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4302, 0, 3, 144, 3852, 322, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4332, 0, 3, 150, 3870, 332, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4362, 0, 3, 156, 3888, 342, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4392, 0, 3, 162, 3906, 352, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4422, 0, 3, 168, 3924, 362, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4452, 0, 3, 174, 3942, 372, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4482, 0, 3, 180, 3960, 382, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4512, 0, 3, 186, 3978, 392, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4542, 0, 3, 192, 3996, 402, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4572, 0, 3, 198, 4014, 412, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4602, 0, 3, 204, 4032, 422, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4632, 0, 3, 210, 4050, 432, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4662, 0, 3, 234, 4086, 452, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4692, 0, 3, 240, 4104, 462, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4722, 0, 3, 246, 4122, 472, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4752, 0, 3, 252, 4140, 482, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4782, 0, 3, 258, 4158, 492, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4812, 0, 3, 264, 4176, 502, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4842, 0, 3, 270, 4194, 512, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4872, 0, 3, 276, 4212, 522, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4902, 0, 3, 282, 4230, 532, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4932, 0, 3, 288, 4248, 542, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4962, 0, 3, 294, 4266, 552, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4992, 0, 3, 300, 4284, 562, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5022, 0, 3, 322, 4332, 602, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5067, 0, 3, 332, 4362, 617, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5112, 0, 3, 342, 4392, 632, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5157, 0, 3, 352, 4422, 647, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5202, 0, 3, 362, 4452, 662, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5247, 0, 3, 372, 4482, 677, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5292, 0, 3, 382, 4512, 692, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5337, 0, 3, 392, 4542, 707, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5382, 0, 3, 402, 4572, 722, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5427, 0, 3, 412, 4602, 737, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5472, 0, 3, 422, 4632, 752, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5517, 0, 3, 452, 4692, 797, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5562, 0, 3, 462, 4722, 812, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5607, 0, 3, 472, 4752, 827, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5652, 0, 3, 482, 4782, 842, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5697, 0, 3, 492, 4812, 857, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5742, 0, 3, 502, 4842, 872, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5787, 0, 3, 512, 4872, 887, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5832, 0, 3, 522, 4902, 902, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5877, 0, 3, 532, 4932, 917, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5922, 0, 3, 542, 4962, 932, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5967, 0, 3, 552, 4992, 947, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6012, 0, 3, 602, 5067, 983, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6075, 0, 3, 617, 5112, 1004, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6138, 0, 3, 632, 5157, 1025, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6201, 0, 3, 647, 5202, 1046, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6264, 0, 3, 662, 5247, 1067, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6327, 0, 3, 677, 5292, 1088, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6390, 0, 3, 692, 5337, 1109, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6453, 0, 3, 707, 5382, 1130, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6516, 0, 3, 722, 5427, 1151, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6579, 0, 3, 737, 5472, 1172, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6642, 0, 3, 797, 5562, 1214, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6705, 0, 3, 812, 5607, 1235, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6768, 0, 3, 827, 5652, 1256, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6831, 0, 3, 842, 5697, 1277, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6894, 0, 3, 857, 5742, 1298, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6957, 0, 3, 872, 5787, 1319, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7020, 0, 3, 887, 5832, 1340, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7083, 0, 3, 902, 5877, 1361, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7146, 0, 3, 917, 5922, 1382, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7209, 0, 3, 932, 5967, 1403, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7272, 0, 3, 983, 6075, 1480, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7356, 0, 3, 1004, 6138, 1508, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7440, 0, 3, 1025, 6201, 1536, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7524, 0, 3, 1046, 6264, 1564, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7608, 0, 3, 1067, 6327, 1592, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7692, 0, 3, 1088, 6390, 1620, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7776, 0, 3, 1109, 6453, 1648, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7860, 0, 3, 1130, 6516, 1676, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7944, 0, 3, 1151, 6579, 1704, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8028, 0, 3, 1214, 6705, 1788, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8112, 0, 3, 1235, 6768, 1816, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8196, 0, 3, 1256, 6831, 1844, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8280, 0, 3, 1277, 6894, 1872, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8364, 0, 3, 1298, 6957, 1900, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8448, 0, 3, 1319, 7020, 1928, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8532, 0, 3, 1340, 7083, 1956, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8616, 0, 3, 1361, 7146, 1984, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8700, 0, 3, 1382, 7209, 2012, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8784, 0, 3, 1480, 7356, 2076, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8892, 0, 3, 1508, 7440, 2112, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9000, 0, 3, 1536, 7524, 2148, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9108, 0, 3, 1564, 7608, 2184, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9216, 0, 3, 1592, 7692, 2220, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9324, 0, 3, 1620, 7776, 2256, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9432, 0, 3, 1648, 7860, 2292, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9540, 0, 3, 1676, 7944, 2328, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9648, 0, 3, 1788, 8112, 2400, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9756, 0, 3, 1816, 8196, 2436, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9864, 0, 3, 1844, 8280, 2472, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9972, 0, 3, 1872, 8364, 2508, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10080, 0, 3, 1900, 8448, 2544, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10188, 0, 3, 1928, 8532, 2580, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10296, 0, 3, 1956, 8616, 2616, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10404, 0, 3, 1984, 8700, 2652, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10512, 0, 3, 2076, 8892, 2778, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10647, 0, 3, 2112, 9000, 2823, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10782, 0, 3, 2148, 9108, 2868, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10917, 0, 3, 2184, 9216, 2913, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11052, 0, 3, 2220, 9324, 2958, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11187, 0, 3, 2256, 9432, 3003, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11322, 0, 3, 2292, 9540, 3048, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11457, 0, 3, 2400, 9756, 3183, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11592, 0, 3, 2436, 9864, 3228, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11727, 0, 3, 2472, 9972, 3273, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11862, 0, 3, 2508, 10080, 3318, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11997, 0, 3, 2544, 10188, 3363, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12132, 0, 3, 2580, 10296, 3408, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12267, 0, 3, 2616, 10404, 3453, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 12402, 3, 9, 10, 3501, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12408, 3, 10, 11, 3504, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12414, 3, 11, 12, 3507, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12420, 3, 12, 13, 3510, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12426, 3, 13, 14, 3513, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12432, 3, 14, 15, 3516, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12438, 3, 15, 16, 3519, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12444, 3, 16, 17, 3522, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12450, 3, 17, 18, 3525, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12456, 3, 18, 19, 3528, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12462, 3, 19, 20, 3531, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12468, 3, 20, 21, 3534, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12474, 3, 21, 22, 3537, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12480, 3, 27, 28, 3543, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12486, 3, 28, 29, 3546, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12492, 3, 29, 30, 3549, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12498, 3, 30, 31, 3552, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12504, 3, 31, 32, 3555, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12510, 3, 32, 33, 3558, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12516, 3, 33, 34, 3561, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12522, 3, 34, 35, 3564, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12528, 3, 35, 36, 3567, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12534, 3, 36, 37, 3570, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12540, 3, 37, 38, 3573, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12546, 3, 38, 39, 3576, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 12552, 3, 39, 40, 3579, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 12558, 0, 3, 3498, 12402, 3591, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12576, 0, 3, 3501, 12408, 3600, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12594, 0, 3, 3504, 12414, 3609, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12612, 0, 3, 3507, 12420, 3618, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12630, 0, 3, 3510, 12426, 3627, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12648, 0, 3, 3513, 12432, 3636, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12666, 0, 3, 3516, 12438, 3645, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12684, 0, 3, 3519, 12444, 3654, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12702, 0, 3, 3522, 12450, 3663, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12720, 0, 3, 3525, 12456, 3672, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12738, 0, 3, 3528, 12462, 3681, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12756, 0, 3, 3531, 12468, 3690, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12774, 0, 3, 3534, 12474, 3699, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12792, 0, 3, 3540, 12480, 3717, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12810, 0, 3, 3543, 12486, 3726, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12828, 0, 3, 3546, 12492, 3735, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12846, 0, 3, 3549, 12498, 3744, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12864, 0, 3, 3552, 12504, 3753, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12882, 0, 3, 3555, 12510, 3762, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12900, 0, 3, 3558, 12516, 3771, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12918, 0, 3, 3561, 12522, 3780, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12936, 0, 3, 3564, 12528, 3789, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12954, 0, 3, 3567, 12534, 3798, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12972, 0, 3, 3570, 12540, 3807, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 12990, 0, 3, 3573, 12546, 3816, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13008, 0, 3, 3576, 12552, 3825, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 13026, 0, 3, 3582, 12558, 132, 138,
                                                 3834, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13062, 0, 3, 3591, 12576, 138, 144,
                                                 3852, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13098, 0, 3, 3600, 12594, 144, 150,
                                                 3870, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13134, 0, 3, 3609, 12612, 150, 156,
                                                 3888, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13170, 0, 3, 3618, 12630, 156, 162,
                                                 3906, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13206, 0, 3, 3627, 12648, 162, 168,
                                                 3924, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13242, 0, 3, 3636, 12666, 168, 174,
                                                 3942, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13278, 0, 3, 3645, 12684, 174, 180,
                                                 3960, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13314, 0, 3, 3654, 12702, 180, 186,
                                                 3978, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13350, 0, 3, 3663, 12720, 186, 192,
                                                 3996, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13386, 0, 3, 3672, 12738, 192, 198,
                                                 4014, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13422, 0, 3, 3681, 12756, 198, 204,
                                                 4032, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13458, 0, 3, 3690, 12774, 204, 210,
                                                 4050, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13494, 0, 3, 3708, 12792, 222, 228,
                                                 4068, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13530, 0, 3, 3717, 12810, 228, 234,
                                                 4086, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13566, 0, 3, 3726, 12828, 234, 240,
                                                 4104, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13602, 0, 3, 3735, 12846, 240, 246,
                                                 4122, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13638, 0, 3, 3744, 12864, 246, 252,
                                                 4140, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13674, 0, 3, 3753, 12882, 252, 258,
                                                 4158, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13710, 0, 3, 3762, 12900, 258, 264,
                                                 4176, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13746, 0, 3, 3771, 12918, 264, 270,
                                                 4194, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13782, 0, 3, 3780, 12936, 270, 276,
                                                 4212, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13818, 0, 3, 3789, 12954, 276, 282,
                                                 4230, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13854, 0, 3, 3798, 12972, 282, 288,
                                                 4248, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13890, 0, 3, 3807, 12990, 288, 294,
                                                 4266, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13926, 0, 3, 3816, 13008, 294, 300,
                                                 4284, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 13962, 0, 3, 3852, 13098, 312, 322,
                                                 4332, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14022, 0, 3, 3870, 13134, 322, 332,
                                                 4362, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14082, 0, 3, 3888, 13170, 332, 342,
                                                 4392, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14142, 0, 3, 3906, 13206, 342, 352,
                                                 4422, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14202, 0, 3, 3924, 13242, 352, 362,
                                                 4452, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14262, 0, 3, 3942, 13278, 362, 372,
                                                 4482, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14322, 0, 3, 3960, 13314, 372, 382,
                                                 4512, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14382, 0, 3, 3978, 13350, 382, 392,
                                                 4542, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14442, 0, 3, 3996, 13386, 392, 402,
                                                 4572, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14502, 0, 3, 4014, 13422, 402, 412,
                                                 4602, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14562, 0, 3, 4032, 13458, 412, 422,
                                                 4632, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14622, 0, 3, 4086, 13566, 442, 452,
                                                 4692, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14682, 0, 3, 4104, 13602, 452, 462,
                                                 4722, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14742, 0, 3, 4122, 13638, 462, 472,
                                                 4752, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14802, 0, 3, 4140, 13674, 472, 482,
                                                 4782, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14862, 0, 3, 4158, 13710, 482, 492,
                                                 4812, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14922, 0, 3, 4176, 13746, 492, 502,
                                                 4842, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14982, 0, 3, 4194, 13782, 502, 512,
                                                 4872, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15042, 0, 3, 4212, 13818, 512, 522,
                                                 4902, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15102, 0, 3, 4230, 13854, 522, 532,
                                                 4932, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15162, 0, 3, 4248, 13890, 532, 542,
                                                 4962, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15222, 0, 3, 4266, 13926, 542, 552,
                                                 4992, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15282, 0, 3, 13026, 13062, 4302, 13962,
                                                 572, 587, 5022, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15372, 0, 3, 13062, 13098, 4332, 14022,
                                                 587, 602, 5067, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15462, 0, 3, 13098, 13134, 4362, 14082,
                                                 602, 617, 5112, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15552, 0, 3, 13134, 13170, 4392, 14142,
                                                 617, 632, 5157, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15642, 0, 3, 13170, 13206, 4422, 14202,
                                                 632, 647, 5202, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15732, 0, 3, 13206, 13242, 4452, 14262,
                                                 647, 662, 5247, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15822, 0, 3, 13242, 13278, 4482, 14322,
                                                 662, 677, 5292, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 15912, 0, 3, 13278, 13314, 4512, 14382,
                                                 677, 692, 5337, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16002, 0, 3, 13314, 13350, 4542, 14442,
                                                 692, 707, 5382, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16092, 0, 3, 13350, 13386, 4572, 14502,
                                                 707, 722, 5427, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16182, 0, 3, 13386, 13422, 4602, 14562,
                                                 722, 737, 5472, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16272, 0, 3, 13494, 13530, 4662, 14622,
                                                 767, 782, 5517, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16362, 0, 3, 13530, 13566, 4692, 14682,
                                                 782, 797, 5562, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16452, 0, 3, 13566, 13602, 4722, 14742,
                                                 797, 812, 5607, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16542, 0, 3, 13602, 13638, 4752, 14802,
                                                 812, 827, 5652, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16632, 0, 3, 13638, 13674, 4782, 14862,
                                                 827, 842, 5697, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16722, 0, 3, 13674, 13710, 4812, 14922,
                                                 842, 857, 5742, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16812, 0, 3, 13710, 13746, 4842, 14982,
                                                 857, 872, 5787, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16902, 0, 3, 13746, 13782, 4872, 15042,
                                                 872, 887, 5832, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16992, 0, 3, 13782, 13818, 4902, 15102,
                                                 887, 902, 5877, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17082, 0, 3, 13818, 13854, 4932, 15162,
                                                 902, 917, 5922, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17172, 0, 3, 13854, 13890, 4962, 15222,
                                                 917, 932, 5967, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17262, 0, 3, 13962, 14022, 5067, 15462,
                                                 962, 983, 6075, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17388, 0, 3, 14022, 14082, 5112, 15552,
                                                 983, 1004, 6138, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17514, 0, 3, 14082, 14142, 5157, 15642,
                                                 1004, 1025, 6201, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17640, 0, 3, 14142, 14202, 5202, 15732,
                                                 1025, 1046, 6264, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17766, 0, 3, 14202, 14262, 5247, 15822,
                                                 1046, 1067, 6327, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17892, 0, 3, 14262, 14322, 5292, 15912,
                                                 1067, 1088, 6390, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18018, 0, 3, 14322, 14382, 5337, 16002,
                                                 1088, 1109, 6453, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18144, 0, 3, 14382, 14442, 5382, 16092,
                                                 1109, 1130, 6516, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18270, 0, 3, 14442, 14502, 5427, 16182,
                                                 1130, 1151, 6579, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18396, 0, 3, 14622, 14682, 5562, 16452,
                                                 1193, 1214, 6705, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18522, 0, 3, 14682, 14742, 5607, 16542,
                                                 1214, 1235, 6768, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18648, 0, 3, 14742, 14802, 5652, 16632,
                                                 1235, 1256, 6831, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18774, 0, 3, 14802, 14862, 5697, 16722,
                                                 1256, 1277, 6894, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18900, 0, 3, 14862, 14922, 5742, 16812,
                                                 1277, 1298, 6957, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19026, 0, 3, 14922, 14982, 5787, 16902,
                                                 1298, 1319, 7020, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19152, 0, 3, 14982, 15042, 5832, 16992,
                                                 1319, 1340, 7083, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19278, 0, 3, 15042, 15102, 5877, 17082,
                                                 1340, 1361, 7146, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19404, 0, 3, 15102, 15162, 5922, 17172,
                                                 1361, 1382, 7209, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19530, 0, 3, 15282, 15372, 6012, 17262,
                                                 1424, 1452, 7272, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19698, 0, 3, 15372, 15462, 6075, 17388,
                                                 1452, 1480, 7356, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19866, 0, 3, 15462, 15552, 6138, 17514,
                                                 1480, 1508, 7440, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20034, 0, 3, 15552, 15642, 6201, 17640,
                                                 1508, 1536, 7524, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20202, 0, 3, 15642, 15732, 6264, 17766,
                                                 1536, 1564, 7608, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20370, 0, 3, 15732, 15822, 6327, 17892,
                                                 1564, 1592, 7692, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20538, 0, 3, 15822, 15912, 6390, 18018,
                                                 1592, 1620, 7776, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20706, 0, 3, 15912, 16002, 6453, 18144,
                                                 1620, 1648, 7860, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20874, 0, 3, 16002, 16092, 6516, 18270,
                                                 1648, 1676, 7944, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21042, 0, 3, 16272, 16362, 6642, 18396,
                                                 1732, 1760, 8028, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21210, 0, 3, 16362, 16452, 6705, 18522,
                                                 1760, 1788, 8112, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21378, 0, 3, 16452, 16542, 6768, 18648,
                                                 1788, 1816, 8196, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21546, 0, 3, 16542, 16632, 6831, 18774,
                                                 1816, 1844, 8280, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21714, 0, 3, 16632, 16722, 6894, 18900,
                                                 1844, 1872, 8364, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21882, 0, 3, 16722, 16812, 6957, 19026,
                                                 1872, 1900, 8448, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 22050, 0, 3, 16812, 16902, 7020, 19152,
                                                 1900, 1928, 8532, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 22218, 0, 3, 16902, 16992, 7083, 19278,
                                                 1928, 1956, 8616, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 22386, 0, 3, 16992, 17082, 7146, 19404,
                                                 1956, 1984, 8700, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22554, 0, 3, 17262, 17388, 7356, 19866,
                                                 2040, 2076, 8892, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22770, 0, 3, 17388, 17514, 7440, 20034,
                                                 2076, 2112, 9000, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22986, 0, 3, 17514, 17640, 7524, 20202,
                                                 2112, 2148, 9108, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 23202, 0, 3, 17640, 17766, 7608, 20370,
                                                 2148, 2184, 9216, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 23418, 0, 3, 17766, 17892, 7692, 20538,
                                                 2184, 2220, 9324, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 23634, 0, 3, 17892, 18018, 7776, 20706,
                                                 2220, 2256, 9432, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 23850, 0, 3, 18018, 18144, 7860, 20874,
                                                 2256, 2292, 9540, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24066, 0, 3, 18396, 18522, 8112, 21378,
                                                 2364, 2400, 9756, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24282, 0, 3, 18522, 18648, 8196, 21546,
                                                 2400, 2436, 9864, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24498, 0, 3, 18648, 18774, 8280, 21714,
                                                 2436, 2472, 9972, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24714, 0, 3, 18774, 18900, 8364, 21882,
                                                 2472, 2508, 10080, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24930, 0, 3, 18900, 19026, 8448, 22050,
                                                 2508, 2544, 10188, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 25146, 0, 3, 19026, 19152, 8532, 22218,
                                                 2544, 2580, 10296, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 25362, 0, 3, 19152, 19278, 8616, 22386,
                                                 2580, 2616, 10404, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 25578, 0, 3, 19530, 19698, 8784, 22554,
                                                 2688, 2733, 10512, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 25848, 0, 3, 19698, 19866, 8892, 22770,
                                                 2733, 2778, 10647, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 26118, 0, 3, 19866, 20034, 9000, 22986,
                                                 2778, 2823, 10782, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 26388, 0, 3, 20034, 20202, 9108, 23202,
                                                 2823, 2868, 10917, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 26658, 0, 3, 20202, 20370, 9216, 23418,
                                                 2868, 2913, 11052, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 26928, 0, 3, 20370, 20538, 9324, 23634,
                                                 2913, 2958, 11187, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 27198, 0, 3, 20538, 20706, 9432, 23850,
                                                 2958, 3003, 11322, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 27468, 0, 3, 21042, 21210, 9648, 24066,
                                                 3093, 3138, 11457, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 27738, 0, 3, 21210, 21378, 9756, 24282,
                                                 3138, 3183, 11592, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 28008, 0, 3, 21378, 21546, 9864, 24498,
                                                 3183, 3228, 11727, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 28278, 0, 3, 21546, 21714, 9972, 24714,
                                                 3228, 3273, 11862, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 28548, 0, 3, 21714, 21882, 10080, 24930,
                                                 3273, 3318, 11997, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 28818, 0, 3, 21882, 22050, 10188, 25146,
                                                 3318, 3363, 12132, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 29088, 0, 3, 22050, 22218, 10296, 25362,
                                                 3363, 3408, 12267, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29358, 3, 3498, 3501, 12408, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29368, 3, 3501, 3504, 12414, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29378, 3, 3504, 3507, 12420, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29388, 3, 3507, 3510, 12426, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29398, 3, 3510, 3513, 12432, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29408, 3, 3513, 3516, 12438, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29418, 3, 3516, 3519, 12444, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29428, 3, 3519, 3522, 12450, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29438, 3, 3522, 3525, 12456, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29448, 3, 3525, 3528, 12462, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29458, 3, 3528, 3531, 12468, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29468, 3, 3531, 3534, 12474, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29478, 3, 3540, 3543, 12486, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29488, 3, 3543, 3546, 12492, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29498, 3, 3546, 3549, 12498, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29508, 3, 3549, 3552, 12504, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29518, 3, 3552, 3555, 12510, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29528, 3, 3555, 3558, 12516, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29538, 3, 3558, 3561, 12522, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29548, 3, 3561, 3564, 12528, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29558, 3, 3564, 3567, 12534, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29568, 3, 3567, 3570, 12540, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29578, 3, 3570, 3573, 12546, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 29588, 3, 3573, 3576, 12552, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 29598, 0, 3, 12402, 29358, 12576, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29628, 0, 3, 12408, 29368, 12594, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29658, 0, 3, 12414, 29378, 12612, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29688, 0, 3, 12420, 29388, 12630, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29718, 0, 3, 12426, 29398, 12648, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29748, 0, 3, 12432, 29408, 12666, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29778, 0, 3, 12438, 29418, 12684, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29808, 0, 3, 12444, 29428, 12702, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29838, 0, 3, 12450, 29438, 12720, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29868, 0, 3, 12456, 29448, 12738, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29898, 0, 3, 12462, 29458, 12756, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29928, 0, 3, 12468, 29468, 12774, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29958, 0, 3, 12480, 29478, 12810, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 29988, 0, 3, 12486, 29488, 12828, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30018, 0, 3, 12492, 29498, 12846, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30048, 0, 3, 12498, 29508, 12864, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30078, 0, 3, 12504, 29518, 12882, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30108, 0, 3, 12510, 29528, 12900, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30138, 0, 3, 12516, 29538, 12918, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30168, 0, 3, 12522, 29548, 12936, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30198, 0, 3, 12528, 29558, 12954, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30228, 0, 3, 12534, 29568, 12972, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30258, 0, 3, 12540, 29578, 12990, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 30288, 0, 3, 12546, 29588, 13008, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 30318, 0, 3, 12576, 29628, 3834, 3852,
                                                 13098, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30378, 0, 3, 12594, 29658, 3852, 3870,
                                                 13134, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30438, 0, 3, 12612, 29688, 3870, 3888,
                                                 13170, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30498, 0, 3, 12630, 29718, 3888, 3906,
                                                 13206, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30558, 0, 3, 12648, 29748, 3906, 3924,
                                                 13242, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30618, 0, 3, 12666, 29778, 3924, 3942,
                                                 13278, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30678, 0, 3, 12684, 29808, 3942, 3960,
                                                 13314, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30738, 0, 3, 12702, 29838, 3960, 3978,
                                                 13350, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30798, 0, 3, 12720, 29868, 3978, 3996,
                                                 13386, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30858, 0, 3, 12738, 29898, 3996, 4014,
                                                 13422, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30918, 0, 3, 12756, 29928, 4014, 4032,
                                                 13458, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 30978, 0, 3, 12810, 29988, 4068, 4086,
                                                 13566, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31038, 0, 3, 12828, 30018, 4086, 4104,
                                                 13602, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31098, 0, 3, 12846, 30048, 4104, 4122,
                                                 13638, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31158, 0, 3, 12864, 30078, 4122, 4140,
                                                 13674, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31218, 0, 3, 12882, 30108, 4140, 4158,
                                                 13710, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31278, 0, 3, 12900, 30138, 4158, 4176,
                                                 13746, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31338, 0, 3, 12918, 30168, 4176, 4194,
                                                 13782, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31398, 0, 3, 12936, 30198, 4194, 4212,
                                                 13818, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31458, 0, 3, 12954, 30228, 4212, 4230,
                                                 13854, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31518, 0, 3, 12972, 30258, 4230, 4248,
                                                 13890, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 31578, 0, 3, 12990, 30288, 4248, 4266,
                                                 13926, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 31638, 0, 3, 13098, 30378, 4302, 4332,
                                                 14022, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 31738, 0, 3, 13134, 30438, 4332, 4362,
                                                 14082, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 31838, 0, 3, 13170, 30498, 4362, 4392,
                                                 14142, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 31938, 0, 3, 13206, 30558, 4392, 4422,
                                                 14202, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32038, 0, 3, 13242, 30618, 4422, 4452,
                                                 14262, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32138, 0, 3, 13278, 30678, 4452, 4482,
                                                 14322, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32238, 0, 3, 13314, 30738, 4482, 4512,
                                                 14382, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32338, 0, 3, 13350, 30798, 4512, 4542,
                                                 14442, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32438, 0, 3, 13386, 30858, 4542, 4572,
                                                 14502, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32538, 0, 3, 13422, 30918, 4572, 4602,
                                                 14562, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32638, 0, 3, 13566, 31038, 4662, 4692,
                                                 14682, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32738, 0, 3, 13602, 31098, 4692, 4722,
                                                 14742, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32838, 0, 3, 13638, 31158, 4722, 4752,
                                                 14802, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 32938, 0, 3, 13674, 31218, 4752, 4782,
                                                 14862, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33038, 0, 3, 13710, 31278, 4782, 4812,
                                                 14922, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33138, 0, 3, 13746, 31338, 4812, 4842,
                                                 14982, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33238, 0, 3, 13782, 31398, 4842, 4872,
                                                 15042, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33338, 0, 3, 13818, 31458, 4872, 4902,
                                                 15102, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33438, 0, 3, 13854, 31518, 4902, 4932,
                                                 15162, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33538, 0, 3, 13890, 31578, 4932, 4962,
                                                 15222, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 33638, 0, 3, 30318, 30378, 14022, 31738,
                                                 5022, 5067, 15462, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 33788, 0, 3, 30378, 30438, 14082, 31838,
                                                 5067, 5112, 15552, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 33938, 0, 3, 30438, 30498, 14142, 31938,
                                                 5112, 5157, 15642, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 34088, 0, 3, 30498, 30558, 14202, 32038,
                                                 5157, 5202, 15732, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 34238, 0, 3, 30558, 30618, 14262, 32138,
                                                 5202, 5247, 15822, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 34388, 0, 3, 30618, 30678, 14322, 32238,
                                                 5247, 5292, 15912, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 34538, 0, 3, 30678, 30738, 14382, 32338,
                                                 5292, 5337, 16002, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 34688, 0, 3, 30738, 30798, 14442, 32438,
                                                 5337, 5382, 16092, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 34838, 0, 3, 30798, 30858, 14502, 32538,
                                                 5382, 5427, 16182, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 34988, 0, 3, 30978, 31038, 14682, 32738,
                                                 5517, 5562, 16452, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35138, 0, 3, 31038, 31098, 14742, 32838,
                                                 5562, 5607, 16542, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35288, 0, 3, 31098, 31158, 14802, 32938,
                                                 5607, 5652, 16632, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35438, 0, 3, 31158, 31218, 14862, 33038,
                                                 5652, 5697, 16722, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35588, 0, 3, 31218, 31278, 14922, 33138,
                                                 5697, 5742, 16812, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35738, 0, 3, 31278, 31338, 14982, 33238,
                                                 5742, 5787, 16902, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35888, 0, 3, 31338, 31398, 15042, 33338,
                                                 5787, 5832, 16992, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36038, 0, 3, 31398, 31458, 15102, 33438,
                                                 5832, 5877, 17082, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36188, 0, 3, 31458, 31518, 15162, 33538,
                                                 5877, 5922, 17172, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 36338, 0, 3, 31638, 31738, 15462, 33788,
                                                 6012, 6075, 17388, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 36548, 0, 3, 31738, 31838, 15552, 33938,
                                                 6075, 6138, 17514, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 36758, 0, 3, 31838, 31938, 15642, 34088,
                                                 6138, 6201, 17640, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 36968, 0, 3, 31938, 32038, 15732, 34238,
                                                 6201, 6264, 17766, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 37178, 0, 3, 32038, 32138, 15822, 34388,
                                                 6264, 6327, 17892, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 37388, 0, 3, 32138, 32238, 15912, 34538,
                                                 6327, 6390, 18018, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 37598, 0, 3, 32238, 32338, 16002, 34688,
                                                 6390, 6453, 18144, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 37808, 0, 3, 32338, 32438, 16092, 34838,
                                                 6453, 6516, 18270, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38018, 0, 3, 32638, 32738, 16452, 35138,
                                                 6642, 6705, 18522, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38228, 0, 3, 32738, 32838, 16542, 35288,
                                                 6705, 6768, 18648, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38438, 0, 3, 32838, 32938, 16632, 35438,
                                                 6768, 6831, 18774, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38648, 0, 3, 32938, 33038, 16722, 35588,
                                                 6831, 6894, 18900, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38858, 0, 3, 33038, 33138, 16812, 35738,
                                                 6894, 6957, 19026, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 39068, 0, 3, 33138, 33238, 16902, 35888,
                                                 6957, 7020, 19152, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 39278, 0, 3, 33238, 33338, 16992, 36038,
                                                 7020, 7083, 19278, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 39488, 0, 3, 33338, 33438, 17082, 36188,
                                                 7083, 7146, 19404, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 39698, 0, 3, 33638, 33788, 17388, 36548,
                                                 7272, 7356, 19866, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 39978, 0, 3, 33788, 33938, 17514, 36758,
                                                 7356, 7440, 20034, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 40258, 0, 3, 33938, 34088, 17640, 36968,
                                                 7440, 7524, 20202, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 40538, 0, 3, 34088, 34238, 17766, 37178,
                                                 7524, 7608, 20370, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 40818, 0, 3, 34238, 34388, 17892, 37388,
                                                 7608, 7692, 20538, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 41098, 0, 3, 34388, 34538, 18018, 37598,
                                                 7692, 7776, 20706, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 41378, 0, 3, 34538, 34688, 18144, 37808,
                                                 7776, 7860, 20874, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 41658, 0, 3, 34988, 35138, 18522, 38228,
                                                 8028, 8112, 21378, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 41938, 0, 3, 35138, 35288, 18648, 38438,
                                                 8112, 8196, 21546, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 42218, 0, 3, 35288, 35438, 18774, 38648,
                                                 8196, 8280, 21714, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 42498, 0, 3, 35438, 35588, 18900, 38858,
                                                 8280, 8364, 21882, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 42778, 0, 3, 35588, 35738, 19026, 39068,
                                                 8364, 8448, 22050, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 43058, 0, 3, 35738, 35888, 19152, 39278,
                                                 8448, 8532, 22218, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 43338, 0, 3, 35888, 36038, 19278, 39488,
                                                 8532, 8616, 22386, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 43618, 0, 3, 36338, 36548, 19866, 39978,
                                                 8784, 8892, 22770, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 43978, 0, 3, 36548, 36758, 20034, 40258,
                                                 8892, 9000, 22986, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 44338, 0, 3, 36758, 36968, 20202, 40538,
                                                 9000, 9108, 23202, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 44698, 0, 3, 36968, 37178, 20370, 40818,
                                                 9108, 9216, 23418, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 45058, 0, 3, 37178, 37388, 20538, 41098,
                                                 9216, 9324, 23634, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 45418, 0, 3, 37388, 37598, 20706, 41378,
                                                 9324, 9432, 23850, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 45778, 0, 3, 38018, 38228, 21378, 41938,
                                                 9648, 9756, 24282, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 46138, 0, 3, 38228, 38438, 21546, 42218,
                                                 9756, 9864, 24498, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 46498, 0, 3, 38438, 38648, 21714, 42498,
                                                 9864, 9972, 24714, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 46858, 0, 3, 38648, 38858, 21882, 42778,
                                                 9972, 10080, 24930, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 47218, 0, 3, 38858, 39068, 22050, 43058,
                                                 10080, 10188, 25146, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 47578, 0, 3, 39068, 39278, 22218, 43338,
                                                 10188, 10296, 25362, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 47938, 0, 3, 39698, 39978, 22770, 43978,
                                                 10512, 10647, 26118, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 48388, 0, 3, 39978, 40258, 22986, 44338,
                                                 10647, 10782, 26388, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 48838, 0, 3, 40258, 40538, 23202, 44698,
                                                 10782, 10917, 26658, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 49288, 0, 3, 40538, 40818, 23418, 45058,
                                                 10917, 11052, 26928, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 49738, 0, 3, 40818, 41098, 23634, 45418,
                                                 11052, 11187, 27198, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 50188, 0, 3, 41658, 41938, 24282, 46138,
                                                 11457, 11592, 28008, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 50638, 0, 3, 41938, 42218, 24498, 46498,
                                                 11592, 11727, 28278, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 51088, 0, 3, 42218, 42498, 24714, 46858,
                                                 11727, 11862, 28548, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 51538, 0, 3, 42498, 42778, 24930, 47218,
                                                 11862, 11997, 28818, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 51988, 0, 3, 42778, 43058, 25146, 47578,
                                                 11997, 12132, 29088, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52438, 3, 12402, 12408, 29368, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52453, 3, 12408, 12414, 29378, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52468, 3, 12414, 12420, 29388, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52483, 3, 12420, 12426, 29398, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52498, 3, 12426, 12432, 29408, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52513, 3, 12432, 12438, 29418, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52528, 3, 12438, 12444, 29428, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52543, 3, 12444, 12450, 29438, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52558, 3, 12450, 12456, 29448, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52573, 3, 12456, 12462, 29458, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52588, 3, 12462, 12468, 29468, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52603, 3, 12480, 12486, 29488, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52618, 3, 12486, 12492, 29498, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52633, 3, 12492, 12498, 29508, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52648, 3, 12498, 12504, 29518, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52663, 3, 12504, 12510, 29528, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52678, 3, 12510, 12516, 29538, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52693, 3, 12516, 12522, 29548, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52708, 3, 12522, 12528, 29558, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52723, 3, 12528, 12534, 29568, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52738, 3, 12534, 12540, 29578, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 52753, 3, 12540, 12546, 29588, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 52768, 0, 3, 29358, 52438, 29628, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 52813, 0, 3, 29368, 52453, 29658, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 52858, 0, 3, 29378, 52468, 29688, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 52903, 0, 3, 29388, 52483, 29718, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 52948, 0, 3, 29398, 52498, 29748, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 52993, 0, 3, 29408, 52513, 29778, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53038, 0, 3, 29418, 52528, 29808, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53083, 0, 3, 29428, 52543, 29838, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53128, 0, 3, 29438, 52558, 29868, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53173, 0, 3, 29448, 52573, 29898, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53218, 0, 3, 29458, 52588, 29928, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53263, 0, 3, 29478, 52603, 29988, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53308, 0, 3, 29488, 52618, 30018, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53353, 0, 3, 29498, 52633, 30048, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53398, 0, 3, 29508, 52648, 30078, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53443, 0, 3, 29518, 52663, 30108, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53488, 0, 3, 29528, 52678, 30138, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53533, 0, 3, 29538, 52693, 30168, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53578, 0, 3, 29548, 52708, 30198, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53623, 0, 3, 29558, 52723, 30228, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53668, 0, 3, 29568, 52738, 30258, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 53713, 0, 3, 29578, 52753, 30288, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 53758, 0, 3, 29598, 52768, 13026, 13062,
                                                 30318, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 53848, 0, 3, 29628, 52813, 13062, 13098,
                                                 30378, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 53938, 0, 3, 29658, 52858, 13098, 13134,
                                                 30438, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54028, 0, 3, 29688, 52903, 13134, 13170,
                                                 30498, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54118, 0, 3, 29718, 52948, 13170, 13206,
                                                 30558, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54208, 0, 3, 29748, 52993, 13206, 13242,
                                                 30618, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54298, 0, 3, 29778, 53038, 13242, 13278,
                                                 30678, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54388, 0, 3, 29808, 53083, 13278, 13314,
                                                 30738, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54478, 0, 3, 29838, 53128, 13314, 13350,
                                                 30798, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54568, 0, 3, 29868, 53173, 13350, 13386,
                                                 30858, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54658, 0, 3, 29898, 53218, 13386, 13422,
                                                 30918, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54748, 0, 3, 29958, 53263, 13494, 13530,
                                                 30978, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54838, 0, 3, 29988, 53308, 13530, 13566,
                                                 31038, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 54928, 0, 3, 30018, 53353, 13566, 13602,
                                                 31098, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55018, 0, 3, 30048, 53398, 13602, 13638,
                                                 31158, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55108, 0, 3, 30078, 53443, 13638, 13674,
                                                 31218, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55198, 0, 3, 30108, 53488, 13674, 13710,
                                                 31278, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55288, 0, 3, 30138, 53533, 13710, 13746,
                                                 31338, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55378, 0, 3, 30168, 53578, 13746, 13782,
                                                 31398, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55468, 0, 3, 30198, 53623, 13782, 13818,
                                                 31458, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55558, 0, 3, 30228, 53668, 13818, 13854,
                                                 31518, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55648, 0, 3, 30258, 53713, 13854, 13890,
                                                 31578, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 55738, 0, 3, 30378, 53938, 13962, 14022,
                                                 31738, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 55888, 0, 3, 30438, 54028, 14022, 14082,
                                                 31838, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 56038, 0, 3, 30498, 54118, 14082, 14142,
                                                 31938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 56188, 0, 3, 30558, 54208, 14142, 14202,
                                                 32038, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 56338, 0, 3, 30618, 54298, 14202, 14262,
                                                 32138, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 56488, 0, 3, 30678, 54388, 14262, 14322,
                                                 32238, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 56638, 0, 3, 30738, 54478, 14322, 14382,
                                                 32338, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 56788, 0, 3, 30798, 54568, 14382, 14442,
                                                 32438, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 56938, 0, 3, 30858, 54658, 14442, 14502,
                                                 32538, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57088, 0, 3, 31038, 54928, 14622, 14682,
                                                 32738, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57238, 0, 3, 31098, 55018, 14682, 14742,
                                                 32838, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57388, 0, 3, 31158, 55108, 14742, 14802,
                                                 32938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57538, 0, 3, 31218, 55198, 14802, 14862,
                                                 33038, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57688, 0, 3, 31278, 55288, 14862, 14922,
                                                 33138, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57838, 0, 3, 31338, 55378, 14922, 14982,
                                                 33238, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57988, 0, 3, 31398, 55468, 14982, 15042,
                                                 33338, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58138, 0, 3, 31458, 55558, 15042, 15102,
                                                 33438, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58288, 0, 3, 31518, 55648, 15102, 15162,
                                                 33538, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 58438, 0, 3, 53758, 53848, 31638, 55738,
                                                 15282, 15372, 33638, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 58663, 0, 3, 53848, 53938, 31738, 55888,
                                                 15372, 15462, 33788, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 58888, 0, 3, 53938, 54028, 31838, 56038,
                                                 15462, 15552, 33938, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 59113, 0, 3, 54028, 54118, 31938, 56188,
                                                 15552, 15642, 34088, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 59338, 0, 3, 54118, 54208, 32038, 56338,
                                                 15642, 15732, 34238, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 59563, 0, 3, 54208, 54298, 32138, 56488,
                                                 15732, 15822, 34388, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 59788, 0, 3, 54298, 54388, 32238, 56638,
                                                 15822, 15912, 34538, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60013, 0, 3, 54388, 54478, 32338, 56788,
                                                 15912, 16002, 34688, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60238, 0, 3, 54478, 54568, 32438, 56938,
                                                 16002, 16092, 34838, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60463, 0, 3, 54748, 54838, 32638, 57088,
                                                 16272, 16362, 34988, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60688, 0, 3, 54838, 54928, 32738, 57238,
                                                 16362, 16452, 35138, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60913, 0, 3, 54928, 55018, 32838, 57388,
                                                 16452, 16542, 35288, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61138, 0, 3, 55018, 55108, 32938, 57538,
                                                 16542, 16632, 35438, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61363, 0, 3, 55108, 55198, 33038, 57688,
                                                 16632, 16722, 35588, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61588, 0, 3, 55198, 55288, 33138, 57838,
                                                 16722, 16812, 35738, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61813, 0, 3, 55288, 55378, 33238, 57988,
                                                 16812, 16902, 35888, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 62038, 0, 3, 55378, 55468, 33338, 58138,
                                                 16902, 16992, 36038, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 62263, 0, 3, 55468, 55558, 33438, 58288,
                                                 16992, 17082, 36188, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 62488, 0, 3, 55738, 55888, 33788, 58888,
                                                 17262, 17388, 36548, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 62803, 0, 3, 55888, 56038, 33938, 59113,
                                                 17388, 17514, 36758, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 63118, 0, 3, 56038, 56188, 34088, 59338,
                                                 17514, 17640, 36968, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 63433, 0, 3, 56188, 56338, 34238, 59563,
                                                 17640, 17766, 37178, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 63748, 0, 3, 56338, 56488, 34388, 59788,
                                                 17766, 17892, 37388, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 64063, 0, 3, 56488, 56638, 34538, 60013,
                                                 17892, 18018, 37598, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 64378, 0, 3, 56638, 56788, 34688, 60238,
                                                 18018, 18144, 37808, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 64693, 0, 3, 57088, 57238, 35138, 60913,
                                                 18396, 18522, 38228, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 65008, 0, 3, 57238, 57388, 35288, 61138,
                                                 18522, 18648, 38438, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 65323, 0, 3, 57388, 57538, 35438, 61363,
                                                 18648, 18774, 38648, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 65638, 0, 3, 57538, 57688, 35588, 61588,
                                                 18774, 18900, 38858, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 65953, 0, 3, 57688, 57838, 35738, 61813,
                                                 18900, 19026, 39068, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 66268, 0, 3, 57838, 57988, 35888, 62038,
                                                 19026, 19152, 39278, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 66583, 0, 3, 57988, 58138, 36038, 62263,
                                                 19152, 19278, 39488, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 66898, 0, 3, 58438, 58663, 36338, 62488,
                                                 19530, 19698, 39698, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 67318, 0, 3, 58663, 58888, 36548, 62803,
                                                 19698, 19866, 39978, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 67738, 0, 3, 58888, 59113, 36758, 63118,
                                                 19866, 20034, 40258, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 68158, 0, 3, 59113, 59338, 36968, 63433,
                                                 20034, 20202, 40538, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 68578, 0, 3, 59338, 59563, 37178, 63748,
                                                 20202, 20370, 40818, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 68998, 0, 3, 59563, 59788, 37388, 64063,
                                                 20370, 20538, 41098, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 69418, 0, 3, 59788, 60013, 37598, 64378,
                                                 20538, 20706, 41378, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 69838, 0, 3, 60463, 60688, 38018, 64693,
                                                 21042, 21210, 41658, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 70258, 0, 3, 60688, 60913, 38228, 65008,
                                                 21210, 21378, 41938, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 70678, 0, 3, 60913, 61138, 38438, 65323,
                                                 21378, 21546, 42218, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 71098, 0, 3, 61138, 61363, 38648, 65638,
                                                 21546, 21714, 42498, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 71518, 0, 3, 61363, 61588, 38858, 65953,
                                                 21714, 21882, 42778, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 71938, 0, 3, 61588, 61813, 39068, 66268,
                                                 21882, 22050, 43058, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 72358, 0, 3, 61813, 62038, 39278, 66583,
                                                 22050, 22218, 43338, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 72778, 0, 3, 62488, 62803, 39978, 67738,
                                                 22554, 22770, 43978, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 73318, 0, 3, 62803, 63118, 40258, 68158,
                                                 22770, 22986, 44338, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 73858, 0, 3, 63118, 63433, 40538, 68578,
                                                 22986, 23202, 44698, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 74398, 0, 3, 63433, 63748, 40818, 68998,
                                                 23202, 23418, 45058, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 74938, 0, 3, 63748, 64063, 41098, 69418,
                                                 23418, 23634, 45418, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 75478, 0, 3, 64693, 65008, 41938, 70678,
                                                 24066, 24282, 46138, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 76018, 0, 3, 65008, 65323, 42218, 71098,
                                                 24282, 24498, 46498, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 76558, 0, 3, 65323, 65638, 42498, 71518,
                                                 24498, 24714, 46858, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 77098, 0, 3, 65638, 65953, 42778, 71938,
                                                 24714, 24930, 47218, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 77638, 0, 3, 65953, 66268, 43058, 72358,
                                                 24930, 25146, 47578, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 78178, 0, 3, 66898, 67318, 43618, 72778,
                                                 25578, 25848, 47938, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 78853, 0, 3, 67318, 67738, 43978, 73318,
                                                 25848, 26118, 48388, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 79528, 0, 3, 67738, 68158, 44338, 73858,
                                                 26118, 26388, 48838, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 80203, 0, 3, 68158, 68578, 44698, 74398,
                                                 26388, 26658, 49288, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 80878, 0, 3, 68578, 68998, 45058, 74938,
                                                 26658, 26928, 49738, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 81553, 0, 3, 69838, 70258, 45778, 75478,
                                                 27468, 27738, 50188, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 82228, 0, 3, 70258, 70678, 46138, 76018,
                                                 27738, 28008, 50638, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 82903, 0, 3, 70678, 71098, 46498, 76558,
                                                 28008, 28278, 51088, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 83578, 0, 3, 71098, 71518, 46858, 77098,
                                                 28278, 28548, 51538, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 84253, 0, 3, 71518, 71938, 47218, 77638,
                                                 28548, 28818, 51988, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 84928, 3, 29358, 29368, 52453, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 84949, 3, 29368, 29378, 52468, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 84970, 3, 29378, 29388, 52483, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 84991, 3, 29388, 29398, 52498, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85012, 3, 29398, 29408, 52513, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85033, 3, 29408, 29418, 52528, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85054, 3, 29418, 29428, 52543, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85075, 3, 29428, 29438, 52558, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85096, 3, 29438, 29448, 52573, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85117, 3, 29448, 29458, 52588, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85138, 3, 29478, 29488, 52618, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85159, 3, 29488, 29498, 52633, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85180, 3, 29498, 29508, 52648, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85201, 3, 29508, 29518, 52663, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85222, 3, 29518, 29528, 52678, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85243, 3, 29528, 29538, 52693, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85264, 3, 29538, 29548, 52708, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85285, 3, 29548, 29558, 52723, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85306, 3, 29558, 29568, 52738, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85327, 3, 29568, 29578, 52753, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 85348, 0, 3, 52438, 84928, 52813, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85411, 0, 3, 52453, 84949, 52858, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85474, 0, 3, 52468, 84970, 52903, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85537, 0, 3, 52483, 84991, 52948, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85600, 0, 3, 52498, 85012, 52993, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85663, 0, 3, 52513, 85033, 53038, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85726, 0, 3, 52528, 85054, 53083, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85789, 0, 3, 52543, 85075, 53128, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85852, 0, 3, 52558, 85096, 53173, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85915, 0, 3, 52573, 85117, 53218, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 85978, 0, 3, 52603, 85138, 53308, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86041, 0, 3, 52618, 85159, 53353, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86104, 0, 3, 52633, 85180, 53398, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86167, 0, 3, 52648, 85201, 53443, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86230, 0, 3, 52663, 85222, 53488, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86293, 0, 3, 52678, 85243, 53533, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86356, 0, 3, 52693, 85264, 53578, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86419, 0, 3, 52708, 85285, 53623, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86482, 0, 3, 52723, 85306, 53668, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86545, 0, 3, 52738, 85327, 53713, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 86608, 0, 3, 52813, 85411, 30318, 30378,
                                                 53938, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 86734, 0, 3, 52858, 85474, 30378, 30438,
                                                 54028, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 86860, 0, 3, 52903, 85537, 30438, 30498,
                                                 54118, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 86986, 0, 3, 52948, 85600, 30498, 30558,
                                                 54208, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87112, 0, 3, 52993, 85663, 30558, 30618,
                                                 54298, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87238, 0, 3, 53038, 85726, 30618, 30678,
                                                 54388, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87364, 0, 3, 53083, 85789, 30678, 30738,
                                                 54478, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87490, 0, 3, 53128, 85852, 30738, 30798,
                                                 54568, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87616, 0, 3, 53173, 85915, 30798, 30858,
                                                 54658, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87742, 0, 3, 53308, 86041, 30978, 31038,
                                                 54928, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87868, 0, 3, 53353, 86104, 31038, 31098,
                                                 55018, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87994, 0, 3, 53398, 86167, 31098, 31158,
                                                 55108, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88120, 0, 3, 53443, 86230, 31158, 31218,
                                                 55198, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88246, 0, 3, 53488, 86293, 31218, 31278,
                                                 55288, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88372, 0, 3, 53533, 86356, 31278, 31338,
                                                 55378, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88498, 0, 3, 53578, 86419, 31338, 31398,
                                                 55468, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88624, 0, 3, 53623, 86482, 31398, 31458,
                                                 55558, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88750, 0, 3, 53668, 86545, 31458, 31518,
                                                 55648, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 88876, 0, 3, 53938, 86734, 31638, 31738,
                                                 55888, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89086, 0, 3, 54028, 86860, 31738, 31838,
                                                 56038, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89296, 0, 3, 54118, 86986, 31838, 31938,
                                                 56188, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89506, 0, 3, 54208, 87112, 31938, 32038,
                                                 56338, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89716, 0, 3, 54298, 87238, 32038, 32138,
                                                 56488, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89926, 0, 3, 54388, 87364, 32138, 32238,
                                                 56638, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90136, 0, 3, 54478, 87490, 32238, 32338,
                                                 56788, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90346, 0, 3, 54568, 87616, 32338, 32438,
                                                 56938, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90556, 0, 3, 54928, 87868, 32638, 32738,
                                                 57238, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90766, 0, 3, 55018, 87994, 32738, 32838,
                                                 57388, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90976, 0, 3, 55108, 88120, 32838, 32938,
                                                 57538, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91186, 0, 3, 55198, 88246, 32938, 33038,
                                                 57688, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91396, 0, 3, 55288, 88372, 33038, 33138,
                                                 57838, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91606, 0, 3, 55378, 88498, 33138, 33238,
                                                 57988, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91816, 0, 3, 55468, 88624, 33238, 33338,
                                                 58138, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 92026, 0, 3, 55558, 88750, 33338, 33438,
                                                 58288, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 92236, 0, 3, 86608, 86734, 55888, 89086,
                                                 33638, 33788, 58888, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 92551, 0, 3, 86734, 86860, 56038, 89296,
                                                 33788, 33938, 59113, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 92866, 0, 3, 86860, 86986, 56188, 89506,
                                                 33938, 34088, 59338, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 93181, 0, 3, 86986, 87112, 56338, 89716,
                                                 34088, 34238, 59563, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 93496, 0, 3, 87112, 87238, 56488, 89926,
                                                 34238, 34388, 59788, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 93811, 0, 3, 87238, 87364, 56638, 90136,
                                                 34388, 34538, 60013, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 94126, 0, 3, 87364, 87490, 56788, 90346,
                                                 34538, 34688, 60238, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 94441, 0, 3, 87742, 87868, 57238, 90766,
                                                 34988, 35138, 60913, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 94756, 0, 3, 87868, 87994, 57388, 90976,
                                                 35138, 35288, 61138, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 95071, 0, 3, 87994, 88120, 57538, 91186,
                                                 35288, 35438, 61363, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 95386, 0, 3, 88120, 88246, 57688, 91396,
                                                 35438, 35588, 61588, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 95701, 0, 3, 88246, 88372, 57838, 91606,
                                                 35588, 35738, 61813, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 96016, 0, 3, 88372, 88498, 57988, 91816,
                                                 35738, 35888, 62038, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 96331, 0, 3, 88498, 88624, 58138, 92026,
                                                 35888, 36038, 62263, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 96646, 0, 3, 88876, 89086, 58888, 92551,
                                                 36338, 36548, 62803, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 97087, 0, 3, 89086, 89296, 59113, 92866,
                                                 36548, 36758, 63118, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 97528, 0, 3, 89296, 89506, 59338, 93181,
                                                 36758, 36968, 63433, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 97969, 0, 3, 89506, 89716, 59563, 93496,
                                                 36968, 37178, 63748, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 98410, 0, 3, 89716, 89926, 59788, 93811,
                                                 37178, 37388, 64063, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 98851, 0, 3, 89926, 90136, 60013, 94126,
                                                 37388, 37598, 64378, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 99292, 0, 3, 90556, 90766, 60913, 94756,
                                                 38018, 38228, 65008, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 99733, 0, 3, 90766, 90976, 61138, 95071,
                                                 38228, 38438, 65323, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 100174, 0, 3, 90976, 91186, 61363,
                                                 95386, 38438, 38648, 65638, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 100615, 0, 3, 91186, 91396, 61588,
                                                 95701, 38648, 38858, 65953, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 101056, 0, 3, 91396, 91606, 61813,
                                                 96016, 38858, 39068, 66268, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 101497, 0, 3, 91606, 91816, 62038,
                                                 96331, 39068, 39278, 66583, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 101938, 0, 3, 92236, 92551, 62803,
                                                 97087, 39698, 39978, 67738, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 102526, 0, 3, 92551, 92866, 63118,
                                                 97528, 39978, 40258, 68158, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 103114, 0, 3, 92866, 93181, 63433,
                                                 97969, 40258, 40538, 68578, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 103702, 0, 3, 93181, 93496, 63748,
                                                 98410, 40538, 40818, 68998, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 104290, 0, 3, 93496, 93811, 64063,
                                                 98851, 40818, 41098, 69418, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 104878, 0, 3, 94441, 94756, 65008,
                                                 99733, 41658, 41938, 70678, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 105466, 0, 3, 94756, 95071, 65323,
                                                 100174, 41938, 42218, 71098, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 106054, 0, 3, 95071, 95386, 65638,
                                                 100615, 42218, 42498, 71518, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 106642, 0, 3, 95386, 95701, 65953,
                                                 101056, 42498, 42778, 71938, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 107230, 0, 3, 95701, 96016, 66268,
                                                 101497, 42778, 43058, 72358, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 107818, 0, 3, 96646, 97087, 67738,
                                                 102526, 43618, 43978, 73318, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 108574, 0, 3, 97087, 97528, 68158,
                                                 103114, 43978, 44338, 73858, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 109330, 0, 3, 97528, 97969, 68578,
                                                 103702, 44338, 44698, 74398, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 110086, 0, 3, 97969, 98410, 68998,
                                                 104290, 44698, 45058, 74938, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 110842, 0, 3, 99292, 99733, 70678,
                                                 105466, 45778, 46138, 76018, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 111598, 0, 3, 99733, 100174, 71098,
                                                 106054, 46138, 46498, 76558, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 112354, 0, 3, 100174, 100615, 71518,
                                                 106642, 46498, 46858, 77098, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 113110, 0, 3, 100615, 101056, 71938,
                                                 107230, 46858, 47218, 77638, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 113866, 0, 3, 101938, 102526, 73318,
                                                 108574, 47938, 48388, 79528, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 114811, 0, 3, 102526, 103114, 73858,
                                                 109330, 48388, 48838, 80203, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 115756, 0, 3, 103114, 103702, 74398,
                                                 110086, 48838, 49288, 80878, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 116701, 0, 3, 104878, 105466, 76018,
                                                 111598, 50188, 50638, 82903, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 117646, 0, 3, 105466, 106054, 76558,
                                                 112354, 50638, 51088, 83578, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 118591, 0, 3, 106054, 106642, 77098,
                                                 113110, 51088, 51538, 84253, ncols, alpha, beta,
                                                 p);

            compute_prim_si_electron_repulsion_0(buffer, 119536, 3, 52438, 52453, 84949, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119564, 3, 52453, 52468, 84970, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119592, 3, 52468, 52483, 84991, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119620, 3, 52483, 52498, 85012, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119648, 3, 52498, 52513, 85033, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119676, 3, 52513, 52528, 85054, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119704, 3, 52528, 52543, 85075, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119732, 3, 52543, 52558, 85096, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119760, 3, 52558, 52573, 85117, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119788, 3, 52603, 52618, 85159, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119816, 3, 52618, 52633, 85180, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119844, 3, 52633, 52648, 85201, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119872, 3, 52648, 52663, 85222, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119900, 3, 52663, 52678, 85243, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119928, 3, 52678, 52693, 85264, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119956, 3, 52693, 52708, 85285, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 119984, 3, 52708, 52723, 85306, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 120012, 3, 52723, 52738, 85327, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120040, 0, 3, 84928, 119536, 85411,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120124, 0, 3, 84949, 119564, 85474,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120208, 0, 3, 84970, 119592, 85537,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120292, 0, 3, 84991, 119620, 85600,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120376, 0, 3, 85012, 119648, 85663,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120460, 0, 3, 85033, 119676, 85726,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120544, 0, 3, 85054, 119704, 85789,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120628, 0, 3, 85075, 119732, 85852,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120712, 0, 3, 85096, 119760, 85915,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120796, 0, 3, 85138, 119788, 86041,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120880, 0, 3, 85159, 119816, 86104,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 120964, 0, 3, 85180, 119844, 86167,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 121048, 0, 3, 85201, 119872, 86230,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 121132, 0, 3, 85222, 119900, 86293,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 121216, 0, 3, 85243, 119928, 86356,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 121300, 0, 3, 85264, 119956, 86419,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 121384, 0, 3, 85285, 119984, 86482,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 121468, 0, 3, 85306, 120012, 86545,
                                                 ncols, p);

            compute_prim_di_electron_repulsion_0(buffer, 121552, 0, 3, 85348, 120040, 53758,
                                                 53848, 86608, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 121720, 0, 3, 85411, 120124, 53848,
                                                 53938, 86734, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 121888, 0, 3, 85474, 120208, 53938,
                                                 54028, 86860, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 122056, 0, 3, 85537, 120292, 54028,
                                                 54118, 86986, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 122224, 0, 3, 85600, 120376, 54118,
                                                 54208, 87112, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 122392, 0, 3, 85663, 120460, 54208,
                                                 54298, 87238, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 122560, 0, 3, 85726, 120544, 54298,
                                                 54388, 87364, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 122728, 0, 3, 85789, 120628, 54388,
                                                 54478, 87490, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 122896, 0, 3, 85852, 120712, 54478,
                                                 54568, 87616, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 123064, 0, 3, 85978, 120796, 54748,
                                                 54838, 87742, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 123232, 0, 3, 86041, 120880, 54838,
                                                 54928, 87868, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 123400, 0, 3, 86104, 120964, 54928,
                                                 55018, 87994, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 123568, 0, 3, 86167, 121048, 55018,
                                                 55108, 88120, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 123736, 0, 3, 86230, 121132, 55108,
                                                 55198, 88246, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 123904, 0, 3, 86293, 121216, 55198,
                                                 55288, 88372, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 124072, 0, 3, 86356, 121300, 55288,
                                                 55378, 88498, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 124240, 0, 3, 86419, 121384, 55378,
                                                 55468, 88624, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 124408, 0, 3, 86482, 121468, 55468,
                                                 55558, 88750, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 124576, 0, 3, 86734, 121888, 55738,
                                                 55888, 89086, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 124856, 0, 3, 86860, 122056, 55888,
                                                 56038, 89296, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 125136, 0, 3, 86986, 122224, 56038,
                                                 56188, 89506, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 125416, 0, 3, 87112, 122392, 56188,
                                                 56338, 89716, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 125696, 0, 3, 87238, 122560, 56338,
                                                 56488, 89926, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 125976, 0, 3, 87364, 122728, 56488,
                                                 56638, 90136, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 126256, 0, 3, 87490, 122896, 56638,
                                                 56788, 90346, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 126536, 0, 3, 87868, 123400, 57088,
                                                 57238, 90766, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 126816, 0, 3, 87994, 123568, 57238,
                                                 57388, 90976, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 127096, 0, 3, 88120, 123736, 57388,
                                                 57538, 91186, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 127376, 0, 3, 88246, 123904, 57538,
                                                 57688, 91396, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 127656, 0, 3, 88372, 124072, 57688,
                                                 57838, 91606, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 127936, 0, 3, 88498, 124240, 57838,
                                                 57988, 91816, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 128216, 0, 3, 88624, 124408, 57988,
                                                 58138, 92026, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 128496, 0, 3, 121552, 121720, 88876,
                                                 124576, 58438, 58663, 92236, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 128916, 0, 3, 121720, 121888, 89086,
                                                 124856, 58663, 58888, 92551, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 129336, 0, 3, 121888, 122056, 89296,
                                                 125136, 58888, 59113, 92866, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 129756, 0, 3, 122056, 122224, 89506,
                                                 125416, 59113, 59338, 93181, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 130176, 0, 3, 122224, 122392, 89716,
                                                 125696, 59338, 59563, 93496, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 130596, 0, 3, 122392, 122560, 89926,
                                                 125976, 59563, 59788, 93811, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 131016, 0, 3, 122560, 122728, 90136,
                                                 126256, 59788, 60013, 94126, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 131436, 0, 3, 123064, 123232, 90556,
                                                 126536, 60463, 60688, 94441, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 131856, 0, 3, 123232, 123400, 90766,
                                                 126816, 60688, 60913, 94756, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 132276, 0, 3, 123400, 123568, 90976,
                                                 127096, 60913, 61138, 95071, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 132696, 0, 3, 123568, 123736, 91186,
                                                 127376, 61138, 61363, 95386, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 133116, 0, 3, 123736, 123904, 91396,
                                                 127656, 61363, 61588, 95701, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 133536, 0, 3, 123904, 124072, 91606,
                                                 127936, 61588, 61813, 96016, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 133956, 0, 3, 124072, 124240, 91816,
                                                 128216, 61813, 62038, 96331, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 134376, 0, 3, 124576, 124856, 92551,
                                                 129336, 62488, 62803, 97087, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 134964, 0, 3, 124856, 125136, 92866,
                                                 129756, 62803, 63118, 97528, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 135552, 0, 3, 125136, 125416, 93181,
                                                 130176, 63118, 63433, 97969, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 136140, 0, 3, 125416, 125696, 93496,
                                                 130596, 63433, 63748, 98410, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 136728, 0, 3, 125696, 125976, 93811,
                                                 131016, 63748, 64063, 98851, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 137316, 0, 3, 126536, 126816, 94756,
                                                 132276, 64693, 65008, 99733, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 137904, 0, 3, 126816, 127096, 95071,
                                                 132696, 65008, 65323, 100174, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 138492, 0, 3, 127096, 127376, 95386,
                                                 133116, 65323, 65638, 100615, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 139080, 0, 3, 127376, 127656, 95701,
                                                 133536, 65638, 65953, 101056, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 139668, 0, 3, 127656, 127936, 96016,
                                                 133956, 65953, 66268, 101497, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 140256, 0, 3, 128496, 128916, 96646,
                                                 134376, 66898, 67318, 101938, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 141040, 0, 3, 128916, 129336, 97087,
                                                 134964, 67318, 67738, 102526, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 141824, 0, 3, 129336, 129756, 97528,
                                                 135552, 67738, 68158, 103114, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 142608, 0, 3, 129756, 130176, 97969,
                                                 136140, 68158, 68578, 103702, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 143392, 0, 3, 130176, 130596, 98410,
                                                 136728, 68578, 68998, 104290, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 144176, 0, 3, 131436, 131856, 99292,
                                                 137316, 69838, 70258, 104878, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 144960, 0, 3, 131856, 132276, 99733,
                                                 137904, 70258, 70678, 105466, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 145744, 0, 3, 132276, 132696, 100174,
                                                 138492, 70678, 71098, 106054, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 146528, 0, 3, 132696, 133116, 100615,
                                                 139080, 71098, 71518, 106642, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 147312, 0, 3, 133116, 133536, 101056,
                                                 139668, 71518, 71938, 107230, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 148096, 0, 3, 134376, 134964, 102526,
                                                 141824, 72778, 73318, 108574, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 149104, 0, 3, 134964, 135552, 103114,
                                                 142608, 73318, 73858, 109330, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 150112, 0, 3, 135552, 136140, 103702,
                                                 143392, 73858, 74398, 110086, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 151120, 0, 3, 137316, 137904, 105466,
                                                 145744, 75478, 76018, 111598, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 152128, 0, 3, 137904, 138492, 106054,
                                                 146528, 76018, 76558, 112354, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 153136, 0, 3, 138492, 139080, 106642,
                                                 147312, 76558, 77098, 113110, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 154144, 0, 3, 140256, 141040, 107818,
                                                 148096, 78178, 78853, 113866, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 155404, 0, 3, 141040, 141824, 108574,
                                                 149104, 78853, 79528, 114811, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 156664, 0, 3, 141824, 142608, 109330,
                                                 150112, 79528, 80203, 115756, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 157924, 0, 3, 144176, 144960, 110842,
                                                 151120, 81553, 82228, 116701, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 159184, 0, 3, 144960, 145744, 111598,
                                                 152128, 82228, 82903, 117646, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 160444, 0, 3, 145744, 146528, 112354,
                                                 153136, 82903, 83578, 118591, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161704, 3, 84928, 84949, 119564, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161740, 3, 84949, 84970, 119592, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161776, 3, 84970, 84991, 119620, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161812, 3, 84991, 85012, 119648, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161848, 3, 85012, 85033, 119676, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161884, 3, 85033, 85054, 119704, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161920, 3, 85054, 85075, 119732, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161956, 3, 85075, 85096, 119760, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 161992, 3, 85138, 85159, 119816, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 162028, 3, 85159, 85180, 119844, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 162064, 3, 85180, 85201, 119872, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 162100, 3, 85201, 85222, 119900, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 162136, 3, 85222, 85243, 119928, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 162172, 3, 85243, 85264, 119956, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 162208, 3, 85264, 85285, 119984, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 162244, 3, 85285, 85306, 120012, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 162280, 0, 3, 119536, 161704, 120124,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 162388, 0, 3, 119564, 161740, 120208,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 162496, 0, 3, 119592, 161776, 120292,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 162604, 0, 3, 119620, 161812, 120376,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 162712, 0, 3, 119648, 161848, 120460,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 162820, 0, 3, 119676, 161884, 120544,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 162928, 0, 3, 119704, 161920, 120628,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163036, 0, 3, 119732, 161956, 120712,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163144, 0, 3, 119788, 161992, 120880,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163252, 0, 3, 119816, 162028, 120964,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163360, 0, 3, 119844, 162064, 121048,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163468, 0, 3, 119872, 162100, 121132,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163576, 0, 3, 119900, 162136, 121216,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163684, 0, 3, 119928, 162172, 121300,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163792, 0, 3, 119956, 162208, 121384,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 163900, 0, 3, 119984, 162244, 121468,
                                                 ncols, p);

            compute_prim_dk_electron_repulsion_0(buffer, 164008, 0, 3, 120124, 162388, 86608,
                                                 86734, 121888, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 164224, 0, 3, 120208, 162496, 86734,
                                                 86860, 122056, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 164440, 0, 3, 120292, 162604, 86860,
                                                 86986, 122224, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 164656, 0, 3, 120376, 162712, 86986,
                                                 87112, 122392, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 164872, 0, 3, 120460, 162820, 87112,
                                                 87238, 122560, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 165088, 0, 3, 120544, 162928, 87238,
                                                 87364, 122728, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 165304, 0, 3, 120628, 163036, 87364,
                                                 87490, 122896, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 165520, 0, 3, 120880, 163252, 87742,
                                                 87868, 123400, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 165736, 0, 3, 120964, 163360, 87868,
                                                 87994, 123568, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 165952, 0, 3, 121048, 163468, 87994,
                                                 88120, 123736, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 166168, 0, 3, 121132, 163576, 88120,
                                                 88246, 123904, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 166384, 0, 3, 121216, 163684, 88246,
                                                 88372, 124072, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 166600, 0, 3, 121300, 163792, 88372,
                                                 88498, 124240, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 166816, 0, 3, 121384, 163900, 88498,
                                                 88624, 124408, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 167032, 0, 3, 121888, 164224, 88876,
                                                 89086, 124856, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 167392, 0, 3, 122056, 164440, 89086,
                                                 89296, 125136, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 167752, 0, 3, 122224, 164656, 89296,
                                                 89506, 125416, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 168112, 0, 3, 122392, 164872, 89506,
                                                 89716, 125696, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 168472, 0, 3, 122560, 165088, 89716,
                                                 89926, 125976, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 168832, 0, 3, 122728, 165304, 89926,
                                                 90136, 126256, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 169192, 0, 3, 123400, 165736, 90556,
                                                 90766, 126816, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 169552, 0, 3, 123568, 165952, 90766,
                                                 90976, 127096, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 169912, 0, 3, 123736, 166168, 90976,
                                                 91186, 127376, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 170272, 0, 3, 123904, 166384, 91186,
                                                 91396, 127656, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 170632, 0, 3, 124072, 166600, 91396,
                                                 91606, 127936, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 170992, 0, 3, 124240, 166816, 91606,
                                                 91816, 128216, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 171352, 0, 3, 164008, 164224, 124856,
                                                 167392, 92236, 92551, 129336, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 171892, 0, 3, 164224, 164440, 125136,
                                                 167752, 92551, 92866, 129756, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 172432, 0, 3, 164440, 164656, 125416,
                                                 168112, 92866, 93181, 130176, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 172972, 0, 3, 164656, 164872, 125696,
                                                 168472, 93181, 93496, 130596, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 173512, 0, 3, 164872, 165088, 125976,
                                                 168832, 93496, 93811, 131016, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 174052, 0, 3, 165520, 165736, 126816,
                                                 169552, 94441, 94756, 132276, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 174592, 0, 3, 165736, 165952, 127096,
                                                 169912, 94756, 95071, 132696, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 175132, 0, 3, 165952, 166168, 127376,
                                                 170272, 95071, 95386, 133116, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 175672, 0, 3, 166168, 166384, 127656,
                                                 170632, 95386, 95701, 133536, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 176212, 0, 3, 166384, 166600, 127936,
                                                 170992, 95701, 96016, 133956, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 176752, 0, 3, 167032, 167392, 129336,
                                                 171892, 96646, 97087, 134964, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 177508, 0, 3, 167392, 167752, 129756,
                                                 172432, 97087, 97528, 135552, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 178264, 0, 3, 167752, 168112, 130176,
                                                 172972, 97528, 97969, 136140, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 179020, 0, 3, 168112, 168472, 130596,
                                                 173512, 97969, 98410, 136728, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 179776, 0, 3, 169192, 169552, 132276,
                                                 174592, 99292, 99733, 137904, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 180532, 0, 3, 169552, 169912, 132696,
                                                 175132, 99733, 100174, 138492, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 181288, 0, 3, 169912, 170272, 133116,
                                                 175672, 100174, 100615, 139080, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 182044, 0, 3, 170272, 170632, 133536,
                                                 176212, 100615, 101056, 139668, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 182800, 0, 3, 171352, 171892, 134964,
                                                 177508, 101938, 102526, 141824, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 183808, 0, 3, 171892, 172432, 135552,
                                                 178264, 102526, 103114, 142608, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 184816, 0, 3, 172432, 172972, 136140,
                                                 179020, 103114, 103702, 143392, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 185824, 0, 3, 174052, 174592, 137904,
                                                 180532, 104878, 105466, 145744, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 186832, 0, 3, 174592, 175132, 138492,
                                                 181288, 105466, 106054, 146528, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 187840, 0, 3, 175132, 175672, 139080,
                                                 182044, 106054, 106642, 147312, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 188848, 0, 3, 176752, 177508, 141824,
                                                 183808, 107818, 108574, 149104, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 190144, 0, 3, 177508, 178264, 142608,
                                                 184816, 108574, 109330, 150112, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 191440, 0, 3, 179776, 180532, 145744,
                                                 186832, 110842, 111598, 152128, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 192736, 0, 3, 180532, 181288, 146528,
                                                 187840, 111598, 112354, 153136, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 194032, 0, 3, 182800, 183808, 149104,
                                                 190144, 113866, 114811, 156664, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 195652, 0, 3, 185824, 186832, 152128,
                                                 192736, 116701, 117646, 160444, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197272, 3, 119536, 119564, 161740,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197317, 3, 119564, 119592, 161776,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197362, 3, 119592, 119620, 161812,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197407, 3, 119620, 119648, 161848,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197452, 3, 119648, 119676, 161884,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197497, 3, 119676, 119704, 161920,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197542, 3, 119704, 119732, 161956,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197587, 3, 119788, 119816, 162028,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197632, 3, 119816, 119844, 162064,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197677, 3, 119844, 119872, 162100,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197722, 3, 119872, 119900, 162136,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197767, 3, 119900, 119928, 162172,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197812, 3, 119928, 119956, 162208,
                                                 ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 197857, 3, 119956, 119984, 162244,
                                                 ncols, alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 197902, 0, 3, 161704, 197272, 162388,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 198037, 0, 3, 161740, 197317, 162496,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 198172, 0, 3, 161776, 197362, 162604,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 198307, 0, 3, 161812, 197407, 162712,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 198442, 0, 3, 161848, 197452, 162820,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 198577, 0, 3, 161884, 197497, 162928,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 198712, 0, 3, 161920, 197542, 163036,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 198847, 0, 3, 161992, 197587, 163252,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 198982, 0, 3, 162028, 197632, 163360,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 199117, 0, 3, 162064, 197677, 163468,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 199252, 0, 3, 162100, 197722, 163576,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 199387, 0, 3, 162136, 197767, 163684,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 199522, 0, 3, 162172, 197812, 163792,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 199657, 0, 3, 162208, 197857, 163900,
                                                 ncols, p);

            compute_prim_dl_electron_repulsion_0(buffer, 199792, 0, 3, 162280, 197902, 121552,
                                                 121720, 164008, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 200062, 0, 3, 162388, 198037, 121720,
                                                 121888, 164224, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 200332, 0, 3, 162496, 198172, 121888,
                                                 122056, 164440, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 200602, 0, 3, 162604, 198307, 122056,
                                                 122224, 164656, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 200872, 0, 3, 162712, 198442, 122224,
                                                 122392, 164872, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 201142, 0, 3, 162820, 198577, 122392,
                                                 122560, 165088, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 201412, 0, 3, 162928, 198712, 122560,
                                                 122728, 165304, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 201682, 0, 3, 163144, 198847, 123064,
                                                 123232, 165520, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 201952, 0, 3, 163252, 198982, 123232,
                                                 123400, 165736, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 202222, 0, 3, 163360, 199117, 123400,
                                                 123568, 165952, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 202492, 0, 3, 163468, 199252, 123568,
                                                 123736, 166168, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 202762, 0, 3, 163576, 199387, 123736,
                                                 123904, 166384, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 203032, 0, 3, 163684, 199522, 123904,
                                                 124072, 166600, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 203302, 0, 3, 163792, 199657, 124072,
                                                 124240, 166816, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 203572, 0, 3, 164224, 200332, 124576,
                                                 124856, 167392, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 204022, 0, 3, 164440, 200602, 124856,
                                                 125136, 167752, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 204472, 0, 3, 164656, 200872, 125136,
                                                 125416, 168112, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 204922, 0, 3, 164872, 201142, 125416,
                                                 125696, 168472, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 205372, 0, 3, 165088, 201412, 125696,
                                                 125976, 168832, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 205822, 0, 3, 165736, 202222, 126536,
                                                 126816, 169552, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 206272, 0, 3, 165952, 202492, 126816,
                                                 127096, 169912, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 206722, 0, 3, 166168, 202762, 127096,
                                                 127376, 170272, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 207172, 0, 3, 166384, 203032, 127376,
                                                 127656, 170632, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 207622, 0, 3, 166600, 203302, 127656,
                                                 127936, 170992, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 208072, 0, 3, 199792, 200062, 167032,
                                                 203572, 128496, 128916, 171352, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 208747, 0, 3, 200062, 200332, 167392,
                                                 204022, 128916, 129336, 171892, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 209422, 0, 3, 200332, 200602, 167752,
                                                 204472, 129336, 129756, 172432, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 210097, 0, 3, 200602, 200872, 168112,
                                                 204922, 129756, 130176, 172972, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 210772, 0, 3, 200872, 201142, 168472,
                                                 205372, 130176, 130596, 173512, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 211447, 0, 3, 201682, 201952, 169192,
                                                 205822, 131436, 131856, 174052, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 212122, 0, 3, 201952, 202222, 169552,
                                                 206272, 131856, 132276, 174592, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 212797, 0, 3, 202222, 202492, 169912,
                                                 206722, 132276, 132696, 175132, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 213472, 0, 3, 202492, 202762, 170272,
                                                 207172, 132696, 133116, 175672, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 214147, 0, 3, 202762, 203032, 170632,
                                                 207622, 133116, 133536, 176212, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 214822, 0, 3, 203572, 204022, 171892,
                                                 209422, 134376, 134964, 177508, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 215767, 0, 3, 204022, 204472, 172432,
                                                 210097, 134964, 135552, 178264, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 216712, 0, 3, 204472, 204922, 172972,
                                                 210772, 135552, 136140, 179020, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 217657, 0, 3, 205822, 206272, 174592,
                                                 212797, 137316, 137904, 180532, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 218602, 0, 3, 206272, 206722, 175132,
                                                 213472, 137904, 138492, 181288, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 219547, 0, 3, 206722, 207172, 175672,
                                                 214147, 138492, 139080, 182044, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 220492, 0, 3, 208072, 208747, 176752,
                                                 214822, 140256, 141040, 182800, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 221752, 0, 3, 208747, 209422, 177508,
                                                 215767, 141040, 141824, 183808, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 223012, 0, 3, 209422, 210097, 178264,
                                                 216712, 141824, 142608, 184816, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 224272, 0, 3, 211447, 212122, 179776,
                                                 217657, 144176, 144960, 185824, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 225532, 0, 3, 212122, 212797, 180532,
                                                 218602, 144960, 145744, 186832, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 226792, 0, 3, 212797, 213472, 181288,
                                                 219547, 145744, 146528, 187840, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 228052, 0, 3, 214822, 215767, 183808,
                                                 223012, 148096, 149104, 190144, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 229672, 0, 3, 217657, 218602, 186832,
                                                 226792, 151120, 152128, 192736, ncols, alpha,
                                                 beta, p);

            compute_prim_ll_electron_repulsion_0(buffer, 231292, 0, 3, 220492, 221752, 188848,
                                                 228052, 154144, 155404, 194032, ncols, alpha,
                                                 beta, p);

            compute_prim_ll_electron_repulsion_0(buffer, 233317, 0, 3, 224272, 225532, 191440,
                                                 229672, 157924, 159184, 195652, ncols, alpha,
                                                 beta, p);

            simdfunc::contract_primitives(buffer, 235342, 231292, 4050, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 239392, 237367, 45, 1, nmax);

    simdtrf::transform_l_outer_tri(values, nvalues, buffer, 239392, nmax);

    simdtrf::transform_l_inner(buffer, 239392, 235342, 45, 1, nmax);

    simdtrf::transform_l_outer_tri(values + 289 * nvalues, nvalues, buffer, 239392, nmax);
}

}  // namespace simdt2ceri
