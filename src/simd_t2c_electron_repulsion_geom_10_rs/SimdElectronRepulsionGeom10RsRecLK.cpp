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


#include "SimdElectronRepulsionGeom10RsRecLK.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFK.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGK.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHK.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIK.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKK.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLK.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMH.hpp"
#include "SimdElectronRepulsionVrrRecMI.hpp"
#include "SimdElectronRepulsionVrrRecMK.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPK.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSK.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryL1.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_lk_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_lk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 247213, 236818, 9720, nvalues);

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
                                                9, 10, 11, 12, 13, 14, 15, 16}, ncols, fj, mu,
                                                omega);

            simdfunc::compute_boys_function(buffer, coordinates, 23, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13, 14, 15, 16}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 82, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 85, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 88, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 91, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 94, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 97, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 100, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 103, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 106, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 109, 0, 33, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 112, 0, 34, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 115, 0, 35, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 118, 0, 36, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 121, 0, 37, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 124, 0, 38, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 127, 0, 39, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 7, 8, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 8, 9, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 9, 10, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 10, 11, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 11, 12, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 12, 13, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 13, 14, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 14, 15, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 178, 0, 15, 16, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 184, 0, 16, 17, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 190, 0, 17, 18, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 196, 0, 18, 19, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 202, 0, 19, 20, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 208, 0, 20, 21, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 214, 0, 24, 25, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 220, 0, 25, 26, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 226, 0, 26, 27, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 232, 0, 27, 28, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 238, 0, 28, 29, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 244, 0, 29, 30, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 250, 0, 30, 31, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 256, 0, 31, 32, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 262, 0, 32, 33, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 268, 0, 33, 34, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 274, 0, 34, 35, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 280, 0, 35, 36, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 286, 0, 36, 37, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 292, 0, 37, 38, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 40, 43, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 43, 46, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 318, 0, 46, 49, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 328, 0, 49, 52, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 338, 0, 52, 55, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 348, 0, 55, 58, 166, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 358, 0, 58, 61, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 368, 0, 61, 64, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 378, 0, 64, 67, 184, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 388, 0, 67, 70, 190, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 398, 0, 70, 73, 196, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 408, 0, 73, 76, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 418, 0, 76, 79, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 428, 0, 85, 88, 220, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 438, 0, 88, 91, 226, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 448, 0, 91, 94, 232, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 458, 0, 94, 97, 238, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 468, 0, 97, 100, 244, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 478, 0, 100, 103, 250, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 488, 0, 103, 106, 256, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 498, 0, 106, 109, 262, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 508, 0, 109, 112, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 518, 0, 112, 115, 274, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 528, 0, 115, 118, 280, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 538, 0, 118, 121, 286, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 548, 0, 121, 124, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 558, 0, 130, 136, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 573, 0, 136, 142, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 588, 0, 142, 148, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 603, 0, 148, 154, 338, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 618, 0, 154, 160, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 633, 0, 160, 166, 358, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 648, 0, 166, 172, 368, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 663, 0, 172, 178, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 678, 0, 178, 184, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 693, 0, 184, 190, 398, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 708, 0, 190, 196, 408, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 723, 0, 196, 202, 418, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 738, 0, 214, 220, 438, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 753, 0, 220, 226, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 768, 0, 226, 232, 458, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 783, 0, 232, 238, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 798, 0, 238, 244, 478, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 813, 0, 244, 250, 488, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 828, 0, 250, 256, 498, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 843, 0, 256, 262, 508, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 858, 0, 262, 268, 518, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 873, 0, 268, 274, 528, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 888, 0, 274, 280, 538, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 903, 0, 280, 286, 548, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 918, 0, 298, 308, 573, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 939, 0, 308, 318, 588, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 960, 0, 318, 328, 603, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 981, 0, 328, 338, 618, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1002, 0, 338, 348, 633, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1023, 0, 348, 358, 648, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1044, 0, 358, 368, 663, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1065, 0, 368, 378, 678, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1086, 0, 378, 388, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1107, 0, 388, 398, 708, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1128, 0, 398, 408, 723, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1149, 0, 428, 438, 753, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1170, 0, 438, 448, 768, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1191, 0, 448, 458, 783, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1212, 0, 458, 468, 798, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1233, 0, 468, 478, 813, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1254, 0, 478, 488, 828, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1275, 0, 488, 498, 843, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1296, 0, 498, 508, 858, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1317, 0, 508, 518, 873, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1338, 0, 518, 528, 888, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1359, 0, 528, 538, 903, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1380, 0, 558, 573, 939, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1408, 0, 573, 588, 960, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1436, 0, 588, 603, 981, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1464, 0, 603, 618, 1002, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1492, 0, 618, 633, 1023, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1520, 0, 633, 648, 1044, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1548, 0, 648, 663, 1065, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1576, 0, 663, 678, 1086, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1604, 0, 678, 693, 1107, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1632, 0, 693, 708, 1128, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1660, 0, 738, 753, 1170, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1688, 0, 753, 768, 1191, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1716, 0, 768, 783, 1212, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1744, 0, 783, 798, 1233, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1772, 0, 798, 813, 1254, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1800, 0, 813, 828, 1275, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1828, 0, 828, 843, 1296, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1856, 0, 843, 858, 1317, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1884, 0, 858, 873, 1338, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1912, 0, 873, 888, 1359, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1940, 0, 918, 939, 1408, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1976, 0, 939, 960, 1436, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2012, 0, 960, 981, 1464, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2048, 0, 981, 1002, 1492, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2084, 0, 1002, 1023, 1520, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2120, 0, 1023, 1044, 1548, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2156, 0, 1044, 1065, 1576, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2192, 0, 1065, 1086, 1604, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2228, 0, 1086, 1107, 1632, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2264, 0, 1149, 1170, 1688, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2300, 0, 1170, 1191, 1716, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2336, 0, 1191, 1212, 1744, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2372, 0, 1212, 1233, 1772, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2408, 0, 1233, 1254, 1800, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2444, 0, 1254, 1275, 1828, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2480, 0, 1275, 1296, 1856, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2516, 0, 1296, 1317, 1884, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2552, 0, 1317, 1338, 1912, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2588, 0, 1380, 1408, 1976, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2633, 0, 1408, 1436, 2012, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2678, 0, 1436, 1464, 2048, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2723, 0, 1464, 1492, 2084, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2768, 0, 1492, 1520, 2120, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2813, 0, 1520, 1548, 2156, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2858, 0, 1548, 1576, 2192, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2903, 0, 1576, 1604, 2228, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2948, 0, 1660, 1688, 2300, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2993, 0, 1688, 1716, 2336, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3038, 0, 1716, 1744, 2372, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3083, 0, 1744, 1772, 2408, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3128, 0, 1772, 1800, 2444, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3173, 0, 1800, 1828, 2480, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3218, 0, 1828, 1856, 2516, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3263, 0, 1856, 1884, 2552, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3308, 0, 1940, 1976, 2633, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3363, 0, 1976, 2012, 2678, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3418, 0, 2012, 2048, 2723, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3473, 0, 2048, 2084, 2768, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3528, 0, 2084, 2120, 2813, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3583, 0, 2120, 2156, 2858, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3638, 0, 2156, 2192, 2903, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3693, 0, 2264, 2300, 2993, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3748, 0, 2300, 2336, 3038, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3803, 0, 2336, 2372, 3083, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3858, 0, 2372, 2408, 3128, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3913, 0, 2408, 2444, 3173, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3968, 0, 2444, 2480, 3218, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 4023, 0, 2480, 2516, 3263, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 4078, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4081, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4084, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4087, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4090, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4093, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4096, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4099, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4102, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4105, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4108, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4111, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4114, 3, 22, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4117, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4120, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4123, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4126, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4129, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4132, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4135, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4138, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4141, 3, 35, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4144, 3, 36, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4147, 3, 37, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4150, 3, 38, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 4153, 3, 39, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 4156, 3, 9, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4165, 3, 10, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4174, 3, 11, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4183, 3, 12, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4192, 3, 13, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4201, 3, 14, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4210, 3, 15, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4219, 3, 16, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4228, 3, 17, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4237, 3, 18, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4246, 3, 19, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4255, 3, 20, 79, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4264, 3, 21, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4273, 3, 26, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4282, 3, 27, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4291, 3, 28, 97, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4300, 3, 29, 100, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4309, 3, 30, 103, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4318, 3, 31, 106, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4327, 3, 32, 109, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4336, 3, 33, 112, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4345, 3, 34, 115, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4354, 3, 35, 118, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4363, 3, 36, 121, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4372, 3, 37, 124, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4381, 3, 38, 127, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4390, 0, 3, 43, 4156, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4408, 0, 3, 46, 4165, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4426, 0, 3, 49, 4174, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4444, 0, 3, 52, 4183, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4462, 0, 3, 55, 4192, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4480, 0, 3, 58, 4201, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4498, 0, 3, 61, 4210, 172, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4516, 0, 3, 64, 4219, 178, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4534, 0, 3, 67, 4228, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4552, 0, 3, 70, 4237, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4570, 0, 3, 73, 4246, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4588, 0, 3, 76, 4255, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4606, 0, 3, 79, 4264, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4624, 0, 3, 88, 4273, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4642, 0, 3, 91, 4282, 226, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4660, 0, 3, 94, 4291, 232, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4678, 0, 3, 97, 4300, 238, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4696, 0, 3, 100, 4309, 244, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4714, 0, 3, 103, 4318, 250, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4732, 0, 3, 106, 4327, 256, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4750, 0, 3, 109, 4336, 262, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4768, 0, 3, 112, 4345, 268, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4786, 0, 3, 115, 4354, 274, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4804, 0, 3, 118, 4363, 280, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4822, 0, 3, 121, 4372, 286, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4840, 0, 3, 124, 4381, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4858, 0, 3, 130, 4390, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4888, 0, 3, 136, 4408, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4918, 0, 3, 142, 4426, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4948, 0, 3, 148, 4444, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4978, 0, 3, 154, 4462, 338, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5008, 0, 3, 160, 4480, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5038, 0, 3, 166, 4498, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5068, 0, 3, 172, 4516, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5098, 0, 3, 178, 4534, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5128, 0, 3, 184, 4552, 388, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5158, 0, 3, 190, 4570, 398, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5188, 0, 3, 196, 4588, 408, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5218, 0, 3, 202, 4606, 418, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5248, 0, 3, 214, 4624, 428, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5278, 0, 3, 220, 4642, 438, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5308, 0, 3, 226, 4660, 448, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5338, 0, 3, 232, 4678, 458, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5368, 0, 3, 238, 4696, 468, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5398, 0, 3, 244, 4714, 478, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5428, 0, 3, 250, 4732, 488, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5458, 0, 3, 256, 4750, 498, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5488, 0, 3, 262, 4768, 508, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5518, 0, 3, 268, 4786, 518, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5548, 0, 3, 274, 4804, 528, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5578, 0, 3, 280, 4822, 538, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5608, 0, 3, 286, 4840, 548, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5638, 0, 3, 308, 4918, 573, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5683, 0, 3, 318, 4948, 588, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5728, 0, 3, 328, 4978, 603, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5773, 0, 3, 338, 5008, 618, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5818, 0, 3, 348, 5038, 633, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5863, 0, 3, 358, 5068, 648, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5908, 0, 3, 368, 5098, 663, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5953, 0, 3, 378, 5128, 678, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5998, 0, 3, 388, 5158, 693, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6043, 0, 3, 398, 5188, 708, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6088, 0, 3, 408, 5218, 723, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6133, 0, 3, 438, 5308, 753, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6178, 0, 3, 448, 5338, 768, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6223, 0, 3, 458, 5368, 783, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6268, 0, 3, 468, 5398, 798, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6313, 0, 3, 478, 5428, 813, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6358, 0, 3, 488, 5458, 828, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6403, 0, 3, 498, 5488, 843, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6448, 0, 3, 508, 5518, 858, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6493, 0, 3, 518, 5548, 873, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6538, 0, 3, 528, 5578, 888, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6583, 0, 3, 538, 5608, 903, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6628, 0, 3, 558, 5638, 918, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6691, 0, 3, 573, 5683, 939, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6754, 0, 3, 588, 5728, 960, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6817, 0, 3, 603, 5773, 981, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6880, 0, 3, 618, 5818, 1002, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6943, 0, 3, 633, 5863, 1023, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7006, 0, 3, 648, 5908, 1044, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7069, 0, 3, 663, 5953, 1065, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7132, 0, 3, 678, 5998, 1086, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7195, 0, 3, 693, 6043, 1107, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7258, 0, 3, 708, 6088, 1128, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7321, 0, 3, 738, 6133, 1149, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7384, 0, 3, 753, 6178, 1170, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7447, 0, 3, 768, 6223, 1191, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7510, 0, 3, 783, 6268, 1212, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7573, 0, 3, 798, 6313, 1233, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7636, 0, 3, 813, 6358, 1254, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7699, 0, 3, 828, 6403, 1275, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7762, 0, 3, 843, 6448, 1296, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7825, 0, 3, 858, 6493, 1317, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7888, 0, 3, 873, 6538, 1338, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7951, 0, 3, 888, 6583, 1359, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 8014, 0, 3, 939, 6754, 1408, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 8098, 0, 3, 960, 6817, 1436, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 8182, 0, 3, 981, 6880, 1464, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 8266, 0, 3, 1002, 6943, 1492, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8350, 0, 3, 1023, 7006, 1520, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8434, 0, 3, 1044, 7069, 1548, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8518, 0, 3, 1065, 7132, 1576, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8602, 0, 3, 1086, 7195, 1604, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8686, 0, 3, 1107, 7258, 1632, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8770, 0, 3, 1170, 7447, 1688, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8854, 0, 3, 1191, 7510, 1716, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8938, 0, 3, 1212, 7573, 1744, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9022, 0, 3, 1233, 7636, 1772, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9106, 0, 3, 1254, 7699, 1800, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9190, 0, 3, 1275, 7762, 1828, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9274, 0, 3, 1296, 7825, 1856, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9358, 0, 3, 1317, 7888, 1884, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 9442, 0, 3, 1338, 7951, 1912, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9526, 0, 3, 1380, 8014, 1940, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9634, 0, 3, 1408, 8098, 1976, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9742, 0, 3, 1436, 8182, 2012, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9850, 0, 3, 1464, 8266, 2048, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9958, 0, 3, 1492, 8350, 2084, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10066, 0, 3, 1520, 8434, 2120, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10174, 0, 3, 1548, 8518, 2156, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10282, 0, 3, 1576, 8602, 2192, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10390, 0, 3, 1604, 8686, 2228, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10498, 0, 3, 1660, 8770, 2264, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10606, 0, 3, 1688, 8854, 2300, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10714, 0, 3, 1716, 8938, 2336, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10822, 0, 3, 1744, 9022, 2372, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10930, 0, 3, 1772, 9106, 2408, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11038, 0, 3, 1800, 9190, 2444, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11146, 0, 3, 1828, 9274, 2480, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11254, 0, 3, 1856, 9358, 2516, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 11362, 0, 3, 1884, 9442, 2552, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11470, 0, 3, 1976, 9742, 2633, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11605, 0, 3, 2012, 9850, 2678, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11740, 0, 3, 2048, 9958, 2723, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11875, 0, 3, 2084, 10066, 2768, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12010, 0, 3, 2120, 10174, 2813, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12145, 0, 3, 2156, 10282, 2858, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12280, 0, 3, 2192, 10390, 2903, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12415, 0, 3, 2300, 10714, 2993, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12550, 0, 3, 2336, 10822, 3038, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12685, 0, 3, 2372, 10930, 3083, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12820, 0, 3, 2408, 11038, 3128, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 12955, 0, 3, 2444, 11146, 3173, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13090, 0, 3, 2480, 11254, 3218, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 13225, 0, 3, 2516, 11362, 3263, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 13360, 0, 3, 2588, 11470, 3308, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 13525, 0, 3, 2633, 11605, 3363, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 13690, 0, 3, 2678, 11740, 3418, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 13855, 0, 3, 2723, 11875, 3473, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14020, 0, 3, 2768, 12010, 3528, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14185, 0, 3, 2813, 12145, 3583, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14350, 0, 3, 2858, 12280, 3638, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14515, 0, 3, 2948, 12415, 3693, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14680, 0, 3, 2993, 12550, 3748, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 14845, 0, 3, 3038, 12685, 3803, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15010, 0, 3, 3083, 12820, 3858, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15175, 0, 3, 3128, 12955, 3913, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15340, 0, 3, 3173, 13090, 3968, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 15505, 0, 3, 3218, 13225, 4023, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 15670, 3, 9, 10, 4081, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15676, 3, 10, 11, 4084, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15682, 3, 11, 12, 4087, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15688, 3, 12, 13, 4090, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15694, 3, 13, 14, 4093, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15700, 3, 14, 15, 4096, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15706, 3, 15, 16, 4099, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15712, 3, 16, 17, 4102, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15718, 3, 17, 18, 4105, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15724, 3, 18, 19, 4108, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15730, 3, 19, 20, 4111, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15736, 3, 20, 21, 4114, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15742, 3, 26, 27, 4120, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15748, 3, 27, 28, 4123, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15754, 3, 28, 29, 4126, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15760, 3, 29, 30, 4129, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15766, 3, 30, 31, 4132, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15772, 3, 31, 32, 4135, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15778, 3, 32, 33, 4138, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15784, 3, 33, 34, 4141, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15790, 3, 34, 35, 4144, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15796, 3, 35, 36, 4147, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15802, 3, 36, 37, 4150, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 15808, 3, 37, 38, 4153, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 15814, 0, 3, 4078, 15670, 4165, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15832, 0, 3, 4081, 15676, 4174, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15850, 0, 3, 4084, 15682, 4183, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15868, 0, 3, 4087, 15688, 4192, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15886, 0, 3, 4090, 15694, 4201, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15904, 0, 3, 4093, 15700, 4210, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15922, 0, 3, 4096, 15706, 4219, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15940, 0, 3, 4099, 15712, 4228, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15958, 0, 3, 4102, 15718, 4237, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15976, 0, 3, 4105, 15724, 4246, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 15994, 0, 3, 4108, 15730, 4255, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16012, 0, 3, 4111, 15736, 4264, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16030, 0, 3, 4117, 15742, 4282, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16048, 0, 3, 4120, 15748, 4291, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16066, 0, 3, 4123, 15754, 4300, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16084, 0, 3, 4126, 15760, 4309, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16102, 0, 3, 4129, 15766, 4318, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16120, 0, 3, 4132, 15772, 4327, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16138, 0, 3, 4135, 15778, 4336, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16156, 0, 3, 4138, 15784, 4345, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16174, 0, 3, 4141, 15790, 4354, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16192, 0, 3, 4144, 15796, 4363, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16210, 0, 3, 4147, 15802, 4372, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 16228, 0, 3, 4150, 15808, 4381, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 16246, 0, 3, 4156, 15814, 130, 136,
                                                 4408, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16282, 0, 3, 4165, 15832, 136, 142,
                                                 4426, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16318, 0, 3, 4174, 15850, 142, 148,
                                                 4444, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16354, 0, 3, 4183, 15868, 148, 154,
                                                 4462, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16390, 0, 3, 4192, 15886, 154, 160,
                                                 4480, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16426, 0, 3, 4201, 15904, 160, 166,
                                                 4498, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16462, 0, 3, 4210, 15922, 166, 172,
                                                 4516, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16498, 0, 3, 4219, 15940, 172, 178,
                                                 4534, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16534, 0, 3, 4228, 15958, 178, 184,
                                                 4552, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16570, 0, 3, 4237, 15976, 184, 190,
                                                 4570, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16606, 0, 3, 4246, 15994, 190, 196,
                                                 4588, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16642, 0, 3, 4255, 16012, 196, 202,
                                                 4606, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16678, 0, 3, 4273, 16030, 214, 220,
                                                 4642, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16714, 0, 3, 4282, 16048, 220, 226,
                                                 4660, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16750, 0, 3, 4291, 16066, 226, 232,
                                                 4678, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16786, 0, 3, 4300, 16084, 232, 238,
                                                 4696, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16822, 0, 3, 4309, 16102, 238, 244,
                                                 4714, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16858, 0, 3, 4318, 16120, 244, 250,
                                                 4732, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16894, 0, 3, 4327, 16138, 250, 256,
                                                 4750, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16930, 0, 3, 4336, 16156, 256, 262,
                                                 4768, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 16966, 0, 3, 4345, 16174, 262, 268,
                                                 4786, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17002, 0, 3, 4354, 16192, 268, 274,
                                                 4804, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17038, 0, 3, 4363, 16210, 274, 280,
                                                 4822, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 17074, 0, 3, 4372, 16228, 280, 286,
                                                 4840, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17110, 0, 3, 4408, 16282, 298, 308,
                                                 4918, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17170, 0, 3, 4426, 16318, 308, 318,
                                                 4948, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17230, 0, 3, 4444, 16354, 318, 328,
                                                 4978, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17290, 0, 3, 4462, 16390, 328, 338,
                                                 5008, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17350, 0, 3, 4480, 16426, 338, 348,
                                                 5038, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17410, 0, 3, 4498, 16462, 348, 358,
                                                 5068, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17470, 0, 3, 4516, 16498, 358, 368,
                                                 5098, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17530, 0, 3, 4534, 16534, 368, 378,
                                                 5128, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17590, 0, 3, 4552, 16570, 378, 388,
                                                 5158, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17650, 0, 3, 4570, 16606, 388, 398,
                                                 5188, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17710, 0, 3, 4588, 16642, 398, 408,
                                                 5218, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17770, 0, 3, 4642, 16714, 428, 438,
                                                 5308, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17830, 0, 3, 4660, 16750, 438, 448,
                                                 5338, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17890, 0, 3, 4678, 16786, 448, 458,
                                                 5368, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 17950, 0, 3, 4696, 16822, 458, 468,
                                                 5398, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18010, 0, 3, 4714, 16858, 468, 478,
                                                 5428, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18070, 0, 3, 4732, 16894, 478, 488,
                                                 5458, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18130, 0, 3, 4750, 16930, 488, 498,
                                                 5488, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18190, 0, 3, 4768, 16966, 498, 508,
                                                 5518, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18250, 0, 3, 4786, 17002, 508, 518,
                                                 5548, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18310, 0, 3, 4804, 17038, 518, 528,
                                                 5578, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 18370, 0, 3, 4822, 17074, 528, 538,
                                                 5608, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 18430, 0, 3, 16246, 16282, 4918, 17170,
                                                 558, 573, 5683, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 18520, 0, 3, 16282, 16318, 4948, 17230,
                                                 573, 588, 5728, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 18610, 0, 3, 16318, 16354, 4978, 17290,
                                                 588, 603, 5773, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 18700, 0, 3, 16354, 16390, 5008, 17350,
                                                 603, 618, 5818, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 18790, 0, 3, 16390, 16426, 5038, 17410,
                                                 618, 633, 5863, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 18880, 0, 3, 16426, 16462, 5068, 17470,
                                                 633, 648, 5908, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 18970, 0, 3, 16462, 16498, 5098, 17530,
                                                 648, 663, 5953, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19060, 0, 3, 16498, 16534, 5128, 17590,
                                                 663, 678, 5998, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19150, 0, 3, 16534, 16570, 5158, 17650,
                                                 678, 693, 6043, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19240, 0, 3, 16570, 16606, 5188, 17710,
                                                 693, 708, 6088, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19330, 0, 3, 16678, 16714, 5308, 17830,
                                                 738, 753, 6178, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19420, 0, 3, 16714, 16750, 5338, 17890,
                                                 753, 768, 6223, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19510, 0, 3, 16750, 16786, 5368, 17950,
                                                 768, 783, 6268, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19600, 0, 3, 16786, 16822, 5398, 18010,
                                                 783, 798, 6313, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19690, 0, 3, 16822, 16858, 5428, 18070,
                                                 798, 813, 6358, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19780, 0, 3, 16858, 16894, 5458, 18130,
                                                 813, 828, 6403, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19870, 0, 3, 16894, 16930, 5488, 18190,
                                                 828, 843, 6448, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 19960, 0, 3, 16930, 16966, 5518, 18250,
                                                 843, 858, 6493, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20050, 0, 3, 16966, 17002, 5548, 18310,
                                                 858, 873, 6538, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 20140, 0, 3, 17002, 17038, 5578, 18370,
                                                 873, 888, 6583, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 20230, 0, 3, 17110, 17170, 5683, 18520,
                                                 918, 939, 6754, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 20356, 0, 3, 17170, 17230, 5728, 18610,
                                                 939, 960, 6817, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 20482, 0, 3, 17230, 17290, 5773, 18700,
                                                 960, 981, 6880, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 20608, 0, 3, 17290, 17350, 5818, 18790,
                                                 981, 1002, 6943, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 20734, 0, 3, 17350, 17410, 5863, 18880,
                                                 1002, 1023, 7006, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 20860, 0, 3, 17410, 17470, 5908, 18970,
                                                 1023, 1044, 7069, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 20986, 0, 3, 17470, 17530, 5953, 19060,
                                                 1044, 1065, 7132, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 21112, 0, 3, 17530, 17590, 5998, 19150,
                                                 1065, 1086, 7195, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 21238, 0, 3, 17590, 17650, 6043, 19240,
                                                 1086, 1107, 7258, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 21364, 0, 3, 17770, 17830, 6178, 19420,
                                                 1149, 1170, 7447, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 21490, 0, 3, 17830, 17890, 6223, 19510,
                                                 1170, 1191, 7510, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 21616, 0, 3, 17890, 17950, 6268, 19600,
                                                 1191, 1212, 7573, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 21742, 0, 3, 17950, 18010, 6313, 19690,
                                                 1212, 1233, 7636, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 21868, 0, 3, 18010, 18070, 6358, 19780,
                                                 1233, 1254, 7699, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 21994, 0, 3, 18070, 18130, 6403, 19870,
                                                 1254, 1275, 7762, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22120, 0, 3, 18130, 18190, 6448, 19960,
                                                 1275, 1296, 7825, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22246, 0, 3, 18190, 18250, 6493, 20050,
                                                 1296, 1317, 7888, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 22372, 0, 3, 18250, 18310, 6538, 20140,
                                                 1317, 1338, 7951, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 22498, 0, 3, 18430, 18520, 6754, 20356,
                                                 1380, 1408, 8098, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 22666, 0, 3, 18520, 18610, 6817, 20482,
                                                 1408, 1436, 8182, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 22834, 0, 3, 18610, 18700, 6880, 20608,
                                                 1436, 1464, 8266, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 23002, 0, 3, 18700, 18790, 6943, 20734,
                                                 1464, 1492, 8350, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 23170, 0, 3, 18790, 18880, 7006, 20860,
                                                 1492, 1520, 8434, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 23338, 0, 3, 18880, 18970, 7069, 20986,
                                                 1520, 1548, 8518, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 23506, 0, 3, 18970, 19060, 7132, 21112,
                                                 1548, 1576, 8602, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 23674, 0, 3, 19060, 19150, 7195, 21238,
                                                 1576, 1604, 8686, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 23842, 0, 3, 19330, 19420, 7447, 21490,
                                                 1660, 1688, 8854, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 24010, 0, 3, 19420, 19510, 7510, 21616,
                                                 1688, 1716, 8938, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 24178, 0, 3, 19510, 19600, 7573, 21742,
                                                 1716, 1744, 9022, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 24346, 0, 3, 19600, 19690, 7636, 21868,
                                                 1744, 1772, 9106, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 24514, 0, 3, 19690, 19780, 7699, 21994,
                                                 1772, 1800, 9190, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 24682, 0, 3, 19780, 19870, 7762, 22120,
                                                 1800, 1828, 9274, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 24850, 0, 3, 19870, 19960, 7825, 22246,
                                                 1828, 1856, 9358, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 25018, 0, 3, 19960, 20050, 7888, 22372,
                                                 1856, 1884, 9442, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 25186, 0, 3, 20230, 20356, 8098, 22666,
                                                 1940, 1976, 9742, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 25402, 0, 3, 20356, 20482, 8182, 22834,
                                                 1976, 2012, 9850, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 25618, 0, 3, 20482, 20608, 8266, 23002,
                                                 2012, 2048, 9958, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 25834, 0, 3, 20608, 20734, 8350, 23170,
                                                 2048, 2084, 10066, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 26050, 0, 3, 20734, 20860, 8434, 23338,
                                                 2084, 2120, 10174, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 26266, 0, 3, 20860, 20986, 8518, 23506,
                                                 2120, 2156, 10282, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 26482, 0, 3, 20986, 21112, 8602, 23674,
                                                 2156, 2192, 10390, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 26698, 0, 3, 21364, 21490, 8854, 24010,
                                                 2264, 2300, 10714, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 26914, 0, 3, 21490, 21616, 8938, 24178,
                                                 2300, 2336, 10822, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 27130, 0, 3, 21616, 21742, 9022, 24346,
                                                 2336, 2372, 10930, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 27346, 0, 3, 21742, 21868, 9106, 24514,
                                                 2372, 2408, 11038, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 27562, 0, 3, 21868, 21994, 9190, 24682,
                                                 2408, 2444, 11146, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 27778, 0, 3, 21994, 22120, 9274, 24850,
                                                 2444, 2480, 11254, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 27994, 0, 3, 22120, 22246, 9358, 25018,
                                                 2480, 2516, 11362, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 28210, 0, 3, 22498, 22666, 9742, 25402,
                                                 2588, 2633, 11605, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 28480, 0, 3, 22666, 22834, 9850, 25618,
                                                 2633, 2678, 11740, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 28750, 0, 3, 22834, 23002, 9958, 25834,
                                                 2678, 2723, 11875, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 29020, 0, 3, 23002, 23170, 10066, 26050,
                                                 2723, 2768, 12010, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 29290, 0, 3, 23170, 23338, 10174, 26266,
                                                 2768, 2813, 12145, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 29560, 0, 3, 23338, 23506, 10282, 26482,
                                                 2813, 2858, 12280, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 29830, 0, 3, 23842, 24010, 10714, 26914,
                                                 2948, 2993, 12550, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 30100, 0, 3, 24010, 24178, 10822, 27130,
                                                 2993, 3038, 12685, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 30370, 0, 3, 24178, 24346, 10930, 27346,
                                                 3038, 3083, 12820, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 30640, 0, 3, 24346, 24514, 11038, 27562,
                                                 3083, 3128, 12955, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 30910, 0, 3, 24514, 24682, 11146, 27778,
                                                 3128, 3173, 13090, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 31180, 0, 3, 24682, 24850, 11254, 27994,
                                                 3173, 3218, 13225, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 31450, 0, 3, 25186, 25402, 11605, 28480,
                                                 3308, 3363, 13690, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 31780, 0, 3, 25402, 25618, 11740, 28750,
                                                 3363, 3418, 13855, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 32110, 0, 3, 25618, 25834, 11875, 29020,
                                                 3418, 3473, 14020, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 32440, 0, 3, 25834, 26050, 12010, 29290,
                                                 3473, 3528, 14185, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 32770, 0, 3, 26050, 26266, 12145, 29560,
                                                 3528, 3583, 14350, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 33100, 0, 3, 26698, 26914, 12550, 30100,
                                                 3693, 3748, 14845, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 33430, 0, 3, 26914, 27130, 12685, 30370,
                                                 3748, 3803, 15010, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 33760, 0, 3, 27130, 27346, 12820, 30640,
                                                 3803, 3858, 15175, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 34090, 0, 3, 27346, 27562, 12955, 30910,
                                                 3858, 3913, 15340, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 34420, 0, 3, 27562, 27778, 13090, 31180,
                                                 3913, 3968, 15505, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34750, 3, 4078, 4081, 15676, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34760, 3, 4081, 4084, 15682, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34770, 3, 4084, 4087, 15688, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34780, 3, 4087, 4090, 15694, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34790, 3, 4090, 4093, 15700, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34800, 3, 4093, 4096, 15706, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34810, 3, 4096, 4099, 15712, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34820, 3, 4099, 4102, 15718, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34830, 3, 4102, 4105, 15724, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34840, 3, 4105, 4108, 15730, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34850, 3, 4108, 4111, 15736, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34860, 3, 4117, 4120, 15748, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34870, 3, 4120, 4123, 15754, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34880, 3, 4123, 4126, 15760, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34890, 3, 4126, 4129, 15766, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34900, 3, 4129, 4132, 15772, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34910, 3, 4132, 4135, 15778, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34920, 3, 4135, 4138, 15784, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34930, 3, 4138, 4141, 15790, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34940, 3, 4141, 4144, 15796, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34950, 3, 4144, 4147, 15802, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 34960, 3, 4147, 4150, 15808, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 34970, 0, 3, 15670, 34750, 15832, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35000, 0, 3, 15676, 34760, 15850, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35030, 0, 3, 15682, 34770, 15868, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35060, 0, 3, 15688, 34780, 15886, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35090, 0, 3, 15694, 34790, 15904, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35120, 0, 3, 15700, 34800, 15922, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35150, 0, 3, 15706, 34810, 15940, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35180, 0, 3, 15712, 34820, 15958, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35210, 0, 3, 15718, 34830, 15976, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35240, 0, 3, 15724, 34840, 15994, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35270, 0, 3, 15730, 34850, 16012, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35300, 0, 3, 15742, 34860, 16048, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35330, 0, 3, 15748, 34870, 16066, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35360, 0, 3, 15754, 34880, 16084, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35390, 0, 3, 15760, 34890, 16102, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35420, 0, 3, 15766, 34900, 16120, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35450, 0, 3, 15772, 34910, 16138, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35480, 0, 3, 15778, 34920, 16156, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35510, 0, 3, 15784, 34930, 16174, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35540, 0, 3, 15790, 34940, 16192, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35570, 0, 3, 15796, 34950, 16210, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 35600, 0, 3, 15802, 34960, 16228, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 35630, 0, 3, 15814, 34970, 4390, 4408,
                                                 16282, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 35690, 0, 3, 15832, 35000, 4408, 4426,
                                                 16318, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 35750, 0, 3, 15850, 35030, 4426, 4444,
                                                 16354, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 35810, 0, 3, 15868, 35060, 4444, 4462,
                                                 16390, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 35870, 0, 3, 15886, 35090, 4462, 4480,
                                                 16426, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 35930, 0, 3, 15904, 35120, 4480, 4498,
                                                 16462, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 35990, 0, 3, 15922, 35150, 4498, 4516,
                                                 16498, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36050, 0, 3, 15940, 35180, 4516, 4534,
                                                 16534, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36110, 0, 3, 15958, 35210, 4534, 4552,
                                                 16570, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36170, 0, 3, 15976, 35240, 4552, 4570,
                                                 16606, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36230, 0, 3, 15994, 35270, 4570, 4588,
                                                 16642, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36290, 0, 3, 16030, 35300, 4624, 4642,
                                                 16714, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36350, 0, 3, 16048, 35330, 4642, 4660,
                                                 16750, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36410, 0, 3, 16066, 35360, 4660, 4678,
                                                 16786, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36470, 0, 3, 16084, 35390, 4678, 4696,
                                                 16822, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36530, 0, 3, 16102, 35420, 4696, 4714,
                                                 16858, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36590, 0, 3, 16120, 35450, 4714, 4732,
                                                 16894, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36650, 0, 3, 16138, 35480, 4732, 4750,
                                                 16930, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36710, 0, 3, 16156, 35510, 4750, 4768,
                                                 16966, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36770, 0, 3, 16174, 35540, 4768, 4786,
                                                 17002, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36830, 0, 3, 16192, 35570, 4786, 4804,
                                                 17038, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 36890, 0, 3, 16210, 35600, 4804, 4822,
                                                 17074, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 36950, 0, 3, 16246, 35630, 4858, 4888,
                                                 17110, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37050, 0, 3, 16282, 35690, 4888, 4918,
                                                 17170, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37150, 0, 3, 16318, 35750, 4918, 4948,
                                                 17230, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37250, 0, 3, 16354, 35810, 4948, 4978,
                                                 17290, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37350, 0, 3, 16390, 35870, 4978, 5008,
                                                 17350, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37450, 0, 3, 16426, 35930, 5008, 5038,
                                                 17410, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37550, 0, 3, 16462, 35990, 5038, 5068,
                                                 17470, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37650, 0, 3, 16498, 36050, 5068, 5098,
                                                 17530, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37750, 0, 3, 16534, 36110, 5098, 5128,
                                                 17590, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37850, 0, 3, 16570, 36170, 5128, 5158,
                                                 17650, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 37950, 0, 3, 16606, 36230, 5158, 5188,
                                                 17710, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38050, 0, 3, 16678, 36290, 5248, 5278,
                                                 17770, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38150, 0, 3, 16714, 36350, 5278, 5308,
                                                 17830, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38250, 0, 3, 16750, 36410, 5308, 5338,
                                                 17890, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38350, 0, 3, 16786, 36470, 5338, 5368,
                                                 17950, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38450, 0, 3, 16822, 36530, 5368, 5398,
                                                 18010, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38550, 0, 3, 16858, 36590, 5398, 5428,
                                                 18070, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38650, 0, 3, 16894, 36650, 5428, 5458,
                                                 18130, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38750, 0, 3, 16930, 36710, 5458, 5488,
                                                 18190, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38850, 0, 3, 16966, 36770, 5488, 5518,
                                                 18250, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 38950, 0, 3, 17002, 36830, 5518, 5548,
                                                 18310, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 39050, 0, 3, 17038, 36890, 5548, 5578,
                                                 18370, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 39150, 0, 3, 35630, 35690, 17170, 37150,
                                                 5638, 5683, 18520, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 39300, 0, 3, 35690, 35750, 17230, 37250,
                                                 5683, 5728, 18610, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 39450, 0, 3, 35750, 35810, 17290, 37350,
                                                 5728, 5773, 18700, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 39600, 0, 3, 35810, 35870, 17350, 37450,
                                                 5773, 5818, 18790, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 39750, 0, 3, 35870, 35930, 17410, 37550,
                                                 5818, 5863, 18880, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 39900, 0, 3, 35930, 35990, 17470, 37650,
                                                 5863, 5908, 18970, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 40050, 0, 3, 35990, 36050, 17530, 37750,
                                                 5908, 5953, 19060, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 40200, 0, 3, 36050, 36110, 17590, 37850,
                                                 5953, 5998, 19150, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 40350, 0, 3, 36110, 36170, 17650, 37950,
                                                 5998, 6043, 19240, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 40500, 0, 3, 36290, 36350, 17830, 38250,
                                                 6133, 6178, 19420, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 40650, 0, 3, 36350, 36410, 17890, 38350,
                                                 6178, 6223, 19510, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 40800, 0, 3, 36410, 36470, 17950, 38450,
                                                 6223, 6268, 19600, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 40950, 0, 3, 36470, 36530, 18010, 38550,
                                                 6268, 6313, 19690, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 41100, 0, 3, 36530, 36590, 18070, 38650,
                                                 6313, 6358, 19780, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 41250, 0, 3, 36590, 36650, 18130, 38750,
                                                 6358, 6403, 19870, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 41400, 0, 3, 36650, 36710, 18190, 38850,
                                                 6403, 6448, 19960, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 41550, 0, 3, 36710, 36770, 18250, 38950,
                                                 6448, 6493, 20050, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 41700, 0, 3, 36770, 36830, 18310, 39050,
                                                 6493, 6538, 20140, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 41850, 0, 3, 36950, 37050, 18430, 39150,
                                                 6628, 6691, 20230, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 42060, 0, 3, 37050, 37150, 18520, 39300,
                                                 6691, 6754, 20356, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 42270, 0, 3, 37150, 37250, 18610, 39450,
                                                 6754, 6817, 20482, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 42480, 0, 3, 37250, 37350, 18700, 39600,
                                                 6817, 6880, 20608, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 42690, 0, 3, 37350, 37450, 18790, 39750,
                                                 6880, 6943, 20734, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 42900, 0, 3, 37450, 37550, 18880, 39900,
                                                 6943, 7006, 20860, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 43110, 0, 3, 37550, 37650, 18970, 40050,
                                                 7006, 7069, 20986, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 43320, 0, 3, 37650, 37750, 19060, 40200,
                                                 7069, 7132, 21112, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 43530, 0, 3, 37750, 37850, 19150, 40350,
                                                 7132, 7195, 21238, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 43740, 0, 3, 38050, 38150, 19330, 40500,
                                                 7321, 7384, 21364, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 43950, 0, 3, 38150, 38250, 19420, 40650,
                                                 7384, 7447, 21490, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 44160, 0, 3, 38250, 38350, 19510, 40800,
                                                 7447, 7510, 21616, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 44370, 0, 3, 38350, 38450, 19600, 40950,
                                                 7510, 7573, 21742, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 44580, 0, 3, 38450, 38550, 19690, 41100,
                                                 7573, 7636, 21868, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 44790, 0, 3, 38550, 38650, 19780, 41250,
                                                 7636, 7699, 21994, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 45000, 0, 3, 38650, 38750, 19870, 41400,
                                                 7699, 7762, 22120, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 45210, 0, 3, 38750, 38850, 19960, 41550,
                                                 7762, 7825, 22246, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 45420, 0, 3, 38850, 38950, 20050, 41700,
                                                 7825, 7888, 22372, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 45630, 0, 3, 39150, 39300, 20356, 42270,
                                                 8014, 8098, 22666, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 45910, 0, 3, 39300, 39450, 20482, 42480,
                                                 8098, 8182, 22834, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 46190, 0, 3, 39450, 39600, 20608, 42690,
                                                 8182, 8266, 23002, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 46470, 0, 3, 39600, 39750, 20734, 42900,
                                                 8266, 8350, 23170, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 46750, 0, 3, 39750, 39900, 20860, 43110,
                                                 8350, 8434, 23338, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 47030, 0, 3, 39900, 40050, 20986, 43320,
                                                 8434, 8518, 23506, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 47310, 0, 3, 40050, 40200, 21112, 43530,
                                                 8518, 8602, 23674, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 47590, 0, 3, 40500, 40650, 21490, 44160,
                                                 8770, 8854, 24010, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 47870, 0, 3, 40650, 40800, 21616, 44370,
                                                 8854, 8938, 24178, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 48150, 0, 3, 40800, 40950, 21742, 44580,
                                                 8938, 9022, 24346, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 48430, 0, 3, 40950, 41100, 21868, 44790,
                                                 9022, 9106, 24514, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 48710, 0, 3, 41100, 41250, 21994, 45000,
                                                 9106, 9190, 24682, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 48990, 0, 3, 41250, 41400, 22120, 45210,
                                                 9190, 9274, 24850, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 49270, 0, 3, 41400, 41550, 22246, 45420,
                                                 9274, 9358, 25018, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 49550, 0, 3, 41850, 42060, 22498, 45630,
                                                 9526, 9634, 25186, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 49910, 0, 3, 42060, 42270, 22666, 45910,
                                                 9634, 9742, 25402, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 50270, 0, 3, 42270, 42480, 22834, 46190,
                                                 9742, 9850, 25618, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 50630, 0, 3, 42480, 42690, 23002, 46470,
                                                 9850, 9958, 25834, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 50990, 0, 3, 42690, 42900, 23170, 46750,
                                                 9958, 10066, 26050, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 51350, 0, 3, 42900, 43110, 23338, 47030,
                                                 10066, 10174, 26266, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 51710, 0, 3, 43110, 43320, 23506, 47310,
                                                 10174, 10282, 26482, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 52070, 0, 3, 43740, 43950, 23842, 47590,
                                                 10498, 10606, 26698, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 52430, 0, 3, 43950, 44160, 24010, 47870,
                                                 10606, 10714, 26914, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 52790, 0, 3, 44160, 44370, 24178, 48150,
                                                 10714, 10822, 27130, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 53150, 0, 3, 44370, 44580, 24346, 48430,
                                                 10822, 10930, 27346, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 53510, 0, 3, 44580, 44790, 24514, 48710,
                                                 10930, 11038, 27562, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 53870, 0, 3, 44790, 45000, 24682, 48990,
                                                 11038, 11146, 27778, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 54230, 0, 3, 45000, 45210, 24850, 49270,
                                                 11146, 11254, 27994, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 54590, 0, 3, 45630, 45910, 25402, 50270,
                                                 11470, 11605, 28480, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 55040, 0, 3, 45910, 46190, 25618, 50630,
                                                 11605, 11740, 28750, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 55490, 0, 3, 46190, 46470, 25834, 50990,
                                                 11740, 11875, 29020, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 55940, 0, 3, 46470, 46750, 26050, 51350,
                                                 11875, 12010, 29290, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 56390, 0, 3, 46750, 47030, 26266, 51710,
                                                 12010, 12145, 29560, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 56840, 0, 3, 47590, 47870, 26914, 52790,
                                                 12415, 12550, 30100, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 57290, 0, 3, 47870, 48150, 27130, 53150,
                                                 12550, 12685, 30370, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 57740, 0, 3, 48150, 48430, 27346, 53510,
                                                 12685, 12820, 30640, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 58190, 0, 3, 48430, 48710, 27562, 53870,
                                                 12820, 12955, 30910, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 58640, 0, 3, 48710, 48990, 27778, 54230,
                                                 12955, 13090, 31180, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 59090, 0, 3, 49550, 49910, 28210, 54590,
                                                 13360, 13525, 31450, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 59640, 0, 3, 49910, 50270, 28480, 55040,
                                                 13525, 13690, 31780, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 60190, 0, 3, 50270, 50630, 28750, 55490,
                                                 13690, 13855, 32110, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 60740, 0, 3, 50630, 50990, 29020, 55940,
                                                 13855, 14020, 32440, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 61290, 0, 3, 50990, 51350, 29290, 56390,
                                                 14020, 14185, 32770, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 61840, 0, 3, 52070, 52430, 29830, 56840,
                                                 14515, 14680, 33100, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 62390, 0, 3, 52430, 52790, 30100, 57290,
                                                 14680, 14845, 33430, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 62940, 0, 3, 52790, 53150, 30370, 57740,
                                                 14845, 15010, 33760, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 63490, 0, 3, 53150, 53510, 30640, 58190,
                                                 15010, 15175, 34090, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 64040, 0, 3, 53510, 53870, 30910, 58640,
                                                 15175, 15340, 34420, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64590, 3, 15670, 15676, 34760, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64605, 3, 15676, 15682, 34770, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64620, 3, 15682, 15688, 34780, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64635, 3, 15688, 15694, 34790, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64650, 3, 15694, 15700, 34800, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64665, 3, 15700, 15706, 34810, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64680, 3, 15706, 15712, 34820, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64695, 3, 15712, 15718, 34830, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64710, 3, 15718, 15724, 34840, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64725, 3, 15724, 15730, 34850, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64740, 3, 15742, 15748, 34870, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64755, 3, 15748, 15754, 34880, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64770, 3, 15754, 15760, 34890, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64785, 3, 15760, 15766, 34900, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64800, 3, 15766, 15772, 34910, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64815, 3, 15772, 15778, 34920, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64830, 3, 15778, 15784, 34930, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64845, 3, 15784, 15790, 34940, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64860, 3, 15790, 15796, 34950, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 64875, 3, 15796, 15802, 34960, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 64890, 0, 3, 34750, 64590, 35000, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 64935, 0, 3, 34760, 64605, 35030, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 64980, 0, 3, 34770, 64620, 35060, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65025, 0, 3, 34780, 64635, 35090, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65070, 0, 3, 34790, 64650, 35120, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65115, 0, 3, 34800, 64665, 35150, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65160, 0, 3, 34810, 64680, 35180, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65205, 0, 3, 34820, 64695, 35210, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65250, 0, 3, 34830, 64710, 35240, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65295, 0, 3, 34840, 64725, 35270, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65340, 0, 3, 34860, 64740, 35330, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65385, 0, 3, 34870, 64755, 35360, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65430, 0, 3, 34880, 64770, 35390, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65475, 0, 3, 34890, 64785, 35420, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65520, 0, 3, 34900, 64800, 35450, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65565, 0, 3, 34910, 64815, 35480, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65610, 0, 3, 34920, 64830, 35510, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65655, 0, 3, 34930, 64845, 35540, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65700, 0, 3, 34940, 64860, 35570, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 65745, 0, 3, 34950, 64875, 35600, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 65790, 0, 3, 34970, 64890, 16246, 16282,
                                                 35690, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 65880, 0, 3, 35000, 64935, 16282, 16318,
                                                 35750, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 65970, 0, 3, 35030, 64980, 16318, 16354,
                                                 35810, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66060, 0, 3, 35060, 65025, 16354, 16390,
                                                 35870, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66150, 0, 3, 35090, 65070, 16390, 16426,
                                                 35930, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66240, 0, 3, 35120, 65115, 16426, 16462,
                                                 35990, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66330, 0, 3, 35150, 65160, 16462, 16498,
                                                 36050, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66420, 0, 3, 35180, 65205, 16498, 16534,
                                                 36110, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66510, 0, 3, 35210, 65250, 16534, 16570,
                                                 36170, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66600, 0, 3, 35240, 65295, 16570, 16606,
                                                 36230, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66690, 0, 3, 35300, 65340, 16678, 16714,
                                                 36350, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66780, 0, 3, 35330, 65385, 16714, 16750,
                                                 36410, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66870, 0, 3, 35360, 65430, 16750, 16786,
                                                 36470, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 66960, 0, 3, 35390, 65475, 16786, 16822,
                                                 36530, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 67050, 0, 3, 35420, 65520, 16822, 16858,
                                                 36590, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 67140, 0, 3, 35450, 65565, 16858, 16894,
                                                 36650, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 67230, 0, 3, 35480, 65610, 16894, 16930,
                                                 36710, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 67320, 0, 3, 35510, 65655, 16930, 16966,
                                                 36770, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 67410, 0, 3, 35540, 65700, 16966, 17002,
                                                 36830, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 67500, 0, 3, 35570, 65745, 17002, 17038,
                                                 36890, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 67590, 0, 3, 35690, 65880, 17110, 17170,
                                                 37150, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 67740, 0, 3, 35750, 65970, 17170, 17230,
                                                 37250, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 67890, 0, 3, 35810, 66060, 17230, 17290,
                                                 37350, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 68040, 0, 3, 35870, 66150, 17290, 17350,
                                                 37450, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 68190, 0, 3, 35930, 66240, 17350, 17410,
                                                 37550, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 68340, 0, 3, 35990, 66330, 17410, 17470,
                                                 37650, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 68490, 0, 3, 36050, 66420, 17470, 17530,
                                                 37750, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 68640, 0, 3, 36110, 66510, 17530, 17590,
                                                 37850, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 68790, 0, 3, 36170, 66600, 17590, 17650,
                                                 37950, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 68940, 0, 3, 36350, 66780, 17770, 17830,
                                                 38250, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 69090, 0, 3, 36410, 66870, 17830, 17890,
                                                 38350, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 69240, 0, 3, 36470, 66960, 17890, 17950,
                                                 38450, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 69390, 0, 3, 36530, 67050, 17950, 18010,
                                                 38550, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 69540, 0, 3, 36590, 67140, 18010, 18070,
                                                 38650, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 69690, 0, 3, 36650, 67230, 18070, 18130,
                                                 38750, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 69840, 0, 3, 36710, 67320, 18130, 18190,
                                                 38850, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 69990, 0, 3, 36770, 67410, 18190, 18250,
                                                 38950, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 70140, 0, 3, 36830, 67500, 18250, 18310,
                                                 39050, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 70290, 0, 3, 65790, 65880, 37150, 67740,
                                                 18430, 18520, 39300, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 70515, 0, 3, 65880, 65970, 37250, 67890,
                                                 18520, 18610, 39450, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 70740, 0, 3, 65970, 66060, 37350, 68040,
                                                 18610, 18700, 39600, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 70965, 0, 3, 66060, 66150, 37450, 68190,
                                                 18700, 18790, 39750, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 71190, 0, 3, 66150, 66240, 37550, 68340,
                                                 18790, 18880, 39900, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 71415, 0, 3, 66240, 66330, 37650, 68490,
                                                 18880, 18970, 40050, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 71640, 0, 3, 66330, 66420, 37750, 68640,
                                                 18970, 19060, 40200, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 71865, 0, 3, 66420, 66510, 37850, 68790,
                                                 19060, 19150, 40350, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 72090, 0, 3, 66690, 66780, 38250, 69090,
                                                 19330, 19420, 40650, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 72315, 0, 3, 66780, 66870, 38350, 69240,
                                                 19420, 19510, 40800, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 72540, 0, 3, 66870, 66960, 38450, 69390,
                                                 19510, 19600, 40950, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 72765, 0, 3, 66960, 67050, 38550, 69540,
                                                 19600, 19690, 41100, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 72990, 0, 3, 67050, 67140, 38650, 69690,
                                                 19690, 19780, 41250, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 73215, 0, 3, 67140, 67230, 38750, 69840,
                                                 19780, 19870, 41400, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 73440, 0, 3, 67230, 67320, 38850, 69990,
                                                 19870, 19960, 41550, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 73665, 0, 3, 67320, 67410, 38950, 70140,
                                                 19960, 20050, 41700, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 73890, 0, 3, 67590, 67740, 39300, 70515,
                                                 20230, 20356, 42270, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 74205, 0, 3, 67740, 67890, 39450, 70740,
                                                 20356, 20482, 42480, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 74520, 0, 3, 67890, 68040, 39600, 70965,
                                                 20482, 20608, 42690, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 74835, 0, 3, 68040, 68190, 39750, 71190,
                                                 20608, 20734, 42900, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 75150, 0, 3, 68190, 68340, 39900, 71415,
                                                 20734, 20860, 43110, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 75465, 0, 3, 68340, 68490, 40050, 71640,
                                                 20860, 20986, 43320, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 75780, 0, 3, 68490, 68640, 40200, 71865,
                                                 20986, 21112, 43530, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 76095, 0, 3, 68940, 69090, 40650, 72315,
                                                 21364, 21490, 44160, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 76410, 0, 3, 69090, 69240, 40800, 72540,
                                                 21490, 21616, 44370, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 76725, 0, 3, 69240, 69390, 40950, 72765,
                                                 21616, 21742, 44580, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 77040, 0, 3, 69390, 69540, 41100, 72990,
                                                 21742, 21868, 44790, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 77355, 0, 3, 69540, 69690, 41250, 73215,
                                                 21868, 21994, 45000, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 77670, 0, 3, 69690, 69840, 41400, 73440,
                                                 21994, 22120, 45210, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 77985, 0, 3, 69840, 69990, 41550, 73665,
                                                 22120, 22246, 45420, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 78300, 0, 3, 70290, 70515, 42270, 74205,
                                                 22498, 22666, 45910, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 78720, 0, 3, 70515, 70740, 42480, 74520,
                                                 22666, 22834, 46190, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 79140, 0, 3, 70740, 70965, 42690, 74835,
                                                 22834, 23002, 46470, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 79560, 0, 3, 70965, 71190, 42900, 75150,
                                                 23002, 23170, 46750, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 79980, 0, 3, 71190, 71415, 43110, 75465,
                                                 23170, 23338, 47030, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 80400, 0, 3, 71415, 71640, 43320, 75780,
                                                 23338, 23506, 47310, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 80820, 0, 3, 72090, 72315, 44160, 76410,
                                                 23842, 24010, 47870, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 81240, 0, 3, 72315, 72540, 44370, 76725,
                                                 24010, 24178, 48150, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 81660, 0, 3, 72540, 72765, 44580, 77040,
                                                 24178, 24346, 48430, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 82080, 0, 3, 72765, 72990, 44790, 77355,
                                                 24346, 24514, 48710, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 82500, 0, 3, 72990, 73215, 45000, 77670,
                                                 24514, 24682, 48990, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 82920, 0, 3, 73215, 73440, 45210, 77985,
                                                 24682, 24850, 49270, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 83340, 0, 3, 73890, 74205, 45910, 78720,
                                                 25186, 25402, 50270, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 83880, 0, 3, 74205, 74520, 46190, 79140,
                                                 25402, 25618, 50630, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 84420, 0, 3, 74520, 74835, 46470, 79560,
                                                 25618, 25834, 50990, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 84960, 0, 3, 74835, 75150, 46750, 79980,
                                                 25834, 26050, 51350, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 85500, 0, 3, 75150, 75465, 47030, 80400,
                                                 26050, 26266, 51710, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 86040, 0, 3, 76095, 76410, 47870, 81240,
                                                 26698, 26914, 52790, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 86580, 0, 3, 76410, 76725, 48150, 81660,
                                                 26914, 27130, 53150, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 87120, 0, 3, 76725, 77040, 48430, 82080,
                                                 27130, 27346, 53510, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 87660, 0, 3, 77040, 77355, 48710, 82500,
                                                 27346, 27562, 53870, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 88200, 0, 3, 77355, 77670, 48990, 82920,
                                                 27562, 27778, 54230, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 88740, 0, 3, 78300, 78720, 50270, 83880,
                                                 28210, 28480, 55040, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 89415, 0, 3, 78720, 79140, 50630, 84420,
                                                 28480, 28750, 55490, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 90090, 0, 3, 79140, 79560, 50990, 84960,
                                                 28750, 29020, 55940, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 90765, 0, 3, 79560, 79980, 51350, 85500,
                                                 29020, 29290, 56390, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 91440, 0, 3, 80820, 81240, 52790, 86580,
                                                 29830, 30100, 57290, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 92115, 0, 3, 81240, 81660, 53150, 87120,
                                                 30100, 30370, 57740, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 92790, 0, 3, 81660, 82080, 53510, 87660,
                                                 30370, 30640, 58190, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 93465, 0, 3, 82080, 82500, 53870, 88200,
                                                 30640, 30910, 58640, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 94140, 0, 3, 83340, 83880, 55040, 89415,
                                                 31450, 31780, 60190, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 94965, 0, 3, 83880, 84420, 55490, 90090,
                                                 31780, 32110, 60740, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 95790, 0, 3, 84420, 84960, 55940, 90765,
                                                 32110, 32440, 61290, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 96615, 0, 3, 86040, 86580, 57290, 92115,
                                                 33100, 33430, 62940, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 97440, 0, 3, 86580, 87120, 57740, 92790,
                                                 33430, 33760, 63490, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 98265, 0, 3, 87120, 87660, 58190, 93465,
                                                 33760, 34090, 64040, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99090, 3, 34750, 34760, 64605, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99111, 3, 34760, 34770, 64620, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99132, 3, 34770, 34780, 64635, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99153, 3, 34780, 34790, 64650, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99174, 3, 34790, 34800, 64665, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99195, 3, 34800, 34810, 64680, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99216, 3, 34810, 34820, 64695, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99237, 3, 34820, 34830, 64710, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99258, 3, 34830, 34840, 64725, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99279, 3, 34860, 34870, 64755, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99300, 3, 34870, 34880, 64770, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99321, 3, 34880, 34890, 64785, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99342, 3, 34890, 34900, 64800, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99363, 3, 34900, 34910, 64815, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99384, 3, 34910, 34920, 64830, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99405, 3, 34920, 34930, 64845, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99426, 3, 34930, 34940, 64860, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 99447, 3, 34940, 34950, 64875, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 99468, 0, 3, 64590, 99090, 64935, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 99531, 0, 3, 64605, 99111, 64980, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 99594, 0, 3, 64620, 99132, 65025, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 99657, 0, 3, 64635, 99153, 65070, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 99720, 0, 3, 64650, 99174, 65115, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 99783, 0, 3, 64665, 99195, 65160, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 99846, 0, 3, 64680, 99216, 65205, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 99909, 0, 3, 64695, 99237, 65250, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 99972, 0, 3, 64710, 99258, 65295, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 100035, 0, 3, 64740, 99279, 65385,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 100098, 0, 3, 64755, 99300, 65430,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 100161, 0, 3, 64770, 99321, 65475,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 100224, 0, 3, 64785, 99342, 65520,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 100287, 0, 3, 64800, 99363, 65565,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 100350, 0, 3, 64815, 99384, 65610,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 100413, 0, 3, 64830, 99405, 65655,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 100476, 0, 3, 64845, 99426, 65700,
                                                 ncols, p);

            compute_prim_ph_electron_repulsion_0(buffer, 100539, 0, 3, 64860, 99447, 65745,
                                                 ncols, p);

            compute_prim_dh_electron_repulsion_0(buffer, 100602, 0, 3, 64890, 99468, 35630,
                                                 35690, 65880, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 100728, 0, 3, 64935, 99531, 35690,
                                                 35750, 65970, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 100854, 0, 3, 64980, 99594, 35750,
                                                 35810, 66060, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 100980, 0, 3, 65025, 99657, 35810,
                                                 35870, 66150, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 101106, 0, 3, 65070, 99720, 35870,
                                                 35930, 66240, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 101232, 0, 3, 65115, 99783, 35930,
                                                 35990, 66330, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 101358, 0, 3, 65160, 99846, 35990,
                                                 36050, 66420, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 101484, 0, 3, 65205, 99909, 36050,
                                                 36110, 66510, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 101610, 0, 3, 65250, 99972, 36110,
                                                 36170, 66600, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 101736, 0, 3, 65340, 100035, 36290,
                                                 36350, 66780, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 101862, 0, 3, 65385, 100098, 36350,
                                                 36410, 66870, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 101988, 0, 3, 65430, 100161, 36410,
                                                 36470, 66960, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 102114, 0, 3, 65475, 100224, 36470,
                                                 36530, 67050, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 102240, 0, 3, 65520, 100287, 36530,
                                                 36590, 67140, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 102366, 0, 3, 65565, 100350, 36590,
                                                 36650, 67230, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 102492, 0, 3, 65610, 100413, 36650,
                                                 36710, 67320, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 102618, 0, 3, 65655, 100476, 36710,
                                                 36770, 67410, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 102744, 0, 3, 65700, 100539, 36770,
                                                 36830, 67500, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 102870, 0, 3, 65790, 100602, 36950,
                                                 37050, 67590, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 103080, 0, 3, 65880, 100728, 37050,
                                                 37150, 67740, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 103290, 0, 3, 65970, 100854, 37150,
                                                 37250, 67890, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 103500, 0, 3, 66060, 100980, 37250,
                                                 37350, 68040, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 103710, 0, 3, 66150, 101106, 37350,
                                                 37450, 68190, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 103920, 0, 3, 66240, 101232, 37450,
                                                 37550, 68340, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 104130, 0, 3, 66330, 101358, 37550,
                                                 37650, 68490, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 104340, 0, 3, 66420, 101484, 37650,
                                                 37750, 68640, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 104550, 0, 3, 66510, 101610, 37750,
                                                 37850, 68790, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 104760, 0, 3, 66690, 101736, 38050,
                                                 38150, 68940, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 104970, 0, 3, 66780, 101862, 38150,
                                                 38250, 69090, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 105180, 0, 3, 66870, 101988, 38250,
                                                 38350, 69240, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 105390, 0, 3, 66960, 102114, 38350,
                                                 38450, 69390, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 105600, 0, 3, 67050, 102240, 38450,
                                                 38550, 69540, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 105810, 0, 3, 67140, 102366, 38550,
                                                 38650, 69690, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 106020, 0, 3, 67230, 102492, 38650,
                                                 38750, 69840, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 106230, 0, 3, 67320, 102618, 38750,
                                                 38850, 69990, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 106440, 0, 3, 67410, 102744, 38850,
                                                 38950, 70140, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 106650, 0, 3, 100602, 100728, 67740,
                                                 103290, 39150, 39300, 70515, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 106965, 0, 3, 100728, 100854, 67890,
                                                 103500, 39300, 39450, 70740, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 107280, 0, 3, 100854, 100980, 68040,
                                                 103710, 39450, 39600, 70965, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 107595, 0, 3, 100980, 101106, 68190,
                                                 103920, 39600, 39750, 71190, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 107910, 0, 3, 101106, 101232, 68340,
                                                 104130, 39750, 39900, 71415, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 108225, 0, 3, 101232, 101358, 68490,
                                                 104340, 39900, 40050, 71640, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 108540, 0, 3, 101358, 101484, 68640,
                                                 104550, 40050, 40200, 71865, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 108855, 0, 3, 101736, 101862, 69090,
                                                 105180, 40500, 40650, 72315, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 109170, 0, 3, 101862, 101988, 69240,
                                                 105390, 40650, 40800, 72540, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 109485, 0, 3, 101988, 102114, 69390,
                                                 105600, 40800, 40950, 72765, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 109800, 0, 3, 102114, 102240, 69540,
                                                 105810, 40950, 41100, 72990, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 110115, 0, 3, 102240, 102366, 69690,
                                                 106020, 41100, 41250, 73215, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 110430, 0, 3, 102366, 102492, 69840,
                                                 106230, 41250, 41400, 73440, ncols, alpha, beta,
                                                 p);

            compute_prim_gh_electron_repulsion_0(buffer, 110745, 0, 3, 102492, 102618, 69990,
                                                 106440, 41400, 41550, 73665, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 111060, 0, 3, 102870, 103080, 70290,
                                                 106650, 41850, 42060, 73890, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 111501, 0, 3, 103080, 103290, 70515,
                                                 106965, 42060, 42270, 74205, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 111942, 0, 3, 103290, 103500, 70740,
                                                 107280, 42270, 42480, 74520, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 112383, 0, 3, 103500, 103710, 70965,
                                                 107595, 42480, 42690, 74835, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 112824, 0, 3, 103710, 103920, 71190,
                                                 107910, 42690, 42900, 75150, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 113265, 0, 3, 103920, 104130, 71415,
                                                 108225, 42900, 43110, 75465, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 113706, 0, 3, 104130, 104340, 71640,
                                                 108540, 43110, 43320, 75780, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 114147, 0, 3, 104760, 104970, 72090,
                                                 108855, 43740, 43950, 76095, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 114588, 0, 3, 104970, 105180, 72315,
                                                 109170, 43950, 44160, 76410, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 115029, 0, 3, 105180, 105390, 72540,
                                                 109485, 44160, 44370, 76725, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 115470, 0, 3, 105390, 105600, 72765,
                                                 109800, 44370, 44580, 77040, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 115911, 0, 3, 105600, 105810, 72990,
                                                 110115, 44580, 44790, 77355, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 116352, 0, 3, 105810, 106020, 73215,
                                                 110430, 44790, 45000, 77670, ncols, alpha, beta,
                                                 p);

            compute_prim_hh_electron_repulsion_0(buffer, 116793, 0, 3, 106020, 106230, 73440,
                                                 110745, 45000, 45210, 77985, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 117234, 0, 3, 106650, 106965, 74205,
                                                 111942, 45630, 45910, 78720, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 117822, 0, 3, 106965, 107280, 74520,
                                                 112383, 45910, 46190, 79140, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 118410, 0, 3, 107280, 107595, 74835,
                                                 112824, 46190, 46470, 79560, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 118998, 0, 3, 107595, 107910, 75150,
                                                 113265, 46470, 46750, 79980, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 119586, 0, 3, 107910, 108225, 75465,
                                                 113706, 46750, 47030, 80400, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 120174, 0, 3, 108855, 109170, 76410,
                                                 115029, 47590, 47870, 81240, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 120762, 0, 3, 109170, 109485, 76725,
                                                 115470, 47870, 48150, 81660, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 121350, 0, 3, 109485, 109800, 77040,
                                                 115911, 48150, 48430, 82080, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 121938, 0, 3, 109800, 110115, 77355,
                                                 116352, 48430, 48710, 82500, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 122526, 0, 3, 110115, 110430, 77670,
                                                 116793, 48710, 48990, 82920, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 123114, 0, 3, 111060, 111501, 78300,
                                                 117234, 49550, 49910, 83340, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 123870, 0, 3, 111501, 111942, 78720,
                                                 117822, 49910, 50270, 83880, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 124626, 0, 3, 111942, 112383, 79140,
                                                 118410, 50270, 50630, 84420, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 125382, 0, 3, 112383, 112824, 79560,
                                                 118998, 50630, 50990, 84960, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 126138, 0, 3, 112824, 113265, 79980,
                                                 119586, 50990, 51350, 85500, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 126894, 0, 3, 114147, 114588, 80820,
                                                 120174, 52070, 52430, 86040, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 127650, 0, 3, 114588, 115029, 81240,
                                                 120762, 52430, 52790, 86580, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 128406, 0, 3, 115029, 115470, 81660,
                                                 121350, 52790, 53150, 87120, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 129162, 0, 3, 115470, 115911, 82080,
                                                 121938, 53150, 53510, 87660, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 129918, 0, 3, 115911, 116352, 82500,
                                                 122526, 53510, 53870, 88200, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 130674, 0, 3, 117234, 117822, 83880,
                                                 124626, 54590, 55040, 89415, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 131619, 0, 3, 117822, 118410, 84420,
                                                 125382, 55040, 55490, 90090, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 132564, 0, 3, 118410, 118998, 84960,
                                                 126138, 55490, 55940, 90765, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 133509, 0, 3, 120174, 120762, 86580,
                                                 128406, 56840, 57290, 92115, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 134454, 0, 3, 120762, 121350, 87120,
                                                 129162, 57290, 57740, 92790, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 135399, 0, 3, 121350, 121938, 87660,
                                                 129918, 57740, 58190, 93465, ncols, alpha, beta,
                                                 p);

            compute_prim_mh_electron_repulsion_0(buffer, 136344, 0, 3, 123114, 123870, 88740,
                                                 130674, 59090, 59640, 94140, ncols, alpha, beta,
                                                 p);

            compute_prim_mh_electron_repulsion_0(buffer, 137499, 0, 3, 123870, 124626, 89415,
                                                 131619, 59640, 60190, 94965, ncols, alpha, beta,
                                                 p);

            compute_prim_mh_electron_repulsion_0(buffer, 138654, 0, 3, 124626, 125382, 90090,
                                                 132564, 60190, 60740, 95790, ncols, alpha, beta,
                                                 p);

            compute_prim_mh_electron_repulsion_0(buffer, 139809, 0, 3, 126894, 127650, 91440,
                                                 133509, 61840, 62390, 96615, ncols, alpha, beta,
                                                 p);

            compute_prim_mh_electron_repulsion_0(buffer, 140964, 0, 3, 127650, 128406, 92115,
                                                 134454, 62390, 62940, 97440, ncols, alpha, beta,
                                                 p);

            compute_prim_mh_electron_repulsion_0(buffer, 142119, 0, 3, 128406, 129162, 92790,
                                                 135399, 62940, 63490, 98265, ncols, alpha, beta,
                                                 p);

            compute_prim_si_electron_repulsion_0(buffer, 143274, 3, 64590, 64605, 99111, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143302, 3, 64605, 64620, 99132, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143330, 3, 64620, 64635, 99153, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143358, 3, 64635, 64650, 99174, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143386, 3, 64650, 64665, 99195, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143414, 3, 64665, 64680, 99216, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143442, 3, 64680, 64695, 99237, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143470, 3, 64695, 64710, 99258, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143498, 3, 64740, 64755, 99300, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143526, 3, 64755, 64770, 99321, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143554, 3, 64770, 64785, 99342, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143582, 3, 64785, 64800, 99363, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143610, 3, 64800, 64815, 99384, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143638, 3, 64815, 64830, 99405, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143666, 3, 64830, 64845, 99426, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 143694, 3, 64845, 64860, 99447, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 143722, 0, 3, 99090, 143274, 99531,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 143806, 0, 3, 99111, 143302, 99594,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 143890, 0, 3, 99132, 143330, 99657,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 143974, 0, 3, 99153, 143358, 99720,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144058, 0, 3, 99174, 143386, 99783,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144142, 0, 3, 99195, 143414, 99846,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144226, 0, 3, 99216, 143442, 99909,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144310, 0, 3, 99237, 143470, 99972,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144394, 0, 3, 99279, 143498, 100098,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144478, 0, 3, 99300, 143526, 100161,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144562, 0, 3, 99321, 143554, 100224,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144646, 0, 3, 99342, 143582, 100287,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144730, 0, 3, 99363, 143610, 100350,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144814, 0, 3, 99384, 143638, 100413,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144898, 0, 3, 99405, 143666, 100476,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 144982, 0, 3, 99426, 143694, 100539,
                                                 ncols, p);

            compute_prim_di_electron_repulsion_0(buffer, 145066, 0, 3, 99468, 143722, 65790,
                                                 65880, 100728, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 145234, 0, 3, 99531, 143806, 65880,
                                                 65970, 100854, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 145402, 0, 3, 99594, 143890, 65970,
                                                 66060, 100980, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 145570, 0, 3, 99657, 143974, 66060,
                                                 66150, 101106, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 145738, 0, 3, 99720, 144058, 66150,
                                                 66240, 101232, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 145906, 0, 3, 99783, 144142, 66240,
                                                 66330, 101358, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 146074, 0, 3, 99846, 144226, 66330,
                                                 66420, 101484, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 146242, 0, 3, 99909, 144310, 66420,
                                                 66510, 101610, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 146410, 0, 3, 100035, 144394, 66690,
                                                 66780, 101862, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 146578, 0, 3, 100098, 144478, 66780,
                                                 66870, 101988, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 146746, 0, 3, 100161, 144562, 66870,
                                                 66960, 102114, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 146914, 0, 3, 100224, 144646, 66960,
                                                 67050, 102240, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 147082, 0, 3, 100287, 144730, 67050,
                                                 67140, 102366, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 147250, 0, 3, 100350, 144814, 67140,
                                                 67230, 102492, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 147418, 0, 3, 100413, 144898, 67230,
                                                 67320, 102618, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 147586, 0, 3, 100476, 144982, 67320,
                                                 67410, 102744, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 147754, 0, 3, 100728, 145234, 67590,
                                                 67740, 103290, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 148034, 0, 3, 100854, 145402, 67740,
                                                 67890, 103500, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 148314, 0, 3, 100980, 145570, 67890,
                                                 68040, 103710, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 148594, 0, 3, 101106, 145738, 68040,
                                                 68190, 103920, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 148874, 0, 3, 101232, 145906, 68190,
                                                 68340, 104130, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 149154, 0, 3, 101358, 146074, 68340,
                                                 68490, 104340, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 149434, 0, 3, 101484, 146242, 68490,
                                                 68640, 104550, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 149714, 0, 3, 101862, 146578, 68940,
                                                 69090, 105180, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 149994, 0, 3, 101988, 146746, 69090,
                                                 69240, 105390, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 150274, 0, 3, 102114, 146914, 69240,
                                                 69390, 105600, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 150554, 0, 3, 102240, 147082, 69390,
                                                 69540, 105810, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 150834, 0, 3, 102366, 147250, 69540,
                                                 69690, 106020, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 151114, 0, 3, 102492, 147418, 69690,
                                                 69840, 106230, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 151394, 0, 3, 102618, 147586, 69840,
                                                 69990, 106440, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 151674, 0, 3, 145066, 145234, 103290,
                                                 148034, 70290, 70515, 106965, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 152094, 0, 3, 145234, 145402, 103500,
                                                 148314, 70515, 70740, 107280, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 152514, 0, 3, 145402, 145570, 103710,
                                                 148594, 70740, 70965, 107595, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 152934, 0, 3, 145570, 145738, 103920,
                                                 148874, 70965, 71190, 107910, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 153354, 0, 3, 145738, 145906, 104130,
                                                 149154, 71190, 71415, 108225, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 153774, 0, 3, 145906, 146074, 104340,
                                                 149434, 71415, 71640, 108540, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 154194, 0, 3, 146410, 146578, 105180,
                                                 149994, 72090, 72315, 109170, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 154614, 0, 3, 146578, 146746, 105390,
                                                 150274, 72315, 72540, 109485, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 155034, 0, 3, 146746, 146914, 105600,
                                                 150554, 72540, 72765, 109800, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 155454, 0, 3, 146914, 147082, 105810,
                                                 150834, 72765, 72990, 110115, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 155874, 0, 3, 147082, 147250, 106020,
                                                 151114, 72990, 73215, 110430, ncols, alpha,
                                                 beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 156294, 0, 3, 147250, 147418, 106230,
                                                 151394, 73215, 73440, 110745, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 156714, 0, 3, 147754, 148034, 106965,
                                                 152094, 73890, 74205, 111942, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 157302, 0, 3, 148034, 148314, 107280,
                                                 152514, 74205, 74520, 112383, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 157890, 0, 3, 148314, 148594, 107595,
                                                 152934, 74520, 74835, 112824, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 158478, 0, 3, 148594, 148874, 107910,
                                                 153354, 74835, 75150, 113265, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 159066, 0, 3, 148874, 149154, 108225,
                                                 153774, 75150, 75465, 113706, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 159654, 0, 3, 149714, 149994, 109170,
                                                 154614, 76095, 76410, 115029, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 160242, 0, 3, 149994, 150274, 109485,
                                                 155034, 76410, 76725, 115470, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 160830, 0, 3, 150274, 150554, 109800,
                                                 155454, 76725, 77040, 115911, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 161418, 0, 3, 150554, 150834, 110115,
                                                 155874, 77040, 77355, 116352, ncols, alpha,
                                                 beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 162006, 0, 3, 150834, 151114, 110430,
                                                 156294, 77355, 77670, 116793, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 162594, 0, 3, 151674, 152094, 111942,
                                                 157302, 78300, 78720, 117822, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 163378, 0, 3, 152094, 152514, 112383,
                                                 157890, 78720, 79140, 118410, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 164162, 0, 3, 152514, 152934, 112824,
                                                 158478, 79140, 79560, 118998, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 164946, 0, 3, 152934, 153354, 113265,
                                                 159066, 79560, 79980, 119586, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 165730, 0, 3, 154194, 154614, 115029,
                                                 160242, 80820, 81240, 120762, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 166514, 0, 3, 154614, 155034, 115470,
                                                 160830, 81240, 81660, 121350, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 167298, 0, 3, 155034, 155454, 115911,
                                                 161418, 81660, 82080, 121938, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 168082, 0, 3, 155454, 155874, 116352,
                                                 162006, 82080, 82500, 122526, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 168866, 0, 3, 156714, 157302, 117822,
                                                 163378, 83340, 83880, 124626, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 169874, 0, 3, 157302, 157890, 118410,
                                                 164162, 83880, 84420, 125382, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 170882, 0, 3, 157890, 158478, 118998,
                                                 164946, 84420, 84960, 126138, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 171890, 0, 3, 159654, 160242, 120762,
                                                 166514, 86040, 86580, 128406, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 172898, 0, 3, 160242, 160830, 121350,
                                                 167298, 86580, 87120, 129162, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 173906, 0, 3, 160830, 161418, 121938,
                                                 168082, 87120, 87660, 129918, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 174914, 0, 3, 162594, 163378, 124626,
                                                 169874, 88740, 89415, 131619, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 176174, 0, 3, 163378, 164162, 125382,
                                                 170882, 89415, 90090, 132564, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 177434, 0, 3, 165730, 166514, 128406,
                                                 172898, 91440, 92115, 134454, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 178694, 0, 3, 166514, 167298, 129162,
                                                 173906, 92115, 92790, 135399, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 179954, 0, 3, 168866, 169874, 131619,
                                                 176174, 94140, 94965, 138654, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 181494, 0, 3, 171890, 172898, 134454,
                                                 178694, 96615, 97440, 142119, ncols, alpha,
                                                 beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183034, 3, 99090, 99111, 143302, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183070, 3, 99111, 99132, 143330, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183106, 3, 99132, 99153, 143358, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183142, 3, 99153, 99174, 143386, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183178, 3, 99174, 99195, 143414, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183214, 3, 99195, 99216, 143442, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183250, 3, 99216, 99237, 143470, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183286, 3, 99279, 99300, 143526, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183322, 3, 99300, 99321, 143554, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183358, 3, 99321, 99342, 143582, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183394, 3, 99342, 99363, 143610, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183430, 3, 99363, 99384, 143638, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183466, 3, 99384, 99405, 143666, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 183502, 3, 99405, 99426, 143694, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 183538, 0, 3, 143274, 183034, 143806,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 183646, 0, 3, 143302, 183070, 143890,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 183754, 0, 3, 143330, 183106, 143974,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 183862, 0, 3, 143358, 183142, 144058,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 183970, 0, 3, 143386, 183178, 144142,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184078, 0, 3, 143414, 183214, 144226,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184186, 0, 3, 143442, 183250, 144310,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184294, 0, 3, 143498, 183286, 144478,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184402, 0, 3, 143526, 183322, 144562,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184510, 0, 3, 143554, 183358, 144646,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184618, 0, 3, 143582, 183394, 144730,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184726, 0, 3, 143610, 183430, 144814,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184834, 0, 3, 143638, 183466, 144898,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 184942, 0, 3, 143666, 183502, 144982,
                                                 ncols, p);

            compute_prim_dk_electron_repulsion_0(buffer, 185050, 0, 3, 143722, 183538, 100602,
                                                 100728, 145234, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 185266, 0, 3, 143806, 183646, 100728,
                                                 100854, 145402, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 185482, 0, 3, 143890, 183754, 100854,
                                                 100980, 145570, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 185698, 0, 3, 143974, 183862, 100980,
                                                 101106, 145738, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 185914, 0, 3, 144058, 183970, 101106,
                                                 101232, 145906, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 186130, 0, 3, 144142, 184078, 101232,
                                                 101358, 146074, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 186346, 0, 3, 144226, 184186, 101358,
                                                 101484, 146242, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 186562, 0, 3, 144394, 184294, 101736,
                                                 101862, 146578, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 186778, 0, 3, 144478, 184402, 101862,
                                                 101988, 146746, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 186994, 0, 3, 144562, 184510, 101988,
                                                 102114, 146914, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 187210, 0, 3, 144646, 184618, 102114,
                                                 102240, 147082, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 187426, 0, 3, 144730, 184726, 102240,
                                                 102366, 147250, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 187642, 0, 3, 144814, 184834, 102366,
                                                 102492, 147418, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 187858, 0, 3, 144898, 184942, 102492,
                                                 102618, 147586, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 188074, 0, 3, 145066, 185050, 102870,
                                                 103080, 147754, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 188434, 0, 3, 145234, 185266, 103080,
                                                 103290, 148034, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 188794, 0, 3, 145402, 185482, 103290,
                                                 103500, 148314, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 189154, 0, 3, 145570, 185698, 103500,
                                                 103710, 148594, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 189514, 0, 3, 145738, 185914, 103710,
                                                 103920, 148874, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 189874, 0, 3, 145906, 186130, 103920,
                                                 104130, 149154, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 190234, 0, 3, 146074, 186346, 104130,
                                                 104340, 149434, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 190594, 0, 3, 146410, 186562, 104760,
                                                 104970, 149714, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 190954, 0, 3, 146578, 186778, 104970,
                                                 105180, 149994, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 191314, 0, 3, 146746, 186994, 105180,
                                                 105390, 150274, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 191674, 0, 3, 146914, 187210, 105390,
                                                 105600, 150554, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 192034, 0, 3, 147082, 187426, 105600,
                                                 105810, 150834, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 192394, 0, 3, 147250, 187642, 105810,
                                                 106020, 151114, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 192754, 0, 3, 147418, 187858, 106020,
                                                 106230, 151394, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 193114, 0, 3, 185050, 185266, 148034,
                                                 188794, 106650, 106965, 152094, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 193654, 0, 3, 185266, 185482, 148314,
                                                 189154, 106965, 107280, 152514, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 194194, 0, 3, 185482, 185698, 148594,
                                                 189514, 107280, 107595, 152934, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 194734, 0, 3, 185698, 185914, 148874,
                                                 189874, 107595, 107910, 153354, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 195274, 0, 3, 185914, 186130, 149154,
                                                 190234, 107910, 108225, 153774, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 195814, 0, 3, 186562, 186778, 149994,
                                                 191314, 108855, 109170, 154614, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 196354, 0, 3, 186778, 186994, 150274,
                                                 191674, 109170, 109485, 155034, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 196894, 0, 3, 186994, 187210, 150554,
                                                 192034, 109485, 109800, 155454, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 197434, 0, 3, 187210, 187426, 150834,
                                                 192394, 109800, 110115, 155874, ncols, alpha,
                                                 beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 197974, 0, 3, 187426, 187642, 151114,
                                                 192754, 110115, 110430, 156294, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 198514, 0, 3, 188074, 188434, 151674,
                                                 193114, 111060, 111501, 156714, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 199270, 0, 3, 188434, 188794, 152094,
                                                 193654, 111501, 111942, 157302, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 200026, 0, 3, 188794, 189154, 152514,
                                                 194194, 111942, 112383, 157890, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 200782, 0, 3, 189154, 189514, 152934,
                                                 194734, 112383, 112824, 158478, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 201538, 0, 3, 189514, 189874, 153354,
                                                 195274, 112824, 113265, 159066, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 202294, 0, 3, 190594, 190954, 154194,
                                                 195814, 114147, 114588, 159654, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 203050, 0, 3, 190954, 191314, 154614,
                                                 196354, 114588, 115029, 160242, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 203806, 0, 3, 191314, 191674, 155034,
                                                 196894, 115029, 115470, 160830, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 204562, 0, 3, 191674, 192034, 155454,
                                                 197434, 115470, 115911, 161418, ncols, alpha,
                                                 beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 205318, 0, 3, 192034, 192394, 155874,
                                                 197974, 115911, 116352, 162006, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 206074, 0, 3, 193114, 193654, 157302,
                                                 200026, 117234, 117822, 163378, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 207082, 0, 3, 193654, 194194, 157890,
                                                 200782, 117822, 118410, 164162, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 208090, 0, 3, 194194, 194734, 158478,
                                                 201538, 118410, 118998, 164946, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 209098, 0, 3, 195814, 196354, 160242,
                                                 203806, 120174, 120762, 166514, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 210106, 0, 3, 196354, 196894, 160830,
                                                 204562, 120762, 121350, 167298, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 211114, 0, 3, 196894, 197434, 161418,
                                                 205318, 121350, 121938, 168082, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 212122, 0, 3, 198514, 199270, 162594,
                                                 206074, 123114, 123870, 168866, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 213418, 0, 3, 199270, 200026, 163378,
                                                 207082, 123870, 124626, 169874, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 214714, 0, 3, 200026, 200782, 164162,
                                                 208090, 124626, 125382, 170882, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 216010, 0, 3, 202294, 203050, 165730,
                                                 209098, 126894, 127650, 171890, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 217306, 0, 3, 203050, 203806, 166514,
                                                 210106, 127650, 128406, 172898, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 218602, 0, 3, 203806, 204562, 167298,
                                                 211114, 128406, 129162, 173906, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 219898, 0, 3, 206074, 207082, 169874,
                                                 214714, 130674, 131619, 176174, ncols, alpha,
                                                 beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 221518, 0, 3, 209098, 210106, 172898,
                                                 218602, 133509, 134454, 178694, ncols, alpha,
                                                 beta, p);

            compute_prim_mk_electron_repulsion_0(buffer, 223138, 0, 3, 212122, 213418, 174914,
                                                 219898, 136344, 137499, 179954, ncols, alpha,
                                                 beta, p);

            compute_prim_mk_electron_repulsion_0(buffer, 225118, 0, 3, 216010, 217306, 177434,
                                                 221518, 139809, 140964, 181494, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_l_x(buffer, 227098, 216010, 225118, 1, 36, ncols, alpha);

            simdgeo::geom_l_y(buffer, 228718, 216010, 225118, 1, 36, ncols, alpha);

            simdgeo::geom_l_z(buffer, 230338, 216010, 225118, 1, 36, ncols, alpha);

            simdgeo::geom_l_x(buffer, 231958, 212122, 223138, 1, 36, ncols, alpha);

            simdgeo::geom_l_y(buffer, 233578, 212122, 223138, 1, 36, ncols, alpha);

            simdgeo::geom_l_z(buffer, 235198, 212122, 223138, 1, 36, ncols, alpha);

            simdfunc::contract_primitives(buffer, 236818, 231958, 4860, ncols);

            simdfunc::contract_primitives(buffer, 241678, 227098, 4860, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 246538, 241678, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 246538, 15, nmax);

    simdtrf::transform_k_inner(buffer, 246538, 243298, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 255 * nvalues, nvalues, buffer, 246538, 15, nmax);

    simdtrf::transform_k_inner(buffer, 246538, 244918, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 510 * nvalues, nvalues, buffer, 246538, 15, nmax);

    simdtrf::transform_k_inner(buffer, 246538, 236818, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 765 * nvalues, nvalues, buffer, 246538, 15, nmax);

    simdtrf::transform_k_inner(buffer, 246538, 238438, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 1020 * nvalues, nvalues, buffer, 246538, 15, nmax);

    simdtrf::transform_k_inner(buffer, 246538, 240058, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 1275 * nvalues, nvalues, buffer, 246538, 15, nmax);
}

}  // namespace simdt2ceri
