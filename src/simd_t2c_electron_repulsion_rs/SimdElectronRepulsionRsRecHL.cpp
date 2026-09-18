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


#include "SimdElectronRepulsionRsRecHL.hpp"

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
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_hl_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 66571, 64324, 1890, nvalues);

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
                                                9, 10, 11, 12, 13}, ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 20, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 82, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 85, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 88, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 91, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 94, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 97, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 100, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 103, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 106, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 109, 0, 33, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 7, 8, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 8, 9, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 9, 10, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 10, 11, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 11, 12, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 12, 13, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 13, 14, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 14, 15, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 15, 16, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 16, 17, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 17, 18, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 178, 0, 21, 22, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 184, 0, 22, 23, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 190, 0, 23, 24, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 196, 0, 24, 25, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 202, 0, 25, 26, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 208, 0, 26, 27, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 214, 0, 27, 28, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 220, 0, 28, 29, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 226, 0, 29, 30, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 232, 0, 30, 31, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 238, 0, 31, 32, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 244, 0, 34, 37, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 254, 0, 37, 40, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 264, 0, 40, 43, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 274, 0, 43, 46, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 284, 0, 46, 49, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 294, 0, 49, 52, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 304, 0, 52, 55, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 314, 0, 55, 58, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 324, 0, 58, 61, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 334, 0, 61, 64, 166, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 344, 0, 64, 67, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 354, 0, 73, 76, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 364, 0, 76, 79, 184, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 374, 0, 79, 82, 190, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 384, 0, 82, 85, 196, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 394, 0, 85, 88, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 404, 0, 88, 91, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 414, 0, 91, 94, 214, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 424, 0, 94, 97, 220, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 434, 0, 97, 100, 226, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 444, 0, 100, 103, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 454, 0, 103, 106, 238, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 464, 0, 112, 118, 264, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 479, 0, 118, 124, 274, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 494, 0, 124, 130, 284, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 509, 0, 130, 136, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 524, 0, 136, 142, 304, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 539, 0, 142, 148, 314, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 554, 0, 148, 154, 324, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 569, 0, 154, 160, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 584, 0, 160, 166, 344, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 599, 0, 178, 184, 374, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 614, 0, 184, 190, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 629, 0, 190, 196, 394, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 644, 0, 196, 202, 404, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 659, 0, 202, 208, 414, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 674, 0, 208, 214, 424, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 689, 0, 214, 220, 434, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 704, 0, 220, 226, 444, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 719, 0, 226, 232, 454, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 734, 0, 244, 254, 464, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 755, 0, 254, 264, 479, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 776, 0, 264, 274, 494, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 797, 0, 274, 284, 509, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 818, 0, 284, 294, 524, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 839, 0, 294, 304, 539, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 860, 0, 304, 314, 554, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 881, 0, 314, 324, 569, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 902, 0, 324, 334, 584, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 923, 0, 354, 364, 599, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 944, 0, 364, 374, 614, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 965, 0, 374, 384, 629, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 986, 0, 384, 394, 644, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1007, 0, 394, 404, 659, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1028, 0, 404, 414, 674, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1049, 0, 414, 424, 689, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1070, 0, 424, 434, 704, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1091, 0, 434, 444, 719, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1112, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1115, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1118, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1121, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1124, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1127, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1130, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1133, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1136, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1139, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1142, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1145, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1148, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1151, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1154, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1157, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1160, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1163, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1166, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1169, 3, 33, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1172, 3, 9, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1181, 3, 10, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1190, 3, 11, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1199, 3, 12, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1208, 3, 13, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1217, 3, 14, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1226, 3, 15, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1235, 3, 16, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1244, 3, 17, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1253, 3, 18, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1262, 3, 23, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1271, 3, 24, 85, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1280, 3, 25, 88, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1289, 3, 26, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1298, 3, 27, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1307, 3, 28, 97, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1316, 3, 29, 100, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1325, 3, 30, 103, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1334, 3, 31, 106, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1343, 3, 32, 109, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1352, 0, 3, 40, 1172, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1370, 0, 3, 43, 1181, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1388, 0, 3, 46, 1190, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1406, 0, 3, 49, 1199, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1424, 0, 3, 52, 1208, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1442, 0, 3, 55, 1217, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1460, 0, 3, 58, 1226, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1478, 0, 3, 61, 1235, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1496, 0, 3, 64, 1244, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1514, 0, 3, 67, 1253, 172, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1532, 0, 3, 79, 1262, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1550, 0, 3, 82, 1271, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1568, 0, 3, 85, 1280, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1586, 0, 3, 88, 1289, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1604, 0, 3, 91, 1298, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1622, 0, 3, 94, 1307, 214, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1640, 0, 3, 97, 1316, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1658, 0, 3, 100, 1325, 226, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1676, 0, 3, 103, 1334, 232, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1694, 0, 3, 106, 1343, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1712, 0, 3, 118, 1370, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1742, 0, 3, 124, 1388, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1772, 0, 3, 130, 1406, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1802, 0, 3, 136, 1424, 294, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1832, 0, 3, 142, 1442, 304, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1862, 0, 3, 148, 1460, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1892, 0, 3, 154, 1478, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1922, 0, 3, 160, 1496, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1952, 0, 3, 166, 1514, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1982, 0, 3, 184, 1550, 374, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2012, 0, 3, 190, 1568, 384, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2042, 0, 3, 196, 1586, 394, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2072, 0, 3, 202, 1604, 404, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2102, 0, 3, 208, 1622, 414, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2132, 0, 3, 214, 1640, 424, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2162, 0, 3, 220, 1658, 434, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2192, 0, 3, 226, 1676, 444, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2222, 0, 3, 232, 1694, 454, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2252, 0, 3, 264, 1742, 479, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2297, 0, 3, 274, 1772, 494, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2342, 0, 3, 284, 1802, 509, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2387, 0, 3, 294, 1832, 524, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2432, 0, 3, 304, 1862, 539, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2477, 0, 3, 314, 1892, 554, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2522, 0, 3, 324, 1922, 569, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2567, 0, 3, 334, 1952, 584, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2612, 0, 3, 374, 2012, 614, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2657, 0, 3, 384, 2042, 629, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2702, 0, 3, 394, 2072, 644, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2747, 0, 3, 404, 2102, 659, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2792, 0, 3, 414, 2132, 674, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2837, 0, 3, 424, 2162, 689, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2882, 0, 3, 434, 2192, 704, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2927, 0, 3, 444, 2222, 719, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2972, 0, 3, 479, 2297, 776, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3035, 0, 3, 494, 2342, 797, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3098, 0, 3, 509, 2387, 818, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3161, 0, 3, 524, 2432, 839, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3224, 0, 3, 539, 2477, 860, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3287, 0, 3, 554, 2522, 881, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3350, 0, 3, 569, 2567, 902, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3413, 0, 3, 614, 2657, 965, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3476, 0, 3, 629, 2702, 986, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3539, 0, 3, 644, 2747, 1007, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3602, 0, 3, 659, 2792, 1028, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3665, 0, 3, 674, 2837, 1049, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3728, 0, 3, 689, 2882, 1070, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3791, 0, 3, 704, 2927, 1091, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3854, 3, 9, 10, 1115, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3860, 3, 10, 11, 1118, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3866, 3, 11, 12, 1121, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3872, 3, 12, 13, 1124, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3878, 3, 13, 14, 1127, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3884, 3, 14, 15, 1130, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3890, 3, 15, 16, 1133, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3896, 3, 16, 17, 1136, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3902, 3, 17, 18, 1139, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3908, 3, 23, 24, 1145, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3914, 3, 24, 25, 1148, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3920, 3, 25, 26, 1151, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3926, 3, 26, 27, 1154, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3932, 3, 27, 28, 1157, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3938, 3, 28, 29, 1160, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3944, 3, 29, 30, 1163, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3950, 3, 30, 31, 1166, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 3956, 3, 31, 32, 1169, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 3962, 0, 3, 1112, 3854, 1181, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3980, 0, 3, 1115, 3860, 1190, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 3998, 0, 3, 1118, 3866, 1199, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4016, 0, 3, 1121, 3872, 1208, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4034, 0, 3, 1124, 3878, 1217, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4052, 0, 3, 1127, 3884, 1226, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4070, 0, 3, 1130, 3890, 1235, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4088, 0, 3, 1133, 3896, 1244, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4106, 0, 3, 1136, 3902, 1253, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4124, 0, 3, 1142, 3908, 1271, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4142, 0, 3, 1145, 3914, 1280, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4160, 0, 3, 1148, 3920, 1289, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4178, 0, 3, 1151, 3926, 1298, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4196, 0, 3, 1154, 3932, 1307, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4214, 0, 3, 1157, 3938, 1316, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4232, 0, 3, 1160, 3944, 1325, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4250, 0, 3, 1163, 3950, 1334, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4268, 0, 3, 1166, 3956, 1343, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 4286, 0, 3, 1172, 3962, 112, 118, 1370,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4322, 0, 3, 1181, 3980, 118, 124, 1388,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4358, 0, 3, 1190, 3998, 124, 130, 1406,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4394, 0, 3, 1199, 4016, 130, 136, 1424,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4430, 0, 3, 1208, 4034, 136, 142, 1442,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4466, 0, 3, 1217, 4052, 142, 148, 1460,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4502, 0, 3, 1226, 4070, 148, 154, 1478,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4538, 0, 3, 1235, 4088, 154, 160, 1496,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4574, 0, 3, 1244, 4106, 160, 166, 1514,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4610, 0, 3, 1262, 4124, 178, 184, 1550,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4646, 0, 3, 1271, 4142, 184, 190, 1568,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4682, 0, 3, 1280, 4160, 190, 196, 1586,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4718, 0, 3, 1289, 4178, 196, 202, 1604,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4754, 0, 3, 1298, 4196, 202, 208, 1622,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4790, 0, 3, 1307, 4214, 208, 214, 1640,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4826, 0, 3, 1316, 4232, 214, 220, 1658,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4862, 0, 3, 1325, 4250, 220, 226, 1676,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 4898, 0, 3, 1334, 4268, 226, 232, 1694,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4934, 0, 3, 1352, 4286, 244, 254, 1712,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4994, 0, 3, 1370, 4322, 254, 264, 1742,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5054, 0, 3, 1388, 4358, 264, 274, 1772,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5114, 0, 3, 1406, 4394, 274, 284, 1802,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5174, 0, 3, 1424, 4430, 284, 294, 1832,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5234, 0, 3, 1442, 4466, 294, 304, 1862,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5294, 0, 3, 1460, 4502, 304, 314, 1892,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5354, 0, 3, 1478, 4538, 314, 324, 1922,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5414, 0, 3, 1496, 4574, 324, 334, 1952,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5474, 0, 3, 1532, 4610, 354, 364, 1982,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5534, 0, 3, 1550, 4646, 364, 374, 2012,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5594, 0, 3, 1568, 4682, 374, 384, 2042,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5654, 0, 3, 1586, 4718, 384, 394, 2072,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5714, 0, 3, 1604, 4754, 394, 404, 2102,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5774, 0, 3, 1622, 4790, 404, 414, 2132,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5834, 0, 3, 1640, 4826, 414, 424, 2162,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5894, 0, 3, 1658, 4862, 424, 434, 2192,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5954, 0, 3, 1676, 4898, 434, 444, 2222,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6014, 0, 3, 4286, 4322, 1742, 5054, 464,
                                                 479, 2297, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6104, 0, 3, 4322, 4358, 1772, 5114, 479,
                                                 494, 2342, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6194, 0, 3, 4358, 4394, 1802, 5174, 494,
                                                 509, 2387, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6284, 0, 3, 4394, 4430, 1832, 5234, 509,
                                                 524, 2432, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6374, 0, 3, 4430, 4466, 1862, 5294, 524,
                                                 539, 2477, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6464, 0, 3, 4466, 4502, 1892, 5354, 539,
                                                 554, 2522, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6554, 0, 3, 4502, 4538, 1922, 5414, 554,
                                                 569, 2567, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6644, 0, 3, 4610, 4646, 2012, 5594, 599,
                                                 614, 2657, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6734, 0, 3, 4646, 4682, 2042, 5654, 614,
                                                 629, 2702, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6824, 0, 3, 4682, 4718, 2072, 5714, 629,
                                                 644, 2747, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6914, 0, 3, 4718, 4754, 2102, 5774, 644,
                                                 659, 2792, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7004, 0, 3, 4754, 4790, 2132, 5834, 659,
                                                 674, 2837, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7094, 0, 3, 4790, 4826, 2162, 5894, 674,
                                                 689, 2882, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7184, 0, 3, 4826, 4862, 2192, 5954, 689,
                                                 704, 2927, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7274, 0, 3, 4934, 4994, 2252, 6014, 734,
                                                 755, 2972, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7400, 0, 3, 4994, 5054, 2297, 6104, 755,
                                                 776, 3035, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7526, 0, 3, 5054, 5114, 2342, 6194, 776,
                                                 797, 3098, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7652, 0, 3, 5114, 5174, 2387, 6284, 797,
                                                 818, 3161, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7778, 0, 3, 5174, 5234, 2432, 6374, 818,
                                                 839, 3224, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7904, 0, 3, 5234, 5294, 2477, 6464, 839,
                                                 860, 3287, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8030, 0, 3, 5294, 5354, 2522, 6554, 860,
                                                 881, 3350, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8156, 0, 3, 5474, 5534, 2612, 6644, 923,
                                                 944, 3413, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8282, 0, 3, 5534, 5594, 2657, 6734, 944,
                                                 965, 3476, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8408, 0, 3, 5594, 5654, 2702, 6824, 965,
                                                 986, 3539, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8534, 0, 3, 5654, 5714, 2747, 6914, 986,
                                                 1007, 3602, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8660, 0, 3, 5714, 5774, 2792, 7004,
                                                 1007, 1028, 3665, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8786, 0, 3, 5774, 5834, 2837, 7094,
                                                 1028, 1049, 3728, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8912, 0, 3, 5834, 5894, 2882, 7184,
                                                 1049, 1070, 3791, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9038, 3, 1112, 1115, 3860, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9048, 3, 1115, 1118, 3866, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9058, 3, 1118, 1121, 3872, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9068, 3, 1121, 1124, 3878, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9078, 3, 1124, 1127, 3884, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9088, 3, 1127, 1130, 3890, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9098, 3, 1130, 1133, 3896, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9108, 3, 1133, 1136, 3902, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9118, 3, 1142, 1145, 3914, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9128, 3, 1145, 1148, 3920, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9138, 3, 1148, 1151, 3926, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9148, 3, 1151, 1154, 3932, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9158, 3, 1154, 1157, 3938, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9168, 3, 1157, 1160, 3944, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9178, 3, 1160, 1163, 3950, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 9188, 3, 1163, 1166, 3956, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 9198, 0, 3, 3854, 9038, 3980, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9228, 0, 3, 3860, 9048, 3998, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9258, 0, 3, 3866, 9058, 4016, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9288, 0, 3, 3872, 9068, 4034, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9318, 0, 3, 3878, 9078, 4052, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9348, 0, 3, 3884, 9088, 4070, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9378, 0, 3, 3890, 9098, 4088, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9408, 0, 3, 3896, 9108, 4106, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9438, 0, 3, 3908, 9118, 4142, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9468, 0, 3, 3914, 9128, 4160, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9498, 0, 3, 3920, 9138, 4178, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9528, 0, 3, 3926, 9148, 4196, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9558, 0, 3, 3932, 9158, 4214, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9588, 0, 3, 3938, 9168, 4232, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9618, 0, 3, 3944, 9178, 4250, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 9648, 0, 3, 3950, 9188, 4268, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 9678, 0, 3, 3962, 9198, 1352, 1370,
                                                 4322, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9738, 0, 3, 3980, 9228, 1370, 1388,
                                                 4358, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9798, 0, 3, 3998, 9258, 1388, 1406,
                                                 4394, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9858, 0, 3, 4016, 9288, 1406, 1424,
                                                 4430, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9918, 0, 3, 4034, 9318, 1424, 1442,
                                                 4466, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 9978, 0, 3, 4052, 9348, 1442, 1460,
                                                 4502, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10038, 0, 3, 4070, 9378, 1460, 1478,
                                                 4538, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10098, 0, 3, 4088, 9408, 1478, 1496,
                                                 4574, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10158, 0, 3, 4124, 9438, 1532, 1550,
                                                 4646, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10218, 0, 3, 4142, 9468, 1550, 1568,
                                                 4682, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10278, 0, 3, 4160, 9498, 1568, 1586,
                                                 4718, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10338, 0, 3, 4178, 9528, 1586, 1604,
                                                 4754, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10398, 0, 3, 4196, 9558, 1604, 1622,
                                                 4790, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10458, 0, 3, 4214, 9588, 1622, 1640,
                                                 4826, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10518, 0, 3, 4232, 9618, 1640, 1658,
                                                 4862, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 10578, 0, 3, 4250, 9648, 1658, 1676,
                                                 4898, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10638, 0, 3, 4322, 9738, 1712, 1742,
                                                 5054, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10738, 0, 3, 4358, 9798, 1742, 1772,
                                                 5114, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10838, 0, 3, 4394, 9858, 1772, 1802,
                                                 5174, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 10938, 0, 3, 4430, 9918, 1802, 1832,
                                                 5234, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11038, 0, 3, 4466, 9978, 1832, 1862,
                                                 5294, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11138, 0, 3, 4502, 10038, 1862, 1892,
                                                 5354, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11238, 0, 3, 4538, 10098, 1892, 1922,
                                                 5414, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11338, 0, 3, 4646, 10218, 1982, 2012,
                                                 5594, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11438, 0, 3, 4682, 10278, 2012, 2042,
                                                 5654, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11538, 0, 3, 4718, 10338, 2042, 2072,
                                                 5714, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11638, 0, 3, 4754, 10398, 2072, 2102,
                                                 5774, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11738, 0, 3, 4790, 10458, 2102, 2132,
                                                 5834, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11838, 0, 3, 4826, 10518, 2132, 2162,
                                                 5894, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 11938, 0, 3, 4862, 10578, 2162, 2192,
                                                 5954, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12038, 0, 3, 9678, 9738, 5054, 10738,
                                                 2252, 2297, 6104, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12188, 0, 3, 9738, 9798, 5114, 10838,
                                                 2297, 2342, 6194, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12338, 0, 3, 9798, 9858, 5174, 10938,
                                                 2342, 2387, 6284, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12488, 0, 3, 9858, 9918, 5234, 11038,
                                                 2387, 2432, 6374, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12638, 0, 3, 9918, 9978, 5294, 11138,
                                                 2432, 2477, 6464, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12788, 0, 3, 9978, 10038, 5354, 11238,
                                                 2477, 2522, 6554, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 12938, 0, 3, 10158, 10218, 5594, 11438,
                                                 2612, 2657, 6734, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13088, 0, 3, 10218, 10278, 5654, 11538,
                                                 2657, 2702, 6824, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13238, 0, 3, 10278, 10338, 5714, 11638,
                                                 2702, 2747, 6914, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13388, 0, 3, 10338, 10398, 5774, 11738,
                                                 2747, 2792, 7004, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13538, 0, 3, 10398, 10458, 5834, 11838,
                                                 2792, 2837, 7094, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13688, 0, 3, 10458, 10518, 5894, 11938,
                                                 2837, 2882, 7184, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 13838, 0, 3, 10638, 10738, 6104, 12188,
                                                 2972, 3035, 7526, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14048, 0, 3, 10738, 10838, 6194, 12338,
                                                 3035, 3098, 7652, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14258, 0, 3, 10838, 10938, 6284, 12488,
                                                 3098, 3161, 7778, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14468, 0, 3, 10938, 11038, 6374, 12638,
                                                 3161, 3224, 7904, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14678, 0, 3, 11038, 11138, 6464, 12788,
                                                 3224, 3287, 8030, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14888, 0, 3, 11338, 11438, 6734, 13088,
                                                 3413, 3476, 8408, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15098, 0, 3, 11438, 11538, 6824, 13238,
                                                 3476, 3539, 8534, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15308, 0, 3, 11538, 11638, 6914, 13388,
                                                 3539, 3602, 8660, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15518, 0, 3, 11638, 11738, 7004, 13538,
                                                 3602, 3665, 8786, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15728, 0, 3, 11738, 11838, 7094, 13688,
                                                 3665, 3728, 8912, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15938, 3, 3854, 3860, 9048, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15953, 3, 3860, 3866, 9058, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15968, 3, 3866, 3872, 9068, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15983, 3, 3872, 3878, 9078, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 15998, 3, 3878, 3884, 9088, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16013, 3, 3884, 3890, 9098, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16028, 3, 3890, 3896, 9108, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16043, 3, 3908, 3914, 9128, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16058, 3, 3914, 3920, 9138, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16073, 3, 3920, 3926, 9148, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16088, 3, 3926, 3932, 9158, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16103, 3, 3932, 3938, 9168, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16118, 3, 3938, 3944, 9178, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 16133, 3, 3944, 3950, 9188, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 16148, 0, 3, 9038, 15938, 9228, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16193, 0, 3, 9048, 15953, 9258, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16238, 0, 3, 9058, 15968, 9288, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16283, 0, 3, 9068, 15983, 9318, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16328, 0, 3, 9078, 15998, 9348, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16373, 0, 3, 9088, 16013, 9378, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16418, 0, 3, 9098, 16028, 9408, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16463, 0, 3, 9118, 16043, 9468, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16508, 0, 3, 9128, 16058, 9498, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16553, 0, 3, 9138, 16073, 9528, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16598, 0, 3, 9148, 16088, 9558, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16643, 0, 3, 9158, 16103, 9588, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16688, 0, 3, 9168, 16118, 9618, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 16733, 0, 3, 9178, 16133, 9648, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 16778, 0, 3, 9198, 16148, 4286, 4322,
                                                 9738, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16868, 0, 3, 9228, 16193, 4322, 4358,
                                                 9798, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 16958, 0, 3, 9258, 16238, 4358, 4394,
                                                 9858, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17048, 0, 3, 9288, 16283, 4394, 4430,
                                                 9918, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17138, 0, 3, 9318, 16328, 4430, 4466,
                                                 9978, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17228, 0, 3, 9348, 16373, 4466, 4502,
                                                 10038, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17318, 0, 3, 9378, 16418, 4502, 4538,
                                                 10098, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17408, 0, 3, 9438, 16463, 4610, 4646,
                                                 10218, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17498, 0, 3, 9468, 16508, 4646, 4682,
                                                 10278, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17588, 0, 3, 9498, 16553, 4682, 4718,
                                                 10338, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17678, 0, 3, 9528, 16598, 4718, 4754,
                                                 10398, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17768, 0, 3, 9558, 16643, 4754, 4790,
                                                 10458, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17858, 0, 3, 9588, 16688, 4790, 4826,
                                                 10518, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 17948, 0, 3, 9618, 16733, 4826, 4862,
                                                 10578, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18038, 0, 3, 9678, 16778, 4934, 4994,
                                                 10638, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18188, 0, 3, 9738, 16868, 4994, 5054,
                                                 10738, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18338, 0, 3, 9798, 16958, 5054, 5114,
                                                 10838, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18488, 0, 3, 9858, 17048, 5114, 5174,
                                                 10938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18638, 0, 3, 9918, 17138, 5174, 5234,
                                                 11038, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18788, 0, 3, 9978, 17228, 5234, 5294,
                                                 11138, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 18938, 0, 3, 10038, 17318, 5294, 5354,
                                                 11238, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19088, 0, 3, 10158, 17408, 5474, 5534,
                                                 11338, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19238, 0, 3, 10218, 17498, 5534, 5594,
                                                 11438, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19388, 0, 3, 10278, 17588, 5594, 5654,
                                                 11538, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19538, 0, 3, 10338, 17678, 5654, 5714,
                                                 11638, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19688, 0, 3, 10398, 17768, 5714, 5774,
                                                 11738, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19838, 0, 3, 10458, 17858, 5774, 5834,
                                                 11838, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 19988, 0, 3, 10518, 17948, 5834, 5894,
                                                 11938, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20138, 0, 3, 16778, 16868, 10738, 18338,
                                                 6014, 6104, 12188, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20363, 0, 3, 16868, 16958, 10838, 18488,
                                                 6104, 6194, 12338, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20588, 0, 3, 16958, 17048, 10938, 18638,
                                                 6194, 6284, 12488, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 20813, 0, 3, 17048, 17138, 11038, 18788,
                                                 6284, 6374, 12638, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21038, 0, 3, 17138, 17228, 11138, 18938,
                                                 6374, 6464, 12788, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21263, 0, 3, 17408, 17498, 11438, 19388,
                                                 6644, 6734, 13088, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21488, 0, 3, 17498, 17588, 11538, 19538,
                                                 6734, 6824, 13238, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21713, 0, 3, 17588, 17678, 11638, 19688,
                                                 6824, 6914, 13388, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 21938, 0, 3, 17678, 17768, 11738, 19838,
                                                 6914, 7004, 13538, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22163, 0, 3, 17768, 17858, 11838, 19988,
                                                 7004, 7094, 13688, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22388, 0, 3, 18038, 18188, 12038, 20138,
                                                 7274, 7400, 13838, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 22703, 0, 3, 18188, 18338, 12188, 20363,
                                                 7400, 7526, 14048, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23018, 0, 3, 18338, 18488, 12338, 20588,
                                                 7526, 7652, 14258, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23333, 0, 3, 18488, 18638, 12488, 20813,
                                                 7652, 7778, 14468, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23648, 0, 3, 18638, 18788, 12638, 21038,
                                                 7778, 7904, 14678, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23963, 0, 3, 19088, 19238, 12938, 21263,
                                                 8156, 8282, 14888, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24278, 0, 3, 19238, 19388, 13088, 21488,
                                                 8282, 8408, 15098, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24593, 0, 3, 19388, 19538, 13238, 21713,
                                                 8408, 8534, 15308, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24908, 0, 3, 19538, 19688, 13388, 21938,
                                                 8534, 8660, 15518, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 25223, 0, 3, 19688, 19838, 13538, 22163,
                                                 8660, 8786, 15728, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25538, 3, 9038, 9048, 15953, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25559, 3, 9048, 9058, 15968, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25580, 3, 9058, 9068, 15983, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25601, 3, 9068, 9078, 15998, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25622, 3, 9078, 9088, 16013, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25643, 3, 9088, 9098, 16028, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25664, 3, 9118, 9128, 16058, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25685, 3, 9128, 9138, 16073, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25706, 3, 9138, 9148, 16088, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25727, 3, 9148, 9158, 16103, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25748, 3, 9158, 9168, 16118, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 25769, 3, 9168, 9178, 16133, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 25790, 0, 3, 15938, 25538, 16193, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25853, 0, 3, 15953, 25559, 16238, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25916, 0, 3, 15968, 25580, 16283, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 25979, 0, 3, 15983, 25601, 16328, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 26042, 0, 3, 15998, 25622, 16373, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 26105, 0, 3, 16013, 25643, 16418, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 26168, 0, 3, 16043, 25664, 16508, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 26231, 0, 3, 16058, 25685, 16553, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 26294, 0, 3, 16073, 25706, 16598, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 26357, 0, 3, 16088, 25727, 16643, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 26420, 0, 3, 16103, 25748, 16688, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 26483, 0, 3, 16118, 25769, 16733, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 26546, 0, 3, 16148, 25790, 9678, 9738,
                                                 16868, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 26672, 0, 3, 16193, 25853, 9738, 9798,
                                                 16958, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 26798, 0, 3, 16238, 25916, 9798, 9858,
                                                 17048, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 26924, 0, 3, 16283, 25979, 9858, 9918,
                                                 17138, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 27050, 0, 3, 16328, 26042, 9918, 9978,
                                                 17228, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 27176, 0, 3, 16373, 26105, 9978, 10038,
                                                 17318, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 27302, 0, 3, 16463, 26168, 10158, 10218,
                                                 17498, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 27428, 0, 3, 16508, 26231, 10218, 10278,
                                                 17588, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 27554, 0, 3, 16553, 26294, 10278, 10338,
                                                 17678, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 27680, 0, 3, 16598, 26357, 10338, 10398,
                                                 17768, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 27806, 0, 3, 16643, 26420, 10398, 10458,
                                                 17858, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 27932, 0, 3, 16688, 26483, 10458, 10518,
                                                 17948, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 28058, 0, 3, 16868, 26672, 10638, 10738,
                                                 18338, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 28268, 0, 3, 16958, 26798, 10738, 10838,
                                                 18488, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 28478, 0, 3, 17048, 26924, 10838, 10938,
                                                 18638, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 28688, 0, 3, 17138, 27050, 10938, 11038,
                                                 18788, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 28898, 0, 3, 17228, 27176, 11038, 11138,
                                                 18938, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 29108, 0, 3, 17498, 27428, 11338, 11438,
                                                 19388, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 29318, 0, 3, 17588, 27554, 11438, 11538,
                                                 19538, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 29528, 0, 3, 17678, 27680, 11538, 11638,
                                                 19688, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 29738, 0, 3, 17768, 27806, 11638, 11738,
                                                 19838, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 29948, 0, 3, 17858, 27932, 11738, 11838,
                                                 19988, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 30158, 0, 3, 26546, 26672, 18338, 28268,
                                                 12038, 12188, 20363, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 30473, 0, 3, 26672, 26798, 18488, 28478,
                                                 12188, 12338, 20588, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 30788, 0, 3, 26798, 26924, 18638, 28688,
                                                 12338, 12488, 20813, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 31103, 0, 3, 26924, 27050, 18788, 28898,
                                                 12488, 12638, 21038, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 31418, 0, 3, 27302, 27428, 19388, 29318,
                                                 12938, 13088, 21488, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 31733, 0, 3, 27428, 27554, 19538, 29528,
                                                 13088, 13238, 21713, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 32048, 0, 3, 27554, 27680, 19688, 29738,
                                                 13238, 13388, 21938, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 32363, 0, 3, 27680, 27806, 19838, 29948,
                                                 13388, 13538, 22163, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 32678, 0, 3, 28058, 28268, 20363, 30473,
                                                 13838, 14048, 23018, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 33119, 0, 3, 28268, 28478, 20588, 30788,
                                                 14048, 14258, 23333, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 33560, 0, 3, 28478, 28688, 20813, 31103,
                                                 14258, 14468, 23648, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 34001, 0, 3, 29108, 29318, 21488, 31733,
                                                 14888, 15098, 24593, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 34442, 0, 3, 29318, 29528, 21713, 32048,
                                                 15098, 15308, 24908, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 34883, 0, 3, 29528, 29738, 21938, 32363,
                                                 15308, 15518, 25223, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35324, 3, 15938, 15953, 25559, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35352, 3, 15953, 15968, 25580, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35380, 3, 15968, 15983, 25601, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35408, 3, 15983, 15998, 25622, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35436, 3, 15998, 16013, 25643, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35464, 3, 16043, 16058, 25685, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35492, 3, 16058, 16073, 25706, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35520, 3, 16073, 16088, 25727, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35548, 3, 16088, 16103, 25748, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 35576, 3, 16103, 16118, 25769, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 35604, 0, 3, 25538, 35324, 25853, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 35688, 0, 3, 25559, 35352, 25916, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 35772, 0, 3, 25580, 35380, 25979, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 35856, 0, 3, 25601, 35408, 26042, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 35940, 0, 3, 25622, 35436, 26105, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 36024, 0, 3, 25664, 35464, 26231, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 36108, 0, 3, 25685, 35492, 26294, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 36192, 0, 3, 25706, 35520, 26357, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 36276, 0, 3, 25727, 35548, 26420, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 36360, 0, 3, 25748, 35576, 26483, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 36444, 0, 3, 25790, 35604, 16778, 16868,
                                                 26672, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 36612, 0, 3, 25853, 35688, 16868, 16958,
                                                 26798, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 36780, 0, 3, 25916, 35772, 16958, 17048,
                                                 26924, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 36948, 0, 3, 25979, 35856, 17048, 17138,
                                                 27050, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 37116, 0, 3, 26042, 35940, 17138, 17228,
                                                 27176, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 37284, 0, 3, 26168, 36024, 17408, 17498,
                                                 27428, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 37452, 0, 3, 26231, 36108, 17498, 17588,
                                                 27554, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 37620, 0, 3, 26294, 36192, 17588, 17678,
                                                 27680, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 37788, 0, 3, 26357, 36276, 17678, 17768,
                                                 27806, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 37956, 0, 3, 26420, 36360, 17768, 17858,
                                                 27932, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 38124, 0, 3, 26546, 36444, 18038, 18188,
                                                 28058, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 38404, 0, 3, 26672, 36612, 18188, 18338,
                                                 28268, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 38684, 0, 3, 26798, 36780, 18338, 18488,
                                                 28478, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 38964, 0, 3, 26924, 36948, 18488, 18638,
                                                 28688, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 39244, 0, 3, 27050, 37116, 18638, 18788,
                                                 28898, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 39524, 0, 3, 27302, 37284, 19088, 19238,
                                                 29108, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 39804, 0, 3, 27428, 37452, 19238, 19388,
                                                 29318, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 40084, 0, 3, 27554, 37620, 19388, 19538,
                                                 29528, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 40364, 0, 3, 27680, 37788, 19538, 19688,
                                                 29738, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 40644, 0, 3, 27806, 37956, 19688, 19838,
                                                 29948, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 40924, 0, 3, 36444, 36612, 28268, 38684,
                                                 20138, 20363, 30473, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 41344, 0, 3, 36612, 36780, 28478, 38964,
                                                 20363, 20588, 30788, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 41764, 0, 3, 36780, 36948, 28688, 39244,
                                                 20588, 20813, 31103, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 42184, 0, 3, 37284, 37452, 29318, 40084,
                                                 21263, 21488, 31733, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 42604, 0, 3, 37452, 37620, 29528, 40364,
                                                 21488, 21713, 32048, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 43024, 0, 3, 37620, 37788, 29738, 40644,
                                                 21713, 21938, 32363, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 43444, 0, 3, 38124, 38404, 30158, 40924,
                                                 22388, 22703, 32678, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 44032, 0, 3, 38404, 38684, 30473, 41344,
                                                 22703, 23018, 33119, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 44620, 0, 3, 38684, 38964, 30788, 41764,
                                                 23018, 23333, 33560, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 45208, 0, 3, 39524, 39804, 31418, 42184,
                                                 23963, 24278, 34001, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 45796, 0, 3, 39804, 40084, 31733, 42604,
                                                 24278, 24593, 34442, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 46384, 0, 3, 40084, 40364, 32048, 43024,
                                                 24593, 24908, 34883, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 46972, 3, 25538, 25559, 35352, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 47008, 3, 25559, 25580, 35380, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 47044, 3, 25580, 25601, 35408, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 47080, 3, 25601, 25622, 35436, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 47116, 3, 25664, 25685, 35492, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 47152, 3, 25685, 25706, 35520, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 47188, 3, 25706, 25727, 35548, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 47224, 3, 25727, 25748, 35576, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 47260, 0, 3, 35324, 46972, 35688, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 47368, 0, 3, 35352, 47008, 35772, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 47476, 0, 3, 35380, 47044, 35856, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 47584, 0, 3, 35408, 47080, 35940, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 47692, 0, 3, 35464, 47116, 36108, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 47800, 0, 3, 35492, 47152, 36192, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 47908, 0, 3, 35520, 47188, 36276, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 48016, 0, 3, 35548, 47224, 36360, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 48124, 0, 3, 35604, 47260, 26546, 26672,
                                                 36612, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 48340, 0, 3, 35688, 47368, 26672, 26798,
                                                 36780, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 48556, 0, 3, 35772, 47476, 26798, 26924,
                                                 36948, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 48772, 0, 3, 35856, 47584, 26924, 27050,
                                                 37116, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 48988, 0, 3, 36024, 47692, 27302, 27428,
                                                 37452, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 49204, 0, 3, 36108, 47800, 27428, 27554,
                                                 37620, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 49420, 0, 3, 36192, 47908, 27554, 27680,
                                                 37788, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 49636, 0, 3, 36276, 48016, 27680, 27806,
                                                 37956, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 49852, 0, 3, 36612, 48340, 28058, 28268,
                                                 38684, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 50212, 0, 3, 36780, 48556, 28268, 28478,
                                                 38964, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 50572, 0, 3, 36948, 48772, 28478, 28688,
                                                 39244, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 50932, 0, 3, 37452, 49204, 29108, 29318,
                                                 40084, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 51292, 0, 3, 37620, 49420, 29318, 29528,
                                                 40364, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 51652, 0, 3, 37788, 49636, 29528, 29738,
                                                 40644, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 52012, 0, 3, 48124, 48340, 38684, 50212,
                                                 30158, 30473, 41344, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 52552, 0, 3, 48340, 48556, 38964, 50572,
                                                 30473, 30788, 41764, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 53092, 0, 3, 48988, 49204, 40084, 51292,
                                                 31418, 31733, 42604, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 53632, 0, 3, 49204, 49420, 40364, 51652,
                                                 31733, 32048, 43024, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 54172, 0, 3, 49852, 50212, 41344, 52552,
                                                 32678, 33119, 44620, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 54928, 0, 3, 50932, 51292, 42604, 53632,
                                                 34001, 34442, 46384, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 55684, 3, 35324, 35352, 47008, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 55729, 3, 35352, 35380, 47044, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 55774, 3, 35380, 35408, 47080, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 55819, 3, 35464, 35492, 47152, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 55864, 3, 35492, 35520, 47188, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 55909, 3, 35520, 35548, 47224, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 55954, 0, 3, 46972, 55684, 47368, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 56089, 0, 3, 47008, 55729, 47476, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 56224, 0, 3, 47044, 55774, 47584, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 56359, 0, 3, 47116, 55819, 47800, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 56494, 0, 3, 47152, 55864, 47908, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 56629, 0, 3, 47188, 55909, 48016, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 56764, 0, 3, 47260, 55954, 36444, 36612,
                                                 48340, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 57034, 0, 3, 47368, 56089, 36612, 36780,
                                                 48556, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 57304, 0, 3, 47476, 56224, 36780, 36948,
                                                 48772, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 57574, 0, 3, 47692, 56359, 37284, 37452,
                                                 49204, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 57844, 0, 3, 47800, 56494, 37452, 37620,
                                                 49420, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 58114, 0, 3, 47908, 56629, 37620, 37788,
                                                 49636, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 58384, 0, 3, 48124, 56764, 38124, 38404,
                                                 49852, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 58834, 0, 3, 48340, 57034, 38404, 38684,
                                                 50212, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 59284, 0, 3, 48556, 57304, 38684, 38964,
                                                 50572, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 59734, 0, 3, 48988, 57574, 39524, 39804,
                                                 50932, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 60184, 0, 3, 49204, 57844, 39804, 40084,
                                                 51292, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 60634, 0, 3, 49420, 58114, 40084, 40364,
                                                 51652, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 61084, 0, 3, 56764, 57034, 50212, 59284,
                                                 40924, 41344, 52552, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 61759, 0, 3, 57574, 57844, 51292, 60634,
                                                 42184, 42604, 53632, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 62434, 0, 3, 58384, 58834, 52012, 61084,
                                                 43444, 44032, 54172, ncols, alpha, beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 63379, 0, 3, 59734, 60184, 53092, 61759,
                                                 45208, 45796, 54928, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 64324, 62434, 1890, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 66214, 65269, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 66214, 17, nmax);

    simdtrf::transform_l_inner(buffer, 66214, 64324, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 187 * nvalues, nvalues, buffer, 66214, 17, nmax);
}

}  // namespace simdt2ceri
