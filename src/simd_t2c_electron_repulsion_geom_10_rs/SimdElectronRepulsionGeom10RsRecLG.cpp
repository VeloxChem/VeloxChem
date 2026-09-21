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


#include "SimdElectronRepulsionGeom10RsRecLG.hpp"

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
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryL1.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_lg_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_lg_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 63597, 59142, 4050, nvalues);

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

            compute_prim_is_electron_repulsion_0(buffer, 1112, 0, 464, 479, 776, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1140, 0, 479, 494, 797, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1168, 0, 494, 509, 818, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1196, 0, 509, 524, 839, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1224, 0, 524, 539, 860, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1252, 0, 539, 554, 881, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1280, 0, 554, 569, 902, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1308, 0, 599, 614, 965, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1336, 0, 614, 629, 986, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1364, 0, 629, 644, 1007, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1392, 0, 644, 659, 1028, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1420, 0, 659, 674, 1049, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1448, 0, 674, 689, 1070, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1476, 0, 689, 704, 1091, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1504, 0, 734, 755, 1112, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1540, 0, 755, 776, 1140, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1576, 0, 776, 797, 1168, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1612, 0, 797, 818, 1196, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1648, 0, 818, 839, 1224, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1684, 0, 839, 860, 1252, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1720, 0, 860, 881, 1280, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1756, 0, 923, 944, 1308, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1792, 0, 944, 965, 1336, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1828, 0, 965, 986, 1364, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1864, 0, 986, 1007, 1392, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1900, 0, 1007, 1028, 1420, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1936, 0, 1028, 1049, 1448, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1972, 0, 1049, 1070, 1476, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2008, 0, 1112, 1140, 1576, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2053, 0, 1140, 1168, 1612, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2098, 0, 1168, 1196, 1648, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2143, 0, 1196, 1224, 1684, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2188, 0, 1224, 1252, 1720, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2233, 0, 1308, 1336, 1828, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2278, 0, 1336, 1364, 1864, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2323, 0, 1364, 1392, 1900, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2368, 0, 1392, 1420, 1936, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2413, 0, 1420, 1448, 1972, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2458, 0, 1504, 1540, 2008, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2513, 0, 1540, 1576, 2053, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2568, 0, 1576, 1612, 2098, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2623, 0, 1612, 1648, 2143, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2678, 0, 1648, 1684, 2188, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2733, 0, 1756, 1792, 2233, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2788, 0, 1792, 1828, 2278, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2843, 0, 1828, 1864, 2323, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2898, 0, 1864, 1900, 2368, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 2953, 0, 1900, 1936, 2413, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 3008, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3011, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3014, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3017, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3020, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3023, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3026, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3029, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3032, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3035, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3038, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3041, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3044, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3047, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3050, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3053, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3056, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3059, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3062, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3065, 3, 33, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 3068, 3, 9, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3077, 3, 10, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3086, 3, 11, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3095, 3, 12, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3104, 3, 13, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3113, 3, 14, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3122, 3, 15, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3131, 3, 16, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3140, 3, 17, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3149, 3, 18, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3158, 3, 23, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3167, 3, 24, 85, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3176, 3, 25, 88, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3185, 3, 26, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3194, 3, 27, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3203, 3, 28, 97, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3212, 3, 29, 100, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3221, 3, 30, 103, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3230, 3, 31, 106, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3239, 3, 32, 109, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3248, 0, 3, 40, 3068, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3266, 0, 3, 43, 3077, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3284, 0, 3, 46, 3086, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3302, 0, 3, 49, 3095, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3320, 0, 3, 52, 3104, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3338, 0, 3, 55, 3113, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3356, 0, 3, 58, 3122, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3374, 0, 3, 61, 3131, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3392, 0, 3, 64, 3140, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3410, 0, 3, 67, 3149, 172, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3428, 0, 3, 79, 3158, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3446, 0, 3, 82, 3167, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3464, 0, 3, 85, 3176, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3482, 0, 3, 88, 3185, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3500, 0, 3, 91, 3194, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3518, 0, 3, 94, 3203, 214, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3536, 0, 3, 97, 3212, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3554, 0, 3, 100, 3221, 226, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3572, 0, 3, 103, 3230, 232, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3590, 0, 3, 106, 3239, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3608, 0, 3, 118, 3266, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3638, 0, 3, 124, 3284, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3668, 0, 3, 130, 3302, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3698, 0, 3, 136, 3320, 294, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3728, 0, 3, 142, 3338, 304, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3758, 0, 3, 148, 3356, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3788, 0, 3, 154, 3374, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3818, 0, 3, 160, 3392, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3848, 0, 3, 166, 3410, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3878, 0, 3, 184, 3446, 374, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3908, 0, 3, 190, 3464, 384, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3938, 0, 3, 196, 3482, 394, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3968, 0, 3, 202, 3500, 404, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3998, 0, 3, 208, 3518, 414, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4028, 0, 3, 214, 3536, 424, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4058, 0, 3, 220, 3554, 434, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4088, 0, 3, 226, 3572, 444, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4118, 0, 3, 232, 3590, 454, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4148, 0, 3, 264, 3638, 479, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4193, 0, 3, 274, 3668, 494, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4238, 0, 3, 284, 3698, 509, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4283, 0, 3, 294, 3728, 524, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4328, 0, 3, 304, 3758, 539, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4373, 0, 3, 314, 3788, 554, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4418, 0, 3, 324, 3818, 569, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4463, 0, 3, 334, 3848, 584, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4508, 0, 3, 374, 3908, 614, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4553, 0, 3, 384, 3938, 629, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4598, 0, 3, 394, 3968, 644, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4643, 0, 3, 404, 3998, 659, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4688, 0, 3, 414, 4028, 674, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4733, 0, 3, 424, 4058, 689, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4778, 0, 3, 434, 4088, 704, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4823, 0, 3, 444, 4118, 719, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4868, 0, 3, 479, 4193, 776, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4931, 0, 3, 494, 4238, 797, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4994, 0, 3, 509, 4283, 818, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5057, 0, 3, 524, 4328, 839, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5120, 0, 3, 539, 4373, 860, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5183, 0, 3, 554, 4418, 881, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5246, 0, 3, 569, 4463, 902, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5309, 0, 3, 614, 4553, 965, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5372, 0, 3, 629, 4598, 986, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5435, 0, 3, 644, 4643, 1007, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5498, 0, 3, 659, 4688, 1028, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5561, 0, 3, 674, 4733, 1049, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5624, 0, 3, 689, 4778, 1070, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5687, 0, 3, 704, 4823, 1091, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5750, 0, 3, 776, 4931, 1140, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5834, 0, 3, 797, 4994, 1168, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5918, 0, 3, 818, 5057, 1196, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6002, 0, 3, 839, 5120, 1224, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6086, 0, 3, 860, 5183, 1252, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6170, 0, 3, 881, 5246, 1280, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6254, 0, 3, 965, 5372, 1336, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6338, 0, 3, 986, 5435, 1364, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6422, 0, 3, 1007, 5498, 1392, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6506, 0, 3, 1028, 5561, 1420, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6590, 0, 3, 1049, 5624, 1448, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6674, 0, 3, 1070, 5687, 1476, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6758, 0, 3, 1140, 5834, 1576, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6866, 0, 3, 1168, 5918, 1612, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6974, 0, 3, 1196, 6002, 1648, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7082, 0, 3, 1224, 6086, 1684, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7190, 0, 3, 1252, 6170, 1720, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7298, 0, 3, 1336, 6338, 1828, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7406, 0, 3, 1364, 6422, 1864, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7514, 0, 3, 1392, 6506, 1900, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7622, 0, 3, 1420, 6590, 1936, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7730, 0, 3, 1448, 6674, 1972, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7838, 0, 3, 1576, 6866, 2053, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7973, 0, 3, 1612, 6974, 2098, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8108, 0, 3, 1648, 7082, 2143, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8243, 0, 3, 1684, 7190, 2188, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8378, 0, 3, 1828, 7406, 2278, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8513, 0, 3, 1864, 7514, 2323, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8648, 0, 3, 1900, 7622, 2368, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8783, 0, 3, 1936, 7730, 2413, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 8918, 0, 3, 2053, 7973, 2568, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 9083, 0, 3, 2098, 8108, 2623, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 9248, 0, 3, 2143, 8243, 2678, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 9413, 0, 3, 2278, 8513, 2843, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 9578, 0, 3, 2323, 8648, 2898, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 9743, 0, 3, 2368, 8783, 2953, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 9908, 3, 9, 10, 3011, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9914, 3, 10, 11, 3014, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9920, 3, 11, 12, 3017, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9926, 3, 12, 13, 3020, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9932, 3, 13, 14, 3023, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9938, 3, 14, 15, 3026, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9944, 3, 15, 16, 3029, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9950, 3, 16, 17, 3032, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9956, 3, 17, 18, 3035, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9962, 3, 23, 24, 3041, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9968, 3, 24, 25, 3044, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9974, 3, 25, 26, 3047, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9980, 3, 26, 27, 3050, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9986, 3, 27, 28, 3053, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9992, 3, 28, 29, 3056, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9998, 3, 29, 30, 3059, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 10004, 3, 30, 31, 3062, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 10010, 3, 31, 32, 3065, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 10016, 0, 3, 3008, 9908, 3077, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10034, 0, 3, 3011, 9914, 3086, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10052, 0, 3, 3014, 9920, 3095, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10070, 0, 3, 3017, 9926, 3104, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10088, 0, 3, 3020, 9932, 3113, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10106, 0, 3, 3023, 9938, 3122, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10124, 0, 3, 3026, 9944, 3131, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10142, 0, 3, 3029, 9950, 3140, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10160, 0, 3, 3032, 9956, 3149, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10178, 0, 3, 3038, 9962, 3167, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10196, 0, 3, 3041, 9968, 3176, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10214, 0, 3, 3044, 9974, 3185, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10232, 0, 3, 3047, 9980, 3194, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10250, 0, 3, 3050, 9986, 3203, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10268, 0, 3, 3053, 9992, 3212, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10286, 0, 3, 3056, 9998, 3221, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10304, 0, 3, 3059, 10004, 3230, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10322, 0, 3, 3062, 10010, 3239, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 10340, 0, 3, 3068, 10016, 112, 118,
                                                 3266, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10376, 0, 3, 3077, 10034, 118, 124,
                                                 3284, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10412, 0, 3, 3086, 10052, 124, 130,
                                                 3302, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10448, 0, 3, 3095, 10070, 130, 136,
                                                 3320, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10484, 0, 3, 3104, 10088, 136, 142,
                                                 3338, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10520, 0, 3, 3113, 10106, 142, 148,
                                                 3356, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10556, 0, 3, 3122, 10124, 148, 154,
                                                 3374, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10592, 0, 3, 3131, 10142, 154, 160,
                                                 3392, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10628, 0, 3, 3140, 10160, 160, 166,
                                                 3410, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10664, 0, 3, 3158, 10178, 178, 184,
                                                 3446, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10700, 0, 3, 3167, 10196, 184, 190,
                                                 3464, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10736, 0, 3, 3176, 10214, 190, 196,
                                                 3482, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10772, 0, 3, 3185, 10232, 196, 202,
                                                 3500, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10808, 0, 3, 3194, 10250, 202, 208,
                                                 3518, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10844, 0, 3, 3203, 10268, 208, 214,
                                                 3536, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10880, 0, 3, 3212, 10286, 214, 220,
                                                 3554, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10916, 0, 3, 3221, 10304, 220, 226,
                                                 3572, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10952, 0, 3, 3230, 10322, 226, 232,
                                                 3590, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10988, 0, 3, 3248, 10340, 244, 254,
                                                 3608, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11048, 0, 3, 3266, 10376, 254, 264,
                                                 3638, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11108, 0, 3, 3284, 10412, 264, 274,
                                                 3668, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11168, 0, 3, 3302, 10448, 274, 284,
                                                 3698, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11228, 0, 3, 3320, 10484, 284, 294,
                                                 3728, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11288, 0, 3, 3338, 10520, 294, 304,
                                                 3758, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11348, 0, 3, 3356, 10556, 304, 314,
                                                 3788, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11408, 0, 3, 3374, 10592, 314, 324,
                                                 3818, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11468, 0, 3, 3392, 10628, 324, 334,
                                                 3848, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11528, 0, 3, 3428, 10664, 354, 364,
                                                 3878, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11588, 0, 3, 3446, 10700, 364, 374,
                                                 3908, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11648, 0, 3, 3464, 10736, 374, 384,
                                                 3938, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11708, 0, 3, 3482, 10772, 384, 394,
                                                 3968, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11768, 0, 3, 3500, 10808, 394, 404,
                                                 3998, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11828, 0, 3, 3518, 10844, 404, 414,
                                                 4028, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11888, 0, 3, 3536, 10880, 414, 424,
                                                 4058, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11948, 0, 3, 3554, 10916, 424, 434,
                                                 4088, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 12008, 0, 3, 3572, 10952, 434, 444,
                                                 4118, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12068, 0, 3, 10340, 10376, 3638, 11108,
                                                 464, 479, 4193, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12158, 0, 3, 10376, 10412, 3668, 11168,
                                                 479, 494, 4238, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12248, 0, 3, 10412, 10448, 3698, 11228,
                                                 494, 509, 4283, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12338, 0, 3, 10448, 10484, 3728, 11288,
                                                 509, 524, 4328, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12428, 0, 3, 10484, 10520, 3758, 11348,
                                                 524, 539, 4373, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12518, 0, 3, 10520, 10556, 3788, 11408,
                                                 539, 554, 4418, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12608, 0, 3, 10556, 10592, 3818, 11468,
                                                 554, 569, 4463, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12698, 0, 3, 10664, 10700, 3908, 11648,
                                                 599, 614, 4553, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12788, 0, 3, 10700, 10736, 3938, 11708,
                                                 614, 629, 4598, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12878, 0, 3, 10736, 10772, 3968, 11768,
                                                 629, 644, 4643, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12968, 0, 3, 10772, 10808, 3998, 11828,
                                                 644, 659, 4688, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13058, 0, 3, 10808, 10844, 4028, 11888,
                                                 659, 674, 4733, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13148, 0, 3, 10844, 10880, 4058, 11948,
                                                 674, 689, 4778, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13238, 0, 3, 10880, 10916, 4088, 12008,
                                                 689, 704, 4823, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13328, 0, 3, 10988, 11048, 4148, 12068,
                                                 734, 755, 4868, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13454, 0, 3, 11048, 11108, 4193, 12158,
                                                 755, 776, 4931, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13580, 0, 3, 11108, 11168, 4238, 12248,
                                                 776, 797, 4994, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13706, 0, 3, 11168, 11228, 4283, 12338,
                                                 797, 818, 5057, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13832, 0, 3, 11228, 11288, 4328, 12428,
                                                 818, 839, 5120, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13958, 0, 3, 11288, 11348, 4373, 12518,
                                                 839, 860, 5183, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14084, 0, 3, 11348, 11408, 4418, 12608,
                                                 860, 881, 5246, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14210, 0, 3, 11528, 11588, 4508, 12698,
                                                 923, 944, 5309, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14336, 0, 3, 11588, 11648, 4553, 12788,
                                                 944, 965, 5372, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14462, 0, 3, 11648, 11708, 4598, 12878,
                                                 965, 986, 5435, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14588, 0, 3, 11708, 11768, 4643, 12968,
                                                 986, 1007, 5498, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14714, 0, 3, 11768, 11828, 4688, 13058,
                                                 1007, 1028, 5561, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14840, 0, 3, 11828, 11888, 4733, 13148,
                                                 1028, 1049, 5624, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14966, 0, 3, 11888, 11948, 4778, 13238,
                                                 1049, 1070, 5687, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15092, 0, 3, 12068, 12158, 4931, 13580,
                                                 1112, 1140, 5834, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15260, 0, 3, 12158, 12248, 4994, 13706,
                                                 1140, 1168, 5918, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15428, 0, 3, 12248, 12338, 5057, 13832,
                                                 1168, 1196, 6002, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15596, 0, 3, 12338, 12428, 5120, 13958,
                                                 1196, 1224, 6086, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15764, 0, 3, 12428, 12518, 5183, 14084,
                                                 1224, 1252, 6170, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15932, 0, 3, 12698, 12788, 5372, 14462,
                                                 1308, 1336, 6338, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16100, 0, 3, 12788, 12878, 5435, 14588,
                                                 1336, 1364, 6422, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16268, 0, 3, 12878, 12968, 5498, 14714,
                                                 1364, 1392, 6506, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16436, 0, 3, 12968, 13058, 5561, 14840,
                                                 1392, 1420, 6590, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16604, 0, 3, 13058, 13148, 5624, 14966,
                                                 1420, 1448, 6674, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16772, 0, 3, 13328, 13454, 5750, 15092,
                                                 1504, 1540, 6758, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16988, 0, 3, 13454, 13580, 5834, 15260,
                                                 1540, 1576, 6866, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17204, 0, 3, 13580, 13706, 5918, 15428,
                                                 1576, 1612, 6974, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17420, 0, 3, 13706, 13832, 6002, 15596,
                                                 1612, 1648, 7082, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17636, 0, 3, 13832, 13958, 6086, 15764,
                                                 1648, 1684, 7190, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17852, 0, 3, 14210, 14336, 6254, 15932,
                                                 1756, 1792, 7298, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18068, 0, 3, 14336, 14462, 6338, 16100,
                                                 1792, 1828, 7406, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18284, 0, 3, 14462, 14588, 6422, 16268,
                                                 1828, 1864, 7514, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18500, 0, 3, 14588, 14714, 6506, 16436,
                                                 1864, 1900, 7622, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18716, 0, 3, 14714, 14840, 6590, 16604,
                                                 1900, 1936, 7730, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 18932, 0, 3, 15092, 15260, 6866, 17204,
                                                 2008, 2053, 7973, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 19202, 0, 3, 15260, 15428, 6974, 17420,
                                                 2053, 2098, 8108, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 19472, 0, 3, 15428, 15596, 7082, 17636,
                                                 2098, 2143, 8243, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 19742, 0, 3, 15932, 16100, 7406, 18284,
                                                 2233, 2278, 8513, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 20012, 0, 3, 16100, 16268, 7514, 18500,
                                                 2278, 2323, 8648, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 20282, 0, 3, 16268, 16436, 7622, 18716,
                                                 2323, 2368, 8783, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 20552, 0, 3, 16772, 16988, 7838, 18932,
                                                 2458, 2513, 8918, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 20882, 0, 3, 16988, 17204, 7973, 19202,
                                                 2513, 2568, 9083, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 21212, 0, 3, 17204, 17420, 8108, 19472,
                                                 2568, 2623, 9248, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 21542, 0, 3, 17852, 18068, 8378, 19742,
                                                 2733, 2788, 9413, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 21872, 0, 3, 18068, 18284, 8513, 20012,
                                                 2788, 2843, 9578, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 22202, 0, 3, 18284, 18500, 8648, 20282,
                                                 2843, 2898, 9743, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22532, 3, 3008, 3011, 9914, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22542, 3, 3011, 3014, 9920, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22552, 3, 3014, 3017, 9926, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22562, 3, 3017, 3020, 9932, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22572, 3, 3020, 3023, 9938, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22582, 3, 3023, 3026, 9944, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22592, 3, 3026, 3029, 9950, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22602, 3, 3029, 3032, 9956, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22612, 3, 3038, 3041, 9968, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22622, 3, 3041, 3044, 9974, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22632, 3, 3044, 3047, 9980, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22642, 3, 3047, 3050, 9986, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22652, 3, 3050, 3053, 9992, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22662, 3, 3053, 3056, 9998, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22672, 3, 3056, 3059, 10004, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22682, 3, 3059, 3062, 10010, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 22692, 0, 3, 9908, 22532, 10034, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22722, 0, 3, 9914, 22542, 10052, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22752, 0, 3, 9920, 22552, 10070, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22782, 0, 3, 9926, 22562, 10088, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22812, 0, 3, 9932, 22572, 10106, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22842, 0, 3, 9938, 22582, 10124, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22872, 0, 3, 9944, 22592, 10142, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22902, 0, 3, 9950, 22602, 10160, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22932, 0, 3, 9962, 22612, 10196, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22962, 0, 3, 9968, 22622, 10214, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22992, 0, 3, 9974, 22632, 10232, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23022, 0, 3, 9980, 22642, 10250, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23052, 0, 3, 9986, 22652, 10268, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23082, 0, 3, 9992, 22662, 10286, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23112, 0, 3, 9998, 22672, 10304, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23142, 0, 3, 10004, 22682, 10322, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 23172, 0, 3, 10016, 22692, 3248, 3266,
                                                 10376, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23232, 0, 3, 10034, 22722, 3266, 3284,
                                                 10412, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23292, 0, 3, 10052, 22752, 3284, 3302,
                                                 10448, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23352, 0, 3, 10070, 22782, 3302, 3320,
                                                 10484, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23412, 0, 3, 10088, 22812, 3320, 3338,
                                                 10520, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23472, 0, 3, 10106, 22842, 3338, 3356,
                                                 10556, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23532, 0, 3, 10124, 22872, 3356, 3374,
                                                 10592, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23592, 0, 3, 10142, 22902, 3374, 3392,
                                                 10628, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23652, 0, 3, 10178, 22932, 3428, 3446,
                                                 10700, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23712, 0, 3, 10196, 22962, 3446, 3464,
                                                 10736, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23772, 0, 3, 10214, 22992, 3464, 3482,
                                                 10772, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23832, 0, 3, 10232, 23022, 3482, 3500,
                                                 10808, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23892, 0, 3, 10250, 23052, 3500, 3518,
                                                 10844, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23952, 0, 3, 10268, 23082, 3518, 3536,
                                                 10880, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24012, 0, 3, 10286, 23112, 3536, 3554,
                                                 10916, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24072, 0, 3, 10304, 23142, 3554, 3572,
                                                 10952, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24132, 0, 3, 10376, 23232, 3608, 3638,
                                                 11108, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24232, 0, 3, 10412, 23292, 3638, 3668,
                                                 11168, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24332, 0, 3, 10448, 23352, 3668, 3698,
                                                 11228, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24432, 0, 3, 10484, 23412, 3698, 3728,
                                                 11288, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24532, 0, 3, 10520, 23472, 3728, 3758,
                                                 11348, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24632, 0, 3, 10556, 23532, 3758, 3788,
                                                 11408, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24732, 0, 3, 10592, 23592, 3788, 3818,
                                                 11468, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24832, 0, 3, 10700, 23712, 3878, 3908,
                                                 11648, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24932, 0, 3, 10736, 23772, 3908, 3938,
                                                 11708, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25032, 0, 3, 10772, 23832, 3938, 3968,
                                                 11768, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25132, 0, 3, 10808, 23892, 3968, 3998,
                                                 11828, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25232, 0, 3, 10844, 23952, 3998, 4028,
                                                 11888, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25332, 0, 3, 10880, 24012, 4028, 4058,
                                                 11948, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25432, 0, 3, 10916, 24072, 4058, 4088,
                                                 12008, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25532, 0, 3, 23172, 23232, 11108, 24232,
                                                 4148, 4193, 12158, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25682, 0, 3, 23232, 23292, 11168, 24332,
                                                 4193, 4238, 12248, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25832, 0, 3, 23292, 23352, 11228, 24432,
                                                 4238, 4283, 12338, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25982, 0, 3, 23352, 23412, 11288, 24532,
                                                 4283, 4328, 12428, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26132, 0, 3, 23412, 23472, 11348, 24632,
                                                 4328, 4373, 12518, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26282, 0, 3, 23472, 23532, 11408, 24732,
                                                 4373, 4418, 12608, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26432, 0, 3, 23652, 23712, 11648, 24932,
                                                 4508, 4553, 12788, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26582, 0, 3, 23712, 23772, 11708, 25032,
                                                 4553, 4598, 12878, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26732, 0, 3, 23772, 23832, 11768, 25132,
                                                 4598, 4643, 12968, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26882, 0, 3, 23832, 23892, 11828, 25232,
                                                 4643, 4688, 13058, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 27032, 0, 3, 23892, 23952, 11888, 25332,
                                                 4688, 4733, 13148, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 27182, 0, 3, 23952, 24012, 11948, 25432,
                                                 4733, 4778, 13238, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27332, 0, 3, 24132, 24232, 12158, 25682,
                                                 4868, 4931, 13580, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27542, 0, 3, 24232, 24332, 12248, 25832,
                                                 4931, 4994, 13706, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27752, 0, 3, 24332, 24432, 12338, 25982,
                                                 4994, 5057, 13832, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27962, 0, 3, 24432, 24532, 12428, 26132,
                                                 5057, 5120, 13958, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28172, 0, 3, 24532, 24632, 12518, 26282,
                                                 5120, 5183, 14084, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28382, 0, 3, 24832, 24932, 12788, 26582,
                                                 5309, 5372, 14462, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28592, 0, 3, 24932, 25032, 12878, 26732,
                                                 5372, 5435, 14588, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28802, 0, 3, 25032, 25132, 12968, 26882,
                                                 5435, 5498, 14714, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29012, 0, 3, 25132, 25232, 13058, 27032,
                                                 5498, 5561, 14840, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29222, 0, 3, 25232, 25332, 13148, 27182,
                                                 5561, 5624, 14966, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 29432, 0, 3, 25532, 25682, 13580, 27542,
                                                 5750, 5834, 15260, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 29712, 0, 3, 25682, 25832, 13706, 27752,
                                                 5834, 5918, 15428, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 29992, 0, 3, 25832, 25982, 13832, 27962,
                                                 5918, 6002, 15596, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 30272, 0, 3, 25982, 26132, 13958, 28172,
                                                 6002, 6086, 15764, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 30552, 0, 3, 26432, 26582, 14462, 28592,
                                                 6254, 6338, 16100, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 30832, 0, 3, 26582, 26732, 14588, 28802,
                                                 6338, 6422, 16268, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31112, 0, 3, 26732, 26882, 14714, 29012,
                                                 6422, 6506, 16436, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31392, 0, 3, 26882, 27032, 14840, 29222,
                                                 6506, 6590, 16604, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 31672, 0, 3, 27332, 27542, 15260, 29712,
                                                 6758, 6866, 17204, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 32032, 0, 3, 27542, 27752, 15428, 29992,
                                                 6866, 6974, 17420, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 32392, 0, 3, 27752, 27962, 15596, 30272,
                                                 6974, 7082, 17636, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 32752, 0, 3, 28382, 28592, 16100, 30832,
                                                 7298, 7406, 18284, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 33112, 0, 3, 28592, 28802, 16268, 31112,
                                                 7406, 7514, 18500, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 33472, 0, 3, 28802, 29012, 16436, 31392,
                                                 7514, 7622, 18716, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 33832, 0, 3, 29432, 29712, 17204, 32032,
                                                 7838, 7973, 19202, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 34282, 0, 3, 29712, 29992, 17420, 32392,
                                                 7973, 8108, 19472, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 34732, 0, 3, 30552, 30832, 18284, 33112,
                                                 8378, 8513, 20012, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 35182, 0, 3, 30832, 31112, 18500, 33472,
                                                 8513, 8648, 20282, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 35632, 0, 3, 31672, 32032, 19202, 34282,
                                                 8918, 9083, 21212, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 36182, 0, 3, 32752, 33112, 20012, 35182,
                                                 9413, 9578, 22202, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36732, 3, 9908, 9914, 22542, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36747, 3, 9914, 9920, 22552, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36762, 3, 9920, 9926, 22562, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36777, 3, 9926, 9932, 22572, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36792, 3, 9932, 9938, 22582, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36807, 3, 9938, 9944, 22592, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36822, 3, 9944, 9950, 22602, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36837, 3, 9962, 9968, 22622, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36852, 3, 9968, 9974, 22632, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36867, 3, 9974, 9980, 22642, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36882, 3, 9980, 9986, 22652, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36897, 3, 9986, 9992, 22662, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36912, 3, 9992, 9998, 22672, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36927, 3, 9998, 10004, 22682, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 36942, 0, 3, 22532, 36732, 22722, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 36987, 0, 3, 22542, 36747, 22752, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37032, 0, 3, 22552, 36762, 22782, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37077, 0, 3, 22562, 36777, 22812, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37122, 0, 3, 22572, 36792, 22842, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37167, 0, 3, 22582, 36807, 22872, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37212, 0, 3, 22592, 36822, 22902, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37257, 0, 3, 22612, 36837, 22962, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37302, 0, 3, 22622, 36852, 22992, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37347, 0, 3, 22632, 36867, 23022, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37392, 0, 3, 22642, 36882, 23052, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37437, 0, 3, 22652, 36897, 23082, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37482, 0, 3, 22662, 36912, 23112, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37527, 0, 3, 22672, 36927, 23142, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 37572, 0, 3, 22692, 36942, 10340, 10376,
                                                 23232, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37662, 0, 3, 22722, 36987, 10376, 10412,
                                                 23292, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37752, 0, 3, 22752, 37032, 10412, 10448,
                                                 23352, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37842, 0, 3, 22782, 37077, 10448, 10484,
                                                 23412, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37932, 0, 3, 22812, 37122, 10484, 10520,
                                                 23472, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38022, 0, 3, 22842, 37167, 10520, 10556,
                                                 23532, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38112, 0, 3, 22872, 37212, 10556, 10592,
                                                 23592, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38202, 0, 3, 22932, 37257, 10664, 10700,
                                                 23712, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38292, 0, 3, 22962, 37302, 10700, 10736,
                                                 23772, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38382, 0, 3, 22992, 37347, 10736, 10772,
                                                 23832, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38472, 0, 3, 23022, 37392, 10772, 10808,
                                                 23892, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38562, 0, 3, 23052, 37437, 10808, 10844,
                                                 23952, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38652, 0, 3, 23082, 37482, 10844, 10880,
                                                 24012, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38742, 0, 3, 23112, 37527, 10880, 10916,
                                                 24072, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38832, 0, 3, 23172, 37572, 10988, 11048,
                                                 24132, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 38982, 0, 3, 23232, 37662, 11048, 11108,
                                                 24232, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39132, 0, 3, 23292, 37752, 11108, 11168,
                                                 24332, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39282, 0, 3, 23352, 37842, 11168, 11228,
                                                 24432, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39432, 0, 3, 23412, 37932, 11228, 11288,
                                                 24532, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39582, 0, 3, 23472, 38022, 11288, 11348,
                                                 24632, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39732, 0, 3, 23532, 38112, 11348, 11408,
                                                 24732, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39882, 0, 3, 23652, 38202, 11528, 11588,
                                                 24832, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40032, 0, 3, 23712, 38292, 11588, 11648,
                                                 24932, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40182, 0, 3, 23772, 38382, 11648, 11708,
                                                 25032, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40332, 0, 3, 23832, 38472, 11708, 11768,
                                                 25132, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40482, 0, 3, 23892, 38562, 11768, 11828,
                                                 25232, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40632, 0, 3, 23952, 38652, 11828, 11888,
                                                 25332, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40782, 0, 3, 24012, 38742, 11888, 11948,
                                                 25432, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 40932, 0, 3, 37572, 37662, 24232, 39132,
                                                 12068, 12158, 25682, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 41157, 0, 3, 37662, 37752, 24332, 39282,
                                                 12158, 12248, 25832, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 41382, 0, 3, 37752, 37842, 24432, 39432,
                                                 12248, 12338, 25982, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 41607, 0, 3, 37842, 37932, 24532, 39582,
                                                 12338, 12428, 26132, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 41832, 0, 3, 37932, 38022, 24632, 39732,
                                                 12428, 12518, 26282, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42057, 0, 3, 38202, 38292, 24932, 40182,
                                                 12698, 12788, 26582, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42282, 0, 3, 38292, 38382, 25032, 40332,
                                                 12788, 12878, 26732, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42507, 0, 3, 38382, 38472, 25132, 40482,
                                                 12878, 12968, 26882, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42732, 0, 3, 38472, 38562, 25232, 40632,
                                                 12968, 13058, 27032, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42957, 0, 3, 38562, 38652, 25332, 40782,
                                                 13058, 13148, 27182, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 43182, 0, 3, 38832, 38982, 25532, 40932,
                                                 13328, 13454, 27332, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 43497, 0, 3, 38982, 39132, 25682, 41157,
                                                 13454, 13580, 27542, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 43812, 0, 3, 39132, 39282, 25832, 41382,
                                                 13580, 13706, 27752, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 44127, 0, 3, 39282, 39432, 25982, 41607,
                                                 13706, 13832, 27962, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 44442, 0, 3, 39432, 39582, 26132, 41832,
                                                 13832, 13958, 28172, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 44757, 0, 3, 39882, 40032, 26432, 42057,
                                                 14210, 14336, 28382, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 45072, 0, 3, 40032, 40182, 26582, 42282,
                                                 14336, 14462, 28592, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 45387, 0, 3, 40182, 40332, 26732, 42507,
                                                 14462, 14588, 28802, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 45702, 0, 3, 40332, 40482, 26882, 42732,
                                                 14588, 14714, 29012, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 46017, 0, 3, 40482, 40632, 27032, 42957,
                                                 14714, 14840, 29222, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 46332, 0, 3, 40932, 41157, 27542, 43812,
                                                 15092, 15260, 29712, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 46752, 0, 3, 41157, 41382, 27752, 44127,
                                                 15260, 15428, 29992, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 47172, 0, 3, 41382, 41607, 27962, 44442,
                                                 15428, 15596, 30272, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 47592, 0, 3, 42057, 42282, 28592, 45387,
                                                 15932, 16100, 30832, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 48012, 0, 3, 42282, 42507, 28802, 45702,
                                                 16100, 16268, 31112, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 48432, 0, 3, 42507, 42732, 29012, 46017,
                                                 16268, 16436, 31392, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 48852, 0, 3, 43182, 43497, 29432, 46332,
                                                 16772, 16988, 31672, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 49392, 0, 3, 43497, 43812, 29712, 46752,
                                                 16988, 17204, 32032, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 49932, 0, 3, 43812, 44127, 29992, 47172,
                                                 17204, 17420, 32392, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 50472, 0, 3, 44757, 45072, 30552, 47592,
                                                 17852, 18068, 32752, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 51012, 0, 3, 45072, 45387, 30832, 48012,
                                                 18068, 18284, 33112, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 51552, 0, 3, 45387, 45702, 31112, 48432,
                                                 18284, 18500, 33472, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 52092, 0, 3, 46332, 46752, 32032, 49932,
                                                 18932, 19202, 34282, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 52767, 0, 3, 47592, 48012, 33112, 51552,
                                                 19742, 20012, 35182, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 53442, 0, 3, 48852, 49392, 33832, 52092,
                                                 20552, 20882, 35632, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 54267, 0, 3, 50472, 51012, 34732, 52767,
                                                 21542, 21872, 36182, ncols, alpha, beta, p);

            simdgeo::geom_l_x(buffer, 55092, 50472, 54267, 1, 15, ncols, alpha);

            simdgeo::geom_l_y(buffer, 55767, 50472, 54267, 1, 15, ncols, alpha);

            simdgeo::geom_l_z(buffer, 56442, 50472, 54267, 1, 15, ncols, alpha);

            simdgeo::geom_l_x(buffer, 57117, 48852, 53442, 1, 15, ncols, alpha);

            simdgeo::geom_l_y(buffer, 57792, 48852, 53442, 1, 15, ncols, alpha);

            simdgeo::geom_l_z(buffer, 58467, 48852, 53442, 1, 15, ncols, alpha);

            simdfunc::contract_primitives(buffer, 59142, 57117, 2025, ncols);

            simdfunc::contract_primitives(buffer, 61167, 55092, 2025, ncols);
        }
    }

    simdtrf::transform_g_inner(buffer, 63192, 61167, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 63192, 9, nmax);

    simdtrf::transform_g_inner(buffer, 63192, 61842, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 153 * nvalues, nvalues, buffer, 63192, 9, nmax);

    simdtrf::transform_g_inner(buffer, 63192, 62517, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 306 * nvalues, nvalues, buffer, 63192, 9, nmax);

    simdtrf::transform_g_inner(buffer, 63192, 59142, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 459 * nvalues, nvalues, buffer, 63192, 9, nmax);

    simdtrf::transform_g_inner(buffer, 63192, 59817, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 612 * nvalues, nvalues, buffer, 63192, 9, nmax);

    simdtrf::transform_g_inner(buffer, 63192, 60492, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 765 * nvalues, nvalues, buffer, 63192, 9, nmax);
}

}  // namespace simdt2ceri
