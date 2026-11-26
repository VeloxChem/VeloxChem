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


#include "SimdElectronRepulsionGeom10RsRecLI.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGI.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHI.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecII.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKI.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLI.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecMD.hpp"
#include "SimdElectronRepulsionVrrRecMF.hpp"
#include "SimdElectronRepulsionVrrRecMG.hpp"
#include "SimdElectronRepulsionVrrRecMH.hpp"
#include "SimdElectronRepulsionVrrRecMI.hpp"
#include "SimdElectronRepulsionVrrRecMP.hpp"
#include "SimdElectronRepulsionVrrRecMS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryL1.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_li_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_li_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 165689, 157544, 7560, nvalues);

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
                                                9, 10, 11, 12, 13, 14, 15}, ncols, fj, mu,
                                                omega);

            simdfunc::compute_boys_function(buffer, coordinates, 22, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13, 14, 15}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 74, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 77, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 80, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 83, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 86, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 89, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 92, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 95, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 98, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 101, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 104, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 107, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 110, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 113, 0, 33, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 116, 0, 34, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 119, 0, 35, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 122, 0, 36, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 125, 0, 37, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 7, 8, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 8, 9, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 9, 10, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 10, 11, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 11, 12, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 158, 0, 12, 13, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 164, 0, 13, 14, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 170, 0, 14, 15, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 176, 0, 15, 16, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 182, 0, 16, 17, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 188, 0, 17, 18, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 194, 0, 18, 19, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 200, 0, 19, 20, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 206, 0, 23, 24, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 212, 0, 24, 25, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 218, 0, 25, 26, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 224, 0, 26, 27, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 230, 0, 27, 28, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 236, 0, 28, 29, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 242, 0, 29, 30, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 248, 0, 30, 31, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 254, 0, 31, 32, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 260, 0, 32, 33, 116, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 266, 0, 33, 34, 119, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 272, 0, 34, 35, 122, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 278, 0, 35, 36, 125, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 284, 0, 38, 41, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 294, 0, 41, 44, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 304, 0, 44, 47, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 314, 0, 47, 50, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 324, 0, 50, 53, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 334, 0, 53, 56, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 344, 0, 56, 59, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 354, 0, 59, 62, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 364, 0, 62, 65, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 374, 0, 65, 68, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 384, 0, 68, 71, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 394, 0, 71, 74, 194, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 404, 0, 74, 77, 200, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 414, 0, 83, 86, 206, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 424, 0, 86, 89, 212, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 434, 0, 89, 92, 218, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 444, 0, 92, 95, 224, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 454, 0, 95, 98, 230, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 464, 0, 98, 101, 236, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 474, 0, 101, 104, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 484, 0, 104, 107, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 494, 0, 107, 110, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 504, 0, 110, 113, 260, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 514, 0, 113, 116, 266, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 524, 0, 116, 119, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 534, 0, 119, 122, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 544, 0, 128, 134, 304, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 559, 0, 134, 140, 314, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 574, 0, 140, 146, 324, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 589, 0, 146, 152, 334, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 604, 0, 152, 158, 344, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 619, 0, 158, 164, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 634, 0, 164, 170, 364, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 649, 0, 170, 176, 374, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 664, 0, 176, 182, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 679, 0, 182, 188, 394, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 694, 0, 188, 194, 404, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 709, 0, 206, 212, 434, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 724, 0, 212, 218, 444, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 739, 0, 218, 224, 454, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 754, 0, 224, 230, 464, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 769, 0, 230, 236, 474, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 784, 0, 236, 242, 484, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 799, 0, 242, 248, 494, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 814, 0, 248, 254, 504, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 829, 0, 254, 260, 514, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 844, 0, 260, 266, 524, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 859, 0, 266, 272, 534, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 874, 0, 284, 294, 544, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 895, 0, 294, 304, 559, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 916, 0, 304, 314, 574, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 937, 0, 314, 324, 589, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 958, 0, 324, 334, 604, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 979, 0, 334, 344, 619, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1000, 0, 344, 354, 634, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1021, 0, 354, 364, 649, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1042, 0, 364, 374, 664, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1063, 0, 374, 384, 679, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1084, 0, 384, 394, 694, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1105, 0, 414, 424, 709, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1126, 0, 424, 434, 724, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1147, 0, 434, 444, 739, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1168, 0, 444, 454, 754, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1189, 0, 454, 464, 769, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1210, 0, 464, 474, 784, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1231, 0, 474, 484, 799, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1252, 0, 484, 494, 814, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1273, 0, 494, 504, 829, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1294, 0, 504, 514, 844, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1315, 0, 514, 524, 859, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1336, 0, 544, 559, 916, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1364, 0, 559, 574, 937, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1392, 0, 574, 589, 958, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1420, 0, 589, 604, 979, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1448, 0, 604, 619, 1000, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1476, 0, 619, 634, 1021, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1504, 0, 634, 649, 1042, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1532, 0, 649, 664, 1063, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1560, 0, 664, 679, 1084, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1588, 0, 709, 724, 1147, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1616, 0, 724, 739, 1168, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1644, 0, 739, 754, 1189, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1672, 0, 754, 769, 1210, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1700, 0, 769, 784, 1231, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1728, 0, 784, 799, 1252, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1756, 0, 799, 814, 1273, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1784, 0, 814, 829, 1294, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1812, 0, 829, 844, 1315, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1840, 0, 874, 895, 1336, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1876, 0, 895, 916, 1364, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1912, 0, 916, 937, 1392, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1948, 0, 937, 958, 1420, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1984, 0, 958, 979, 1448, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2020, 0, 979, 1000, 1476, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2056, 0, 1000, 1021, 1504, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2092, 0, 1021, 1042, 1532, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2128, 0, 1042, 1063, 1560, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2164, 0, 1105, 1126, 1588, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2200, 0, 1126, 1147, 1616, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2236, 0, 1147, 1168, 1644, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2272, 0, 1168, 1189, 1672, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2308, 0, 1189, 1210, 1700, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2344, 0, 1210, 1231, 1728, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2380, 0, 1231, 1252, 1756, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2416, 0, 1252, 1273, 1784, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2452, 0, 1273, 1294, 1812, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2488, 0, 1336, 1364, 1912, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2533, 0, 1364, 1392, 1948, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2578, 0, 1392, 1420, 1984, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2623, 0, 1420, 1448, 2020, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2668, 0, 1448, 1476, 2056, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2713, 0, 1476, 1504, 2092, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2758, 0, 1504, 1532, 2128, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2803, 0, 1588, 1616, 2236, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2848, 0, 1616, 1644, 2272, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2893, 0, 1644, 1672, 2308, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2938, 0, 1672, 1700, 2344, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2983, 0, 1700, 1728, 2380, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3028, 0, 1728, 1756, 2416, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 3073, 0, 1756, 1784, 2452, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3118, 0, 1840, 1876, 2488, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3173, 0, 1876, 1912, 2533, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3228, 0, 1912, 1948, 2578, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3283, 0, 1948, 1984, 2623, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3338, 0, 1984, 2020, 2668, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3393, 0, 2020, 2056, 2713, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3448, 0, 2056, 2092, 2758, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3503, 0, 2164, 2200, 2803, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3558, 0, 2200, 2236, 2848, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3613, 0, 2236, 2272, 2893, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3668, 0, 2272, 2308, 2938, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3723, 0, 2308, 2344, 2983, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3778, 0, 2344, 2380, 3028, ncols, alpha,
                                                 beta, p);

            compute_prim_ms_electron_repulsion_0(buffer, 3833, 0, 2380, 2416, 3073, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 3888, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3891, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3894, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3897, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3900, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3903, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3906, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3909, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3912, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3915, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3918, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3921, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3924, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3927, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3930, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3933, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3936, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3939, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3942, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3945, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3948, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3951, 3, 35, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3954, 3, 36, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 3957, 3, 37, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 3960, 3, 9, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3969, 3, 10, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3978, 3, 11, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3987, 3, 12, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3996, 3, 13, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4005, 3, 14, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4014, 3, 15, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4023, 3, 16, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4032, 3, 17, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4041, 3, 18, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4050, 3, 19, 77, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4059, 3, 20, 80, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4068, 3, 25, 92, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4077, 3, 26, 95, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4086, 3, 27, 98, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4095, 3, 28, 101, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4104, 3, 29, 104, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4113, 3, 30, 107, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4122, 3, 31, 110, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4131, 3, 32, 113, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4140, 3, 33, 116, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4149, 3, 34, 119, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4158, 3, 35, 122, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 4167, 3, 36, 125, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4176, 0, 3, 44, 3960, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4194, 0, 3, 47, 3969, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4212, 0, 3, 50, 3978, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4230, 0, 3, 53, 3987, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4248, 0, 3, 56, 3996, 158, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4266, 0, 3, 59, 4005, 164, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4284, 0, 3, 62, 4014, 170, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4302, 0, 3, 65, 4023, 176, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4320, 0, 3, 68, 4032, 182, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4338, 0, 3, 71, 4041, 188, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4356, 0, 3, 74, 4050, 194, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4374, 0, 3, 77, 4059, 200, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4392, 0, 3, 89, 4068, 212, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4410, 0, 3, 92, 4077, 218, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4428, 0, 3, 95, 4086, 224, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4446, 0, 3, 98, 4095, 230, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4464, 0, 3, 101, 4104, 236, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4482, 0, 3, 104, 4113, 242, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4500, 0, 3, 107, 4122, 248, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4518, 0, 3, 110, 4131, 254, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4536, 0, 3, 113, 4140, 260, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4554, 0, 3, 116, 4149, 266, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4572, 0, 3, 119, 4158, 272, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 4590, 0, 3, 122, 4167, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4608, 0, 3, 134, 4194, 304, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4638, 0, 3, 140, 4212, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4668, 0, 3, 146, 4230, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4698, 0, 3, 152, 4248, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4728, 0, 3, 158, 4266, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4758, 0, 3, 164, 4284, 354, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4788, 0, 3, 170, 4302, 364, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4818, 0, 3, 176, 4320, 374, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4848, 0, 3, 182, 4338, 384, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4878, 0, 3, 188, 4356, 394, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4908, 0, 3, 194, 4374, 404, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4938, 0, 3, 212, 4410, 434, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4968, 0, 3, 218, 4428, 444, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4998, 0, 3, 224, 4446, 454, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5028, 0, 3, 230, 4464, 464, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5058, 0, 3, 236, 4482, 474, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5088, 0, 3, 242, 4500, 484, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5118, 0, 3, 248, 4518, 494, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5148, 0, 3, 254, 4536, 504, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5178, 0, 3, 260, 4554, 514, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5208, 0, 3, 266, 4572, 524, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 5238, 0, 3, 272, 4590, 534, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5268, 0, 3, 304, 4638, 559, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5313, 0, 3, 314, 4668, 574, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5358, 0, 3, 324, 4698, 589, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5403, 0, 3, 334, 4728, 604, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5448, 0, 3, 344, 4758, 619, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5493, 0, 3, 354, 4788, 634, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5538, 0, 3, 364, 4818, 649, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5583, 0, 3, 374, 4848, 664, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5628, 0, 3, 384, 4878, 679, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5673, 0, 3, 394, 4908, 694, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5718, 0, 3, 434, 4968, 724, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5763, 0, 3, 444, 4998, 739, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5808, 0, 3, 454, 5028, 754, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5853, 0, 3, 464, 5058, 769, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5898, 0, 3, 474, 5088, 784, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5943, 0, 3, 484, 5118, 799, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 5988, 0, 3, 494, 5148, 814, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6033, 0, 3, 504, 5178, 829, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6078, 0, 3, 514, 5208, 844, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 6123, 0, 3, 524, 5238, 859, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6168, 0, 3, 559, 5313, 916, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6231, 0, 3, 574, 5358, 937, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6294, 0, 3, 589, 5403, 958, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6357, 0, 3, 604, 5448, 979, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6420, 0, 3, 619, 5493, 1000, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6483, 0, 3, 634, 5538, 1021, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6546, 0, 3, 649, 5583, 1042, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6609, 0, 3, 664, 5628, 1063, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6672, 0, 3, 679, 5673, 1084, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6735, 0, 3, 724, 5763, 1147, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6798, 0, 3, 739, 5808, 1168, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6861, 0, 3, 754, 5853, 1189, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6924, 0, 3, 769, 5898, 1210, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 6987, 0, 3, 784, 5943, 1231, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7050, 0, 3, 799, 5988, 1252, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7113, 0, 3, 814, 6033, 1273, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7176, 0, 3, 829, 6078, 1294, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 7239, 0, 3, 844, 6123, 1315, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7302, 0, 3, 916, 6231, 1364, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7386, 0, 3, 937, 6294, 1392, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7470, 0, 3, 958, 6357, 1420, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7554, 0, 3, 979, 6420, 1448, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 7638, 0, 3, 1000, 6483, 1476, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7722, 0, 3, 1021, 6546, 1504, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7806, 0, 3, 1042, 6609, 1532, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7890, 0, 3, 1063, 6672, 1560, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7974, 0, 3, 1147, 6798, 1616, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8058, 0, 3, 1168, 6861, 1644, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8142, 0, 3, 1189, 6924, 1672, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8226, 0, 3, 1210, 6987, 1700, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8310, 0, 3, 1231, 7050, 1728, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8394, 0, 3, 1252, 7113, 1756, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8478, 0, 3, 1273, 7176, 1784, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 8562, 0, 3, 1294, 7239, 1812, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8646, 0, 3, 1364, 7386, 1912, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8754, 0, 3, 1392, 7470, 1948, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8862, 0, 3, 1420, 7554, 1984, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8970, 0, 3, 1448, 7638, 2020, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9078, 0, 3, 1476, 7722, 2056, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9186, 0, 3, 1504, 7806, 2092, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9294, 0, 3, 1532, 7890, 2128, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9402, 0, 3, 1616, 8058, 2236, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9510, 0, 3, 1644, 8142, 2272, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9618, 0, 3, 1672, 8226, 2308, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9726, 0, 3, 1700, 8310, 2344, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9834, 0, 3, 1728, 8394, 2380, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 9942, 0, 3, 1756, 8478, 2416, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 10050, 0, 3, 1784, 8562, 2452, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10158, 0, 3, 1912, 8754, 2533, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10293, 0, 3, 1948, 8862, 2578, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10428, 0, 3, 1984, 8970, 2623, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10563, 0, 3, 2020, 9078, 2668, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10698, 0, 3, 2056, 9186, 2713, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10833, 0, 3, 2092, 9294, 2758, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 10968, 0, 3, 2236, 9510, 2848, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11103, 0, 3, 2272, 9618, 2893, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11238, 0, 3, 2308, 9726, 2938, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11373, 0, 3, 2344, 9834, 2983, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11508, 0, 3, 2380, 9942, 3028, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 11643, 0, 3, 2416, 10050, 3073, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 11778, 0, 3, 2533, 10293, 3228, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 11943, 0, 3, 2578, 10428, 3283, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 12108, 0, 3, 2623, 10563, 3338, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 12273, 0, 3, 2668, 10698, 3393, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 12438, 0, 3, 2713, 10833, 3448, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 12603, 0, 3, 2848, 11103, 3613, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 12768, 0, 3, 2893, 11238, 3668, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 12933, 0, 3, 2938, 11373, 3723, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 13098, 0, 3, 2983, 11508, 3778, ncols,
                                                 p);

            compute_prim_mp_electron_repulsion_0(buffer, 13263, 0, 3, 3028, 11643, 3833, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 13428, 3, 9, 10, 3891, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13434, 3, 10, 11, 3894, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13440, 3, 11, 12, 3897, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13446, 3, 12, 13, 3900, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13452, 3, 13, 14, 3903, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13458, 3, 14, 15, 3906, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13464, 3, 15, 16, 3909, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13470, 3, 16, 17, 3912, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13476, 3, 17, 18, 3915, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13482, 3, 18, 19, 3918, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13488, 3, 19, 20, 3921, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13494, 3, 25, 26, 3927, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13500, 3, 26, 27, 3930, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13506, 3, 27, 28, 3933, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13512, 3, 28, 29, 3936, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13518, 3, 29, 30, 3939, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13524, 3, 30, 31, 3942, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13530, 3, 31, 32, 3945, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13536, 3, 32, 33, 3948, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13542, 3, 33, 34, 3951, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13548, 3, 34, 35, 3954, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 13554, 3, 35, 36, 3957, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 13560, 0, 3, 3888, 13428, 3969, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13578, 0, 3, 3891, 13434, 3978, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13596, 0, 3, 3894, 13440, 3987, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13614, 0, 3, 3897, 13446, 3996, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13632, 0, 3, 3900, 13452, 4005, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13650, 0, 3, 3903, 13458, 4014, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13668, 0, 3, 3906, 13464, 4023, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13686, 0, 3, 3909, 13470, 4032, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13704, 0, 3, 3912, 13476, 4041, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13722, 0, 3, 3915, 13482, 4050, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13740, 0, 3, 3918, 13488, 4059, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13758, 0, 3, 3924, 13494, 4077, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13776, 0, 3, 3927, 13500, 4086, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13794, 0, 3, 3930, 13506, 4095, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13812, 0, 3, 3933, 13512, 4104, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13830, 0, 3, 3936, 13518, 4113, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13848, 0, 3, 3939, 13524, 4122, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13866, 0, 3, 3942, 13530, 4131, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13884, 0, 3, 3945, 13536, 4140, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13902, 0, 3, 3948, 13542, 4149, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13920, 0, 3, 3951, 13548, 4158, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 13938, 0, 3, 3954, 13554, 4167, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 13956, 0, 3, 3960, 13560, 128, 134,
                                                 4194, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 13992, 0, 3, 3969, 13578, 134, 140,
                                                 4212, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14028, 0, 3, 3978, 13596, 140, 146,
                                                 4230, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14064, 0, 3, 3987, 13614, 146, 152,
                                                 4248, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14100, 0, 3, 3996, 13632, 152, 158,
                                                 4266, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14136, 0, 3, 4005, 13650, 158, 164,
                                                 4284, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14172, 0, 3, 4014, 13668, 164, 170,
                                                 4302, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14208, 0, 3, 4023, 13686, 170, 176,
                                                 4320, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14244, 0, 3, 4032, 13704, 176, 182,
                                                 4338, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14280, 0, 3, 4041, 13722, 182, 188,
                                                 4356, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14316, 0, 3, 4050, 13740, 188, 194,
                                                 4374, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14352, 0, 3, 4068, 13758, 206, 212,
                                                 4410, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14388, 0, 3, 4077, 13776, 212, 218,
                                                 4428, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14424, 0, 3, 4086, 13794, 218, 224,
                                                 4446, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14460, 0, 3, 4095, 13812, 224, 230,
                                                 4464, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14496, 0, 3, 4104, 13830, 230, 236,
                                                 4482, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14532, 0, 3, 4113, 13848, 236, 242,
                                                 4500, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14568, 0, 3, 4122, 13866, 242, 248,
                                                 4518, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14604, 0, 3, 4131, 13884, 248, 254,
                                                 4536, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14640, 0, 3, 4140, 13902, 254, 260,
                                                 4554, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14676, 0, 3, 4149, 13920, 260, 266,
                                                 4572, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 14712, 0, 3, 4158, 13938, 266, 272,
                                                 4590, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14748, 0, 3, 4176, 13956, 284, 294,
                                                 4608, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14808, 0, 3, 4194, 13992, 294, 304,
                                                 4638, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14868, 0, 3, 4212, 14028, 304, 314,
                                                 4668, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14928, 0, 3, 4230, 14064, 314, 324,
                                                 4698, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 14988, 0, 3, 4248, 14100, 324, 334,
                                                 4728, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15048, 0, 3, 4266, 14136, 334, 344,
                                                 4758, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15108, 0, 3, 4284, 14172, 344, 354,
                                                 4788, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15168, 0, 3, 4302, 14208, 354, 364,
                                                 4818, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15228, 0, 3, 4320, 14244, 364, 374,
                                                 4848, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15288, 0, 3, 4338, 14280, 374, 384,
                                                 4878, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15348, 0, 3, 4356, 14316, 384, 394,
                                                 4908, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15408, 0, 3, 4392, 14352, 414, 424,
                                                 4938, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15468, 0, 3, 4410, 14388, 424, 434,
                                                 4968, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15528, 0, 3, 4428, 14424, 434, 444,
                                                 4998, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15588, 0, 3, 4446, 14460, 444, 454,
                                                 5028, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15648, 0, 3, 4464, 14496, 454, 464,
                                                 5058, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15708, 0, 3, 4482, 14532, 464, 474,
                                                 5088, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15768, 0, 3, 4500, 14568, 474, 484,
                                                 5118, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15828, 0, 3, 4518, 14604, 484, 494,
                                                 5148, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15888, 0, 3, 4536, 14640, 494, 504,
                                                 5178, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 15948, 0, 3, 4554, 14676, 504, 514,
                                                 5208, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 16008, 0, 3, 4572, 14712, 514, 524,
                                                 5238, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16068, 0, 3, 13956, 13992, 4638, 14868,
                                                 544, 559, 5313, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16158, 0, 3, 13992, 14028, 4668, 14928,
                                                 559, 574, 5358, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16248, 0, 3, 14028, 14064, 4698, 14988,
                                                 574, 589, 5403, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16338, 0, 3, 14064, 14100, 4728, 15048,
                                                 589, 604, 5448, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16428, 0, 3, 14100, 14136, 4758, 15108,
                                                 604, 619, 5493, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16518, 0, 3, 14136, 14172, 4788, 15168,
                                                 619, 634, 5538, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16608, 0, 3, 14172, 14208, 4818, 15228,
                                                 634, 649, 5583, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16698, 0, 3, 14208, 14244, 4848, 15288,
                                                 649, 664, 5628, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16788, 0, 3, 14244, 14280, 4878, 15348,
                                                 664, 679, 5673, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16878, 0, 3, 14352, 14388, 4968, 15528,
                                                 709, 724, 5763, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 16968, 0, 3, 14388, 14424, 4998, 15588,
                                                 724, 739, 5808, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17058, 0, 3, 14424, 14460, 5028, 15648,
                                                 739, 754, 5853, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17148, 0, 3, 14460, 14496, 5058, 15708,
                                                 754, 769, 5898, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17238, 0, 3, 14496, 14532, 5088, 15768,
                                                 769, 784, 5943, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17328, 0, 3, 14532, 14568, 5118, 15828,
                                                 784, 799, 5988, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17418, 0, 3, 14568, 14604, 5148, 15888,
                                                 799, 814, 6033, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17508, 0, 3, 14604, 14640, 5178, 15948,
                                                 814, 829, 6078, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 17598, 0, 3, 14640, 14676, 5208, 16008,
                                                 829, 844, 6123, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17688, 0, 3, 14748, 14808, 5268, 16068,
                                                 874, 895, 6168, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17814, 0, 3, 14808, 14868, 5313, 16158,
                                                 895, 916, 6231, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 17940, 0, 3, 14868, 14928, 5358, 16248,
                                                 916, 937, 6294, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18066, 0, 3, 14928, 14988, 5403, 16338,
                                                 937, 958, 6357, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18192, 0, 3, 14988, 15048, 5448, 16428,
                                                 958, 979, 6420, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18318, 0, 3, 15048, 15108, 5493, 16518,
                                                 979, 1000, 6483, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18444, 0, 3, 15108, 15168, 5538, 16608,
                                                 1000, 1021, 6546, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18570, 0, 3, 15168, 15228, 5583, 16698,
                                                 1021, 1042, 6609, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18696, 0, 3, 15228, 15288, 5628, 16788,
                                                 1042, 1063, 6672, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18822, 0, 3, 15408, 15468, 5718, 16878,
                                                 1105, 1126, 6735, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 18948, 0, 3, 15468, 15528, 5763, 16968,
                                                 1126, 1147, 6798, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19074, 0, 3, 15528, 15588, 5808, 17058,
                                                 1147, 1168, 6861, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19200, 0, 3, 15588, 15648, 5853, 17148,
                                                 1168, 1189, 6924, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19326, 0, 3, 15648, 15708, 5898, 17238,
                                                 1189, 1210, 6987, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19452, 0, 3, 15708, 15768, 5943, 17328,
                                                 1210, 1231, 7050, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19578, 0, 3, 15768, 15828, 5988, 17418,
                                                 1231, 1252, 7113, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19704, 0, 3, 15828, 15888, 6033, 17508,
                                                 1252, 1273, 7176, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 19830, 0, 3, 15888, 15948, 6078, 17598,
                                                 1273, 1294, 7239, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 19956, 0, 3, 16068, 16158, 6231, 17940,
                                                 1336, 1364, 7386, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20124, 0, 3, 16158, 16248, 6294, 18066,
                                                 1364, 1392, 7470, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20292, 0, 3, 16248, 16338, 6357, 18192,
                                                 1392, 1420, 7554, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20460, 0, 3, 16338, 16428, 6420, 18318,
                                                 1420, 1448, 7638, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20628, 0, 3, 16428, 16518, 6483, 18444,
                                                 1448, 1476, 7722, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20796, 0, 3, 16518, 16608, 6546, 18570,
                                                 1476, 1504, 7806, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 20964, 0, 3, 16608, 16698, 6609, 18696,
                                                 1504, 1532, 7890, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21132, 0, 3, 16878, 16968, 6798, 19074,
                                                 1588, 1616, 8058, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21300, 0, 3, 16968, 17058, 6861, 19200,
                                                 1616, 1644, 8142, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21468, 0, 3, 17058, 17148, 6924, 19326,
                                                 1644, 1672, 8226, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21636, 0, 3, 17148, 17238, 6987, 19452,
                                                 1672, 1700, 8310, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21804, 0, 3, 17238, 17328, 7050, 19578,
                                                 1700, 1728, 8394, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 21972, 0, 3, 17328, 17418, 7113, 19704,
                                                 1728, 1756, 8478, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 22140, 0, 3, 17418, 17508, 7176, 19830,
                                                 1756, 1784, 8562, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22308, 0, 3, 17688, 17814, 7302, 19956,
                                                 1840, 1876, 8646, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22524, 0, 3, 17814, 17940, 7386, 20124,
                                                 1876, 1912, 8754, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22740, 0, 3, 17940, 18066, 7470, 20292,
                                                 1912, 1948, 8862, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 22956, 0, 3, 18066, 18192, 7554, 20460,
                                                 1948, 1984, 8970, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 23172, 0, 3, 18192, 18318, 7638, 20628,
                                                 1984, 2020, 9078, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 23388, 0, 3, 18318, 18444, 7722, 20796,
                                                 2020, 2056, 9186, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 23604, 0, 3, 18444, 18570, 7806, 20964,
                                                 2056, 2092, 9294, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 23820, 0, 3, 18822, 18948, 7974, 21132,
                                                 2164, 2200, 9402, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24036, 0, 3, 18948, 19074, 8058, 21300,
                                                 2200, 2236, 9510, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24252, 0, 3, 19074, 19200, 8142, 21468,
                                                 2236, 2272, 9618, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24468, 0, 3, 19200, 19326, 8226, 21636,
                                                 2272, 2308, 9726, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24684, 0, 3, 19326, 19452, 8310, 21804,
                                                 2308, 2344, 9834, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 24900, 0, 3, 19452, 19578, 8394, 21972,
                                                 2344, 2380, 9942, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 25116, 0, 3, 19578, 19704, 8478, 22140,
                                                 2380, 2416, 10050, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 25332, 0, 3, 19956, 20124, 8754, 22740,
                                                 2488, 2533, 10293, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 25602, 0, 3, 20124, 20292, 8862, 22956,
                                                 2533, 2578, 10428, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 25872, 0, 3, 20292, 20460, 8970, 23172,
                                                 2578, 2623, 10563, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 26142, 0, 3, 20460, 20628, 9078, 23388,
                                                 2623, 2668, 10698, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 26412, 0, 3, 20628, 20796, 9186, 23604,
                                                 2668, 2713, 10833, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 26682, 0, 3, 21132, 21300, 9510, 24252,
                                                 2803, 2848, 11103, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 26952, 0, 3, 21300, 21468, 9618, 24468,
                                                 2848, 2893, 11238, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 27222, 0, 3, 21468, 21636, 9726, 24684,
                                                 2893, 2938, 11373, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 27492, 0, 3, 21636, 21804, 9834, 24900,
                                                 2938, 2983, 11508, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 27762, 0, 3, 21804, 21972, 9942, 25116,
                                                 2983, 3028, 11643, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 28032, 0, 3, 22308, 22524, 10158, 25332,
                                                 3118, 3173, 11778, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 28362, 0, 3, 22524, 22740, 10293, 25602,
                                                 3173, 3228, 11943, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 28692, 0, 3, 22740, 22956, 10428, 25872,
                                                 3228, 3283, 12108, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 29022, 0, 3, 22956, 23172, 10563, 26142,
                                                 3283, 3338, 12273, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 29352, 0, 3, 23172, 23388, 10698, 26412,
                                                 3338, 3393, 12438, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 29682, 0, 3, 23820, 24036, 10968, 26682,
                                                 3503, 3558, 12603, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 30012, 0, 3, 24036, 24252, 11103, 26952,
                                                 3558, 3613, 12768, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 30342, 0, 3, 24252, 24468, 11238, 27222,
                                                 3613, 3668, 12933, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 30672, 0, 3, 24468, 24684, 11373, 27492,
                                                 3668, 3723, 13098, ncols, alpha, beta, p);

            compute_prim_md_electron_repulsion_0(buffer, 31002, 0, 3, 24684, 24900, 11508, 27762,
                                                 3723, 3778, 13263, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31332, 3, 3888, 3891, 13434, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31342, 3, 3891, 3894, 13440, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31352, 3, 3894, 3897, 13446, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31362, 3, 3897, 3900, 13452, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31372, 3, 3900, 3903, 13458, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31382, 3, 3903, 3906, 13464, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31392, 3, 3906, 3909, 13470, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31402, 3, 3909, 3912, 13476, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31412, 3, 3912, 3915, 13482, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31422, 3, 3915, 3918, 13488, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31432, 3, 3924, 3927, 13500, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31442, 3, 3927, 3930, 13506, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31452, 3, 3930, 3933, 13512, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31462, 3, 3933, 3936, 13518, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31472, 3, 3936, 3939, 13524, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31482, 3, 3939, 3942, 13530, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31492, 3, 3942, 3945, 13536, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31502, 3, 3945, 3948, 13542, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31512, 3, 3948, 3951, 13548, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 31522, 3, 3951, 3954, 13554, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 31532, 0, 3, 13428, 31332, 13578, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31562, 0, 3, 13434, 31342, 13596, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31592, 0, 3, 13440, 31352, 13614, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31622, 0, 3, 13446, 31362, 13632, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31652, 0, 3, 13452, 31372, 13650, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31682, 0, 3, 13458, 31382, 13668, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31712, 0, 3, 13464, 31392, 13686, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31742, 0, 3, 13470, 31402, 13704, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31772, 0, 3, 13476, 31412, 13722, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31802, 0, 3, 13482, 31422, 13740, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31832, 0, 3, 13494, 31432, 13776, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31862, 0, 3, 13500, 31442, 13794, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31892, 0, 3, 13506, 31452, 13812, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31922, 0, 3, 13512, 31462, 13830, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31952, 0, 3, 13518, 31472, 13848, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 31982, 0, 3, 13524, 31482, 13866, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 32012, 0, 3, 13530, 31492, 13884, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 32042, 0, 3, 13536, 31502, 13902, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 32072, 0, 3, 13542, 31512, 13920, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 32102, 0, 3, 13548, 31522, 13938, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 32132, 0, 3, 13560, 31532, 4176, 4194,
                                                 13992, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32192, 0, 3, 13578, 31562, 4194, 4212,
                                                 14028, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32252, 0, 3, 13596, 31592, 4212, 4230,
                                                 14064, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32312, 0, 3, 13614, 31622, 4230, 4248,
                                                 14100, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32372, 0, 3, 13632, 31652, 4248, 4266,
                                                 14136, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32432, 0, 3, 13650, 31682, 4266, 4284,
                                                 14172, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32492, 0, 3, 13668, 31712, 4284, 4302,
                                                 14208, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32552, 0, 3, 13686, 31742, 4302, 4320,
                                                 14244, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32612, 0, 3, 13704, 31772, 4320, 4338,
                                                 14280, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32672, 0, 3, 13722, 31802, 4338, 4356,
                                                 14316, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32732, 0, 3, 13758, 31832, 4392, 4410,
                                                 14388, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32792, 0, 3, 13776, 31862, 4410, 4428,
                                                 14424, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32852, 0, 3, 13794, 31892, 4428, 4446,
                                                 14460, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32912, 0, 3, 13812, 31922, 4446, 4464,
                                                 14496, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 32972, 0, 3, 13830, 31952, 4464, 4482,
                                                 14532, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 33032, 0, 3, 13848, 31982, 4482, 4500,
                                                 14568, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 33092, 0, 3, 13866, 32012, 4500, 4518,
                                                 14604, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 33152, 0, 3, 13884, 32042, 4518, 4536,
                                                 14640, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 33212, 0, 3, 13902, 32072, 4536, 4554,
                                                 14676, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 33272, 0, 3, 13920, 32102, 4554, 4572,
                                                 14712, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33332, 0, 3, 13992, 32192, 4608, 4638,
                                                 14868, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33432, 0, 3, 14028, 32252, 4638, 4668,
                                                 14928, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33532, 0, 3, 14064, 32312, 4668, 4698,
                                                 14988, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33632, 0, 3, 14100, 32372, 4698, 4728,
                                                 15048, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33732, 0, 3, 14136, 32432, 4728, 4758,
                                                 15108, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33832, 0, 3, 14172, 32492, 4758, 4788,
                                                 15168, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 33932, 0, 3, 14208, 32552, 4788, 4818,
                                                 15228, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34032, 0, 3, 14244, 32612, 4818, 4848,
                                                 15288, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34132, 0, 3, 14280, 32672, 4848, 4878,
                                                 15348, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34232, 0, 3, 14388, 32792, 4938, 4968,
                                                 15528, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34332, 0, 3, 14424, 32852, 4968, 4998,
                                                 15588, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34432, 0, 3, 14460, 32912, 4998, 5028,
                                                 15648, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34532, 0, 3, 14496, 32972, 5028, 5058,
                                                 15708, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34632, 0, 3, 14532, 33032, 5058, 5088,
                                                 15768, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34732, 0, 3, 14568, 33092, 5088, 5118,
                                                 15828, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34832, 0, 3, 14604, 33152, 5118, 5148,
                                                 15888, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 34932, 0, 3, 14640, 33212, 5148, 5178,
                                                 15948, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 35032, 0, 3, 14676, 33272, 5178, 5208,
                                                 16008, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35132, 0, 3, 32132, 32192, 14868, 33432,
                                                 5268, 5313, 16158, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35282, 0, 3, 32192, 32252, 14928, 33532,
                                                 5313, 5358, 16248, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35432, 0, 3, 32252, 32312, 14988, 33632,
                                                 5358, 5403, 16338, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35582, 0, 3, 32312, 32372, 15048, 33732,
                                                 5403, 5448, 16428, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35732, 0, 3, 32372, 32432, 15108, 33832,
                                                 5448, 5493, 16518, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 35882, 0, 3, 32432, 32492, 15168, 33932,
                                                 5493, 5538, 16608, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36032, 0, 3, 32492, 32552, 15228, 34032,
                                                 5538, 5583, 16698, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36182, 0, 3, 32552, 32612, 15288, 34132,
                                                 5583, 5628, 16788, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36332, 0, 3, 32732, 32792, 15528, 34332,
                                                 5718, 5763, 16968, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36482, 0, 3, 32792, 32852, 15588, 34432,
                                                 5763, 5808, 17058, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36632, 0, 3, 32852, 32912, 15648, 34532,
                                                 5808, 5853, 17148, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36782, 0, 3, 32912, 32972, 15708, 34632,
                                                 5853, 5898, 17238, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 36932, 0, 3, 32972, 33032, 15768, 34732,
                                                 5898, 5943, 17328, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 37082, 0, 3, 33032, 33092, 15828, 34832,
                                                 5943, 5988, 17418, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 37232, 0, 3, 33092, 33152, 15888, 34932,
                                                 5988, 6033, 17508, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 37382, 0, 3, 33152, 33212, 15948, 35032,
                                                 6033, 6078, 17598, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 37532, 0, 3, 33332, 33432, 16158, 35282,
                                                 6168, 6231, 17940, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 37742, 0, 3, 33432, 33532, 16248, 35432,
                                                 6231, 6294, 18066, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 37952, 0, 3, 33532, 33632, 16338, 35582,
                                                 6294, 6357, 18192, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38162, 0, 3, 33632, 33732, 16428, 35732,
                                                 6357, 6420, 18318, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38372, 0, 3, 33732, 33832, 16518, 35882,
                                                 6420, 6483, 18444, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38582, 0, 3, 33832, 33932, 16608, 36032,
                                                 6483, 6546, 18570, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 38792, 0, 3, 33932, 34032, 16698, 36182,
                                                 6546, 6609, 18696, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 39002, 0, 3, 34232, 34332, 16968, 36482,
                                                 6735, 6798, 19074, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 39212, 0, 3, 34332, 34432, 17058, 36632,
                                                 6798, 6861, 19200, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 39422, 0, 3, 34432, 34532, 17148, 36782,
                                                 6861, 6924, 19326, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 39632, 0, 3, 34532, 34632, 17238, 36932,
                                                 6924, 6987, 19452, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 39842, 0, 3, 34632, 34732, 17328, 37082,
                                                 6987, 7050, 19578, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 40052, 0, 3, 34732, 34832, 17418, 37232,
                                                 7050, 7113, 19704, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 40262, 0, 3, 34832, 34932, 17508, 37382,
                                                 7113, 7176, 19830, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 40472, 0, 3, 35132, 35282, 17940, 37742,
                                                 7302, 7386, 20124, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 40752, 0, 3, 35282, 35432, 18066, 37952,
                                                 7386, 7470, 20292, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 41032, 0, 3, 35432, 35582, 18192, 38162,
                                                 7470, 7554, 20460, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 41312, 0, 3, 35582, 35732, 18318, 38372,
                                                 7554, 7638, 20628, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 41592, 0, 3, 35732, 35882, 18444, 38582,
                                                 7638, 7722, 20796, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 41872, 0, 3, 35882, 36032, 18570, 38792,
                                                 7722, 7806, 20964, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 42152, 0, 3, 36332, 36482, 19074, 39212,
                                                 7974, 8058, 21300, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 42432, 0, 3, 36482, 36632, 19200, 39422,
                                                 8058, 8142, 21468, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 42712, 0, 3, 36632, 36782, 19326, 39632,
                                                 8142, 8226, 21636, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 42992, 0, 3, 36782, 36932, 19452, 39842,
                                                 8226, 8310, 21804, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 43272, 0, 3, 36932, 37082, 19578, 40052,
                                                 8310, 8394, 21972, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 43552, 0, 3, 37082, 37232, 19704, 40262,
                                                 8394, 8478, 22140, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 43832, 0, 3, 37532, 37742, 20124, 40752,
                                                 8646, 8754, 22740, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 44192, 0, 3, 37742, 37952, 20292, 41032,
                                                 8754, 8862, 22956, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 44552, 0, 3, 37952, 38162, 20460, 41312,
                                                 8862, 8970, 23172, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 44912, 0, 3, 38162, 38372, 20628, 41592,
                                                 8970, 9078, 23388, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 45272, 0, 3, 38372, 38582, 20796, 41872,
                                                 9078, 9186, 23604, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 45632, 0, 3, 39002, 39212, 21300, 42432,
                                                 9402, 9510, 24252, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 45992, 0, 3, 39212, 39422, 21468, 42712,
                                                 9510, 9618, 24468, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 46352, 0, 3, 39422, 39632, 21636, 42992,
                                                 9618, 9726, 24684, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 46712, 0, 3, 39632, 39842, 21804, 43272,
                                                 9726, 9834, 24900, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 47072, 0, 3, 39842, 40052, 21972, 43552,
                                                 9834, 9942, 25116, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 47432, 0, 3, 40472, 40752, 22740, 44192,
                                                 10158, 10293, 25602, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 47882, 0, 3, 40752, 41032, 22956, 44552,
                                                 10293, 10428, 25872, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 48332, 0, 3, 41032, 41312, 23172, 44912,
                                                 10428, 10563, 26142, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 48782, 0, 3, 41312, 41592, 23388, 45272,
                                                 10563, 10698, 26412, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 49232, 0, 3, 42152, 42432, 24252, 45992,
                                                 10968, 11103, 26952, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 49682, 0, 3, 42432, 42712, 24468, 46352,
                                                 11103, 11238, 27222, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 50132, 0, 3, 42712, 42992, 24684, 46712,
                                                 11238, 11373, 27492, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 50582, 0, 3, 42992, 43272, 24900, 47072,
                                                 11373, 11508, 27762, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 51032, 0, 3, 43832, 44192, 25602, 47882,
                                                 11778, 11943, 28692, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 51582, 0, 3, 44192, 44552, 25872, 48332,
                                                 11943, 12108, 29022, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 52132, 0, 3, 44552, 44912, 26142, 48782,
                                                 12108, 12273, 29352, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 52682, 0, 3, 45632, 45992, 26952, 49682,
                                                 12603, 12768, 30342, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 53232, 0, 3, 45992, 46352, 27222, 50132,
                                                 12768, 12933, 30672, ncols, alpha, beta, p);

            compute_prim_mf_electron_repulsion_0(buffer, 53782, 0, 3, 46352, 46712, 27492, 50582,
                                                 12933, 13098, 31002, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54332, 3, 13428, 13434, 31342, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54347, 3, 13434, 13440, 31352, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54362, 3, 13440, 13446, 31362, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54377, 3, 13446, 13452, 31372, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54392, 3, 13452, 13458, 31382, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54407, 3, 13458, 13464, 31392, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54422, 3, 13464, 13470, 31402, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54437, 3, 13470, 13476, 31412, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54452, 3, 13476, 13482, 31422, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54467, 3, 13494, 13500, 31442, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54482, 3, 13500, 13506, 31452, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54497, 3, 13506, 13512, 31462, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54512, 3, 13512, 13518, 31472, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54527, 3, 13518, 13524, 31482, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54542, 3, 13524, 13530, 31492, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54557, 3, 13530, 13536, 31502, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54572, 3, 13536, 13542, 31512, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 54587, 3, 13542, 13548, 31522, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 54602, 0, 3, 31332, 54332, 31562, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 54647, 0, 3, 31342, 54347, 31592, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 54692, 0, 3, 31352, 54362, 31622, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 54737, 0, 3, 31362, 54377, 31652, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 54782, 0, 3, 31372, 54392, 31682, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 54827, 0, 3, 31382, 54407, 31712, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 54872, 0, 3, 31392, 54422, 31742, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 54917, 0, 3, 31402, 54437, 31772, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 54962, 0, 3, 31412, 54452, 31802, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55007, 0, 3, 31432, 54467, 31862, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55052, 0, 3, 31442, 54482, 31892, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55097, 0, 3, 31452, 54497, 31922, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55142, 0, 3, 31462, 54512, 31952, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55187, 0, 3, 31472, 54527, 31982, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55232, 0, 3, 31482, 54542, 32012, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55277, 0, 3, 31492, 54557, 32042, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55322, 0, 3, 31502, 54572, 32072, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 55367, 0, 3, 31512, 54587, 32102, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 55412, 0, 3, 31532, 54602, 13956, 13992,
                                                 32192, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55502, 0, 3, 31562, 54647, 13992, 14028,
                                                 32252, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55592, 0, 3, 31592, 54692, 14028, 14064,
                                                 32312, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55682, 0, 3, 31622, 54737, 14064, 14100,
                                                 32372, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55772, 0, 3, 31652, 54782, 14100, 14136,
                                                 32432, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55862, 0, 3, 31682, 54827, 14136, 14172,
                                                 32492, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 55952, 0, 3, 31712, 54872, 14172, 14208,
                                                 32552, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56042, 0, 3, 31742, 54917, 14208, 14244,
                                                 32612, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56132, 0, 3, 31772, 54962, 14244, 14280,
                                                 32672, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56222, 0, 3, 31832, 55007, 14352, 14388,
                                                 32792, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56312, 0, 3, 31862, 55052, 14388, 14424,
                                                 32852, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56402, 0, 3, 31892, 55097, 14424, 14460,
                                                 32912, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56492, 0, 3, 31922, 55142, 14460, 14496,
                                                 32972, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56582, 0, 3, 31952, 55187, 14496, 14532,
                                                 33032, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56672, 0, 3, 31982, 55232, 14532, 14568,
                                                 33092, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56762, 0, 3, 32012, 55277, 14568, 14604,
                                                 33152, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56852, 0, 3, 32042, 55322, 14604, 14640,
                                                 33212, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 56942, 0, 3, 32072, 55367, 14640, 14676,
                                                 33272, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57032, 0, 3, 32132, 55412, 14748, 14808,
                                                 33332, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57182, 0, 3, 32192, 55502, 14808, 14868,
                                                 33432, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57332, 0, 3, 32252, 55592, 14868, 14928,
                                                 33532, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57482, 0, 3, 32312, 55682, 14928, 14988,
                                                 33632, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57632, 0, 3, 32372, 55772, 14988, 15048,
                                                 33732, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57782, 0, 3, 32432, 55862, 15048, 15108,
                                                 33832, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 57932, 0, 3, 32492, 55952, 15108, 15168,
                                                 33932, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58082, 0, 3, 32552, 56042, 15168, 15228,
                                                 34032, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58232, 0, 3, 32612, 56132, 15228, 15288,
                                                 34132, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58382, 0, 3, 32732, 56222, 15408, 15468,
                                                 34232, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58532, 0, 3, 32792, 56312, 15468, 15528,
                                                 34332, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58682, 0, 3, 32852, 56402, 15528, 15588,
                                                 34432, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58832, 0, 3, 32912, 56492, 15588, 15648,
                                                 34532, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 58982, 0, 3, 32972, 56582, 15648, 15708,
                                                 34632, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 59132, 0, 3, 33032, 56672, 15708, 15768,
                                                 34732, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 59282, 0, 3, 33092, 56762, 15768, 15828,
                                                 34832, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 59432, 0, 3, 33152, 56852, 15828, 15888,
                                                 34932, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 59582, 0, 3, 33212, 56942, 15888, 15948,
                                                 35032, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 59732, 0, 3, 55412, 55502, 33432, 57332,
                                                 16068, 16158, 35282, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 59957, 0, 3, 55502, 55592, 33532, 57482,
                                                 16158, 16248, 35432, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60182, 0, 3, 55592, 55682, 33632, 57632,
                                                 16248, 16338, 35582, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60407, 0, 3, 55682, 55772, 33732, 57782,
                                                 16338, 16428, 35732, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60632, 0, 3, 55772, 55862, 33832, 57932,
                                                 16428, 16518, 35882, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 60857, 0, 3, 55862, 55952, 33932, 58082,
                                                 16518, 16608, 36032, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61082, 0, 3, 55952, 56042, 34032, 58232,
                                                 16608, 16698, 36182, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61307, 0, 3, 56222, 56312, 34332, 58682,
                                                 16878, 16968, 36482, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61532, 0, 3, 56312, 56402, 34432, 58832,
                                                 16968, 17058, 36632, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61757, 0, 3, 56402, 56492, 34532, 58982,
                                                 17058, 17148, 36782, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 61982, 0, 3, 56492, 56582, 34632, 59132,
                                                 17148, 17238, 36932, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 62207, 0, 3, 56582, 56672, 34732, 59282,
                                                 17238, 17328, 37082, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 62432, 0, 3, 56672, 56762, 34832, 59432,
                                                 17328, 17418, 37232, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 62657, 0, 3, 56762, 56852, 34932, 59582,
                                                 17418, 17508, 37382, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 62882, 0, 3, 57032, 57182, 35132, 59732,
                                                 17688, 17814, 37532, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 63197, 0, 3, 57182, 57332, 35282, 59957,
                                                 17814, 17940, 37742, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 63512, 0, 3, 57332, 57482, 35432, 60182,
                                                 17940, 18066, 37952, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 63827, 0, 3, 57482, 57632, 35582, 60407,
                                                 18066, 18192, 38162, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 64142, 0, 3, 57632, 57782, 35732, 60632,
                                                 18192, 18318, 38372, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 64457, 0, 3, 57782, 57932, 35882, 60857,
                                                 18318, 18444, 38582, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 64772, 0, 3, 57932, 58082, 36032, 61082,
                                                 18444, 18570, 38792, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 65087, 0, 3, 58382, 58532, 36332, 61307,
                                                 18822, 18948, 39002, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 65402, 0, 3, 58532, 58682, 36482, 61532,
                                                 18948, 19074, 39212, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 65717, 0, 3, 58682, 58832, 36632, 61757,
                                                 19074, 19200, 39422, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 66032, 0, 3, 58832, 58982, 36782, 61982,
                                                 19200, 19326, 39632, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 66347, 0, 3, 58982, 59132, 36932, 62207,
                                                 19326, 19452, 39842, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 66662, 0, 3, 59132, 59282, 37082, 62432,
                                                 19452, 19578, 40052, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 66977, 0, 3, 59282, 59432, 37232, 62657,
                                                 19578, 19704, 40262, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 67292, 0, 3, 59732, 59957, 37742, 63512,
                                                 19956, 20124, 40752, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 67712, 0, 3, 59957, 60182, 37952, 63827,
                                                 20124, 20292, 41032, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 68132, 0, 3, 60182, 60407, 38162, 64142,
                                                 20292, 20460, 41312, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 68552, 0, 3, 60407, 60632, 38372, 64457,
                                                 20460, 20628, 41592, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 68972, 0, 3, 60632, 60857, 38582, 64772,
                                                 20628, 20796, 41872, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 69392, 0, 3, 61307, 61532, 39212, 65717,
                                                 21132, 21300, 42432, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 69812, 0, 3, 61532, 61757, 39422, 66032,
                                                 21300, 21468, 42712, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 70232, 0, 3, 61757, 61982, 39632, 66347,
                                                 21468, 21636, 42992, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 70652, 0, 3, 61982, 62207, 39842, 66662,
                                                 21636, 21804, 43272, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 71072, 0, 3, 62207, 62432, 40052, 66977,
                                                 21804, 21972, 43552, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 71492, 0, 3, 62882, 63197, 40472, 67292,
                                                 22308, 22524, 43832, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 72032, 0, 3, 63197, 63512, 40752, 67712,
                                                 22524, 22740, 44192, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 72572, 0, 3, 63512, 63827, 41032, 68132,
                                                 22740, 22956, 44552, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 73112, 0, 3, 63827, 64142, 41312, 68552,
                                                 22956, 23172, 44912, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 73652, 0, 3, 64142, 64457, 41592, 68972,
                                                 23172, 23388, 45272, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 74192, 0, 3, 65087, 65402, 42152, 69392,
                                                 23820, 24036, 45632, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 74732, 0, 3, 65402, 65717, 42432, 69812,
                                                 24036, 24252, 45992, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 75272, 0, 3, 65717, 66032, 42712, 70232,
                                                 24252, 24468, 46352, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 75812, 0, 3, 66032, 66347, 42992, 70652,
                                                 24468, 24684, 46712, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 76352, 0, 3, 66347, 66662, 43272, 71072,
                                                 24684, 24900, 47072, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 76892, 0, 3, 67292, 67712, 44192, 72572,
                                                 25332, 25602, 47882, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 77567, 0, 3, 67712, 68132, 44552, 73112,
                                                 25602, 25872, 48332, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 78242, 0, 3, 68132, 68552, 44912, 73652,
                                                 25872, 26142, 48782, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 78917, 0, 3, 69392, 69812, 45992, 75272,
                                                 26682, 26952, 49682, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 79592, 0, 3, 69812, 70232, 46352, 75812,
                                                 26952, 27222, 50132, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 80267, 0, 3, 70232, 70652, 46712, 76352,
                                                 27222, 27492, 50582, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 80942, 0, 3, 71492, 72032, 47432, 76892,
                                                 28032, 28362, 51032, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 81767, 0, 3, 72032, 72572, 47882, 77567,
                                                 28362, 28692, 51582, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 82592, 0, 3, 72572, 73112, 48332, 78242,
                                                 28692, 29022, 52132, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 83417, 0, 3, 74192, 74732, 49232, 78917,
                                                 29682, 30012, 52682, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 84242, 0, 3, 74732, 75272, 49682, 79592,
                                                 30012, 30342, 53232, ncols, alpha, beta, p);

            compute_prim_mg_electron_repulsion_0(buffer, 85067, 0, 3, 75272, 75812, 50132, 80267,
                                                 30342, 30672, 53782, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85892, 3, 31332, 31342, 54347, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85913, 3, 31342, 31352, 54362, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85934, 3, 31352, 31362, 54377, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85955, 3, 31362, 31372, 54392, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85976, 3, 31372, 31382, 54407, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 85997, 3, 31382, 31392, 54422, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86018, 3, 31392, 31402, 54437, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86039, 3, 31402, 31412, 54452, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86060, 3, 31432, 31442, 54482, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86081, 3, 31442, 31452, 54497, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86102, 3, 31452, 31462, 54512, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86123, 3, 31462, 31472, 54527, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86144, 3, 31472, 31482, 54542, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86165, 3, 31482, 31492, 54557, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86186, 3, 31492, 31502, 54572, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 86207, 3, 31502, 31512, 54587, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 86228, 0, 3, 54332, 85892, 54647, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86291, 0, 3, 54347, 85913, 54692, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86354, 0, 3, 54362, 85934, 54737, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86417, 0, 3, 54377, 85955, 54782, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86480, 0, 3, 54392, 85976, 54827, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86543, 0, 3, 54407, 85997, 54872, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86606, 0, 3, 54422, 86018, 54917, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86669, 0, 3, 54437, 86039, 54962, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86732, 0, 3, 54467, 86060, 55052, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86795, 0, 3, 54482, 86081, 55097, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86858, 0, 3, 54497, 86102, 55142, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86921, 0, 3, 54512, 86123, 55187, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 86984, 0, 3, 54527, 86144, 55232, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 87047, 0, 3, 54542, 86165, 55277, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 87110, 0, 3, 54557, 86186, 55322, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 87173, 0, 3, 54572, 86207, 55367, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 87236, 0, 3, 54602, 86228, 32132, 32192,
                                                 55502, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87362, 0, 3, 54647, 86291, 32192, 32252,
                                                 55592, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87488, 0, 3, 54692, 86354, 32252, 32312,
                                                 55682, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87614, 0, 3, 54737, 86417, 32312, 32372,
                                                 55772, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87740, 0, 3, 54782, 86480, 32372, 32432,
                                                 55862, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87866, 0, 3, 54827, 86543, 32432, 32492,
                                                 55952, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 87992, 0, 3, 54872, 86606, 32492, 32552,
                                                 56042, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88118, 0, 3, 54917, 86669, 32552, 32612,
                                                 56132, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88244, 0, 3, 55007, 86732, 32732, 32792,
                                                 56312, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88370, 0, 3, 55052, 86795, 32792, 32852,
                                                 56402, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88496, 0, 3, 55097, 86858, 32852, 32912,
                                                 56492, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88622, 0, 3, 55142, 86921, 32912, 32972,
                                                 56582, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88748, 0, 3, 55187, 86984, 32972, 33032,
                                                 56672, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 88874, 0, 3, 55232, 87047, 33032, 33092,
                                                 56762, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 89000, 0, 3, 55277, 87110, 33092, 33152,
                                                 56852, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 89126, 0, 3, 55322, 87173, 33152, 33212,
                                                 56942, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89252, 0, 3, 55502, 87362, 33332, 33432,
                                                 57332, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89462, 0, 3, 55592, 87488, 33432, 33532,
                                                 57482, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89672, 0, 3, 55682, 87614, 33532, 33632,
                                                 57632, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 89882, 0, 3, 55772, 87740, 33632, 33732,
                                                 57782, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90092, 0, 3, 55862, 87866, 33732, 33832,
                                                 57932, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90302, 0, 3, 55952, 87992, 33832, 33932,
                                                 58082, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90512, 0, 3, 56042, 88118, 33932, 34032,
                                                 58232, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90722, 0, 3, 56312, 88370, 34232, 34332,
                                                 58682, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 90932, 0, 3, 56402, 88496, 34332, 34432,
                                                 58832, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91142, 0, 3, 56492, 88622, 34432, 34532,
                                                 58982, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91352, 0, 3, 56582, 88748, 34532, 34632,
                                                 59132, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91562, 0, 3, 56672, 88874, 34632, 34732,
                                                 59282, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91772, 0, 3, 56762, 89000, 34732, 34832,
                                                 59432, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 91982, 0, 3, 56852, 89126, 34832, 34932,
                                                 59582, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 92192, 0, 3, 87236, 87362, 57332, 89462,
                                                 35132, 35282, 59957, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 92507, 0, 3, 87362, 87488, 57482, 89672,
                                                 35282, 35432, 60182, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 92822, 0, 3, 87488, 87614, 57632, 89882,
                                                 35432, 35582, 60407, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 93137, 0, 3, 87614, 87740, 57782, 90092,
                                                 35582, 35732, 60632, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 93452, 0, 3, 87740, 87866, 57932, 90302,
                                                 35732, 35882, 60857, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 93767, 0, 3, 87866, 87992, 58082, 90512,
                                                 35882, 36032, 61082, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 94082, 0, 3, 88244, 88370, 58682, 90932,
                                                 36332, 36482, 61532, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 94397, 0, 3, 88370, 88496, 58832, 91142,
                                                 36482, 36632, 61757, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 94712, 0, 3, 88496, 88622, 58982, 91352,
                                                 36632, 36782, 61982, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 95027, 0, 3, 88622, 88748, 59132, 91562,
                                                 36782, 36932, 62207, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 95342, 0, 3, 88748, 88874, 59282, 91772,
                                                 36932, 37082, 62432, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 95657, 0, 3, 88874, 89000, 59432, 91982,
                                                 37082, 37232, 62657, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 95972, 0, 3, 89252, 89462, 59957, 92507,
                                                 37532, 37742, 63512, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 96413, 0, 3, 89462, 89672, 60182, 92822,
                                                 37742, 37952, 63827, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 96854, 0, 3, 89672, 89882, 60407, 93137,
                                                 37952, 38162, 64142, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 97295, 0, 3, 89882, 90092, 60632, 93452,
                                                 38162, 38372, 64457, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 97736, 0, 3, 90092, 90302, 60857, 93767,
                                                 38372, 38582, 64772, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 98177, 0, 3, 90722, 90932, 61532, 94397,
                                                 39002, 39212, 65717, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 98618, 0, 3, 90932, 91142, 61757, 94712,
                                                 39212, 39422, 66032, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 99059, 0, 3, 91142, 91352, 61982, 95027,
                                                 39422, 39632, 66347, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 99500, 0, 3, 91352, 91562, 62207, 95342,
                                                 39632, 39842, 66662, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 99941, 0, 3, 91562, 91772, 62432, 95657,
                                                 39842, 40052, 66977, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 100382, 0, 3, 92192, 92507, 63512,
                                                 96413, 40472, 40752, 67712, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 100970, 0, 3, 92507, 92822, 63827,
                                                 96854, 40752, 41032, 68132, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 101558, 0, 3, 92822, 93137, 64142,
                                                 97295, 41032, 41312, 68552, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 102146, 0, 3, 93137, 93452, 64457,
                                                 97736, 41312, 41592, 68972, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 102734, 0, 3, 94082, 94397, 65717,
                                                 98618, 42152, 42432, 69812, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 103322, 0, 3, 94397, 94712, 66032,
                                                 99059, 42432, 42712, 70232, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 103910, 0, 3, 94712, 95027, 66347,
                                                 99500, 42712, 42992, 70652, ncols, alpha, beta,
                                                 p);

            compute_prim_ih_electron_repulsion_0(buffer, 104498, 0, 3, 95027, 95342, 66662,
                                                 99941, 42992, 43272, 71072, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 105086, 0, 3, 95972, 96413, 67712,
                                                 100970, 43832, 44192, 72572, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 105842, 0, 3, 96413, 96854, 68132,
                                                 101558, 44192, 44552, 73112, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 106598, 0, 3, 96854, 97295, 68552,
                                                 102146, 44552, 44912, 73652, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 107354, 0, 3, 98177, 98618, 69812,
                                                 103322, 45632, 45992, 75272, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 108110, 0, 3, 98618, 99059, 70232,
                                                 103910, 45992, 46352, 75812, ncols, alpha, beta,
                                                 p);

            compute_prim_kh_electron_repulsion_0(buffer, 108866, 0, 3, 99059, 99500, 70652,
                                                 104498, 46352, 46712, 76352, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 109622, 0, 3, 100382, 100970, 72572,
                                                 105842, 47432, 47882, 77567, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 110567, 0, 3, 100970, 101558, 73112,
                                                 106598, 47882, 48332, 78242, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 111512, 0, 3, 102734, 103322, 75272,
                                                 108110, 49232, 49682, 79592, ncols, alpha, beta,
                                                 p);

            compute_prim_lh_electron_repulsion_0(buffer, 112457, 0, 3, 103322, 103910, 75812,
                                                 108866, 49682, 50132, 80267, ncols, alpha, beta,
                                                 p);

            compute_prim_mh_electron_repulsion_0(buffer, 113402, 0, 3, 105086, 105842, 77567,
                                                 110567, 51032, 51582, 82592, ncols, alpha, beta,
                                                 p);

            compute_prim_mh_electron_repulsion_0(buffer, 114557, 0, 3, 107354, 108110, 79592,
                                                 112457, 52682, 53232, 85067, ncols, alpha, beta,
                                                 p);

            compute_prim_si_electron_repulsion_0(buffer, 115712, 3, 54332, 54347, 85913, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115740, 3, 54347, 54362, 85934, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115768, 3, 54362, 54377, 85955, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115796, 3, 54377, 54392, 85976, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115824, 3, 54392, 54407, 85997, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115852, 3, 54407, 54422, 86018, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115880, 3, 54422, 54437, 86039, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115908, 3, 54467, 54482, 86081, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115936, 3, 54482, 54497, 86102, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115964, 3, 54497, 54512, 86123, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 115992, 3, 54512, 54527, 86144, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 116020, 3, 54527, 54542, 86165, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 116048, 3, 54542, 54557, 86186, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 116076, 3, 54557, 54572, 86207, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116104, 0, 3, 85892, 115712, 86291,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116188, 0, 3, 85913, 115740, 86354,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116272, 0, 3, 85934, 115768, 86417,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116356, 0, 3, 85955, 115796, 86480,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116440, 0, 3, 85976, 115824, 86543,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116524, 0, 3, 85997, 115852, 86606,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116608, 0, 3, 86018, 115880, 86669,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116692, 0, 3, 86060, 115908, 86795,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116776, 0, 3, 86081, 115936, 86858,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116860, 0, 3, 86102, 115964, 86921,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 116944, 0, 3, 86123, 115992, 86984,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 117028, 0, 3, 86144, 116020, 87047,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 117112, 0, 3, 86165, 116048, 87110,
                                                 ncols, p);

            compute_prim_pi_electron_repulsion_0(buffer, 117196, 0, 3, 86186, 116076, 87173,
                                                 ncols, p);

            compute_prim_di_electron_repulsion_0(buffer, 117280, 0, 3, 86228, 116104, 55412,
                                                 55502, 87362, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 117448, 0, 3, 86291, 116188, 55502,
                                                 55592, 87488, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 117616, 0, 3, 86354, 116272, 55592,
                                                 55682, 87614, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 117784, 0, 3, 86417, 116356, 55682,
                                                 55772, 87740, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 117952, 0, 3, 86480, 116440, 55772,
                                                 55862, 87866, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 118120, 0, 3, 86543, 116524, 55862,
                                                 55952, 87992, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 118288, 0, 3, 86606, 116608, 55952,
                                                 56042, 88118, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 118456, 0, 3, 86732, 116692, 56222,
                                                 56312, 88370, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 118624, 0, 3, 86795, 116776, 56312,
                                                 56402, 88496, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 118792, 0, 3, 86858, 116860, 56402,
                                                 56492, 88622, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 118960, 0, 3, 86921, 116944, 56492,
                                                 56582, 88748, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 119128, 0, 3, 86984, 117028, 56582,
                                                 56672, 88874, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 119296, 0, 3, 87047, 117112, 56672,
                                                 56762, 89000, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 119464, 0, 3, 87110, 117196, 56762,
                                                 56852, 89126, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 119632, 0, 3, 87236, 117280, 57032,
                                                 57182, 89252, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 119912, 0, 3, 87362, 117448, 57182,
                                                 57332, 89462, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 120192, 0, 3, 87488, 117616, 57332,
                                                 57482, 89672, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 120472, 0, 3, 87614, 117784, 57482,
                                                 57632, 89882, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 120752, 0, 3, 87740, 117952, 57632,
                                                 57782, 90092, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 121032, 0, 3, 87866, 118120, 57782,
                                                 57932, 90302, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 121312, 0, 3, 87992, 118288, 57932,
                                                 58082, 90512, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 121592, 0, 3, 88244, 118456, 58382,
                                                 58532, 90722, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 121872, 0, 3, 88370, 118624, 58532,
                                                 58682, 90932, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 122152, 0, 3, 88496, 118792, 58682,
                                                 58832, 91142, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 122432, 0, 3, 88622, 118960, 58832,
                                                 58982, 91352, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 122712, 0, 3, 88748, 119128, 58982,
                                                 59132, 91562, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 122992, 0, 3, 88874, 119296, 59132,
                                                 59282, 91772, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 123272, 0, 3, 89000, 119464, 59282,
                                                 59432, 91982, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 123552, 0, 3, 117280, 117448, 89462,
                                                 120192, 59732, 59957, 92507, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 123972, 0, 3, 117448, 117616, 89672,
                                                 120472, 59957, 60182, 92822, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 124392, 0, 3, 117616, 117784, 89882,
                                                 120752, 60182, 60407, 93137, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 124812, 0, 3, 117784, 117952, 90092,
                                                 121032, 60407, 60632, 93452, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 125232, 0, 3, 117952, 118120, 90302,
                                                 121312, 60632, 60857, 93767, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 125652, 0, 3, 118456, 118624, 90932,
                                                 122152, 61307, 61532, 94397, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 126072, 0, 3, 118624, 118792, 91142,
                                                 122432, 61532, 61757, 94712, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 126492, 0, 3, 118792, 118960, 91352,
                                                 122712, 61757, 61982, 95027, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 126912, 0, 3, 118960, 119128, 91562,
                                                 122992, 61982, 62207, 95342, ncols, alpha, beta,
                                                 p);

            compute_prim_gi_electron_repulsion_0(buffer, 127332, 0, 3, 119128, 119296, 91772,
                                                 123272, 62207, 62432, 95657, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 127752, 0, 3, 119632, 119912, 92192,
                                                 123552, 62882, 63197, 95972, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 128340, 0, 3, 119912, 120192, 92507,
                                                 123972, 63197, 63512, 96413, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 128928, 0, 3, 120192, 120472, 92822,
                                                 124392, 63512, 63827, 96854, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 129516, 0, 3, 120472, 120752, 93137,
                                                 124812, 63827, 64142, 97295, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 130104, 0, 3, 120752, 121032, 93452,
                                                 125232, 64142, 64457, 97736, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 130692, 0, 3, 121592, 121872, 94082,
                                                 125652, 65087, 65402, 98177, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 131280, 0, 3, 121872, 122152, 94397,
                                                 126072, 65402, 65717, 98618, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 131868, 0, 3, 122152, 122432, 94712,
                                                 126492, 65717, 66032, 99059, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 132456, 0, 3, 122432, 122712, 95027,
                                                 126912, 66032, 66347, 99500, ncols, alpha, beta,
                                                 p);

            compute_prim_hi_electron_repulsion_0(buffer, 133044, 0, 3, 122712, 122992, 95342,
                                                 127332, 66347, 66662, 99941, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 133632, 0, 3, 123552, 123972, 96413,
                                                 128928, 67292, 67712, 100970, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 134416, 0, 3, 123972, 124392, 96854,
                                                 129516, 67712, 68132, 101558, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 135200, 0, 3, 124392, 124812, 97295,
                                                 130104, 68132, 68552, 102146, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 135984, 0, 3, 125652, 126072, 98618,
                                                 131868, 69392, 69812, 103322, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 136768, 0, 3, 126072, 126492, 99059,
                                                 132456, 69812, 70232, 103910, ncols, alpha,
                                                 beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 137552, 0, 3, 126492, 126912, 99500,
                                                 133044, 70232, 70652, 104498, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 138336, 0, 3, 127752, 128340, 100382,
                                                 133632, 71492, 72032, 105086, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 139344, 0, 3, 128340, 128928, 100970,
                                                 134416, 72032, 72572, 105842, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 140352, 0, 3, 128928, 129516, 101558,
                                                 135200, 72572, 73112, 106598, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 141360, 0, 3, 130692, 131280, 102734,
                                                 135984, 74192, 74732, 107354, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 142368, 0, 3, 131280, 131868, 103322,
                                                 136768, 74732, 75272, 108110, ncols, alpha,
                                                 beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 143376, 0, 3, 131868, 132456, 103910,
                                                 137552, 75272, 75812, 108866, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 144384, 0, 3, 133632, 134416, 105842,
                                                 140352, 76892, 77567, 110567, ncols, alpha,
                                                 beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 145644, 0, 3, 135984, 136768, 108110,
                                                 143376, 78917, 79592, 112457, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 146904, 0, 3, 138336, 139344, 109622,
                                                 144384, 80942, 81767, 113402, ncols, alpha,
                                                 beta, p);

            compute_prim_mi_electron_repulsion_0(buffer, 148444, 0, 3, 141360, 142368, 111512,
                                                 145644, 83417, 84242, 114557, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_l_x(buffer, 149984, 141360, 148444, 1, 28, ncols, alpha);

            simdgeo::geom_l_y(buffer, 151244, 141360, 148444, 1, 28, ncols, alpha);

            simdgeo::geom_l_z(buffer, 152504, 141360, 148444, 1, 28, ncols, alpha);

            simdgeo::geom_l_x(buffer, 153764, 138336, 146904, 1, 28, ncols, alpha);

            simdgeo::geom_l_y(buffer, 155024, 138336, 146904, 1, 28, ncols, alpha);

            simdgeo::geom_l_z(buffer, 156284, 138336, 146904, 1, 28, ncols, alpha);

            simdfunc::contract_primitives(buffer, 157544, 153764, 3780, ncols);

            simdfunc::contract_primitives(buffer, 161324, 149984, 3780, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 165104, 161324, 45, 1, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 165104, 13, nmax);

    simdtrf::transform_i_inner(buffer, 165104, 162584, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 221 * nvalues, nvalues, buffer, 165104, 13, nmax);

    simdtrf::transform_i_inner(buffer, 165104, 163844, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 442 * nvalues, nvalues, buffer, 165104, 13, nmax);

    simdtrf::transform_i_inner(buffer, 165104, 157544, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 663 * nvalues, nvalues, buffer, 165104, 13, nmax);

    simdtrf::transform_i_inner(buffer, 165104, 158804, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 884 * nvalues, nvalues, buffer, 165104, 13, nmax);

    simdtrf::transform_i_inner(buffer, 165104, 160064, 45, 1, nmax);

    simdtrf::transform_l_outer(values + 1105 * nvalues, nvalues, buffer, 165104, 13, nmax);
}

}  // namespace simdt2ceri
