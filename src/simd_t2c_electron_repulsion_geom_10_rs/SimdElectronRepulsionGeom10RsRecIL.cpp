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


#include "SimdElectronRepulsionGeom10RsRecIL.hpp"

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
#include "SimdGeometryI1.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_il_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_il_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 175200, 167164, 7560, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 2488, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2491, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2494, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2497, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2500, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2503, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2506, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2509, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2512, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2515, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2518, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2521, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2524, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2527, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2530, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2533, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2536, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2539, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2542, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2545, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2548, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2551, 3, 35, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2554, 3, 36, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2557, 3, 37, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2560, 3, 9, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2569, 3, 10, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2578, 3, 11, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2587, 3, 12, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2596, 3, 13, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2605, 3, 14, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2614, 3, 15, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2623, 3, 16, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2632, 3, 17, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2641, 3, 18, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2650, 3, 19, 77, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2659, 3, 20, 80, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2668, 3, 25, 92, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2677, 3, 26, 95, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2686, 3, 27, 98, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2695, 3, 28, 101, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2704, 3, 29, 104, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2713, 3, 30, 107, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2722, 3, 31, 110, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2731, 3, 32, 113, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2740, 3, 33, 116, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2749, 3, 34, 119, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2758, 3, 35, 122, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2767, 3, 36, 125, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2776, 0, 3, 44, 2560, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2794, 0, 3, 47, 2569, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2812, 0, 3, 50, 2578, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2830, 0, 3, 53, 2587, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2848, 0, 3, 56, 2596, 158, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2866, 0, 3, 59, 2605, 164, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2884, 0, 3, 62, 2614, 170, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2902, 0, 3, 65, 2623, 176, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2920, 0, 3, 68, 2632, 182, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2938, 0, 3, 71, 2641, 188, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2956, 0, 3, 74, 2650, 194, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2974, 0, 3, 77, 2659, 200, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2992, 0, 3, 89, 2668, 212, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3010, 0, 3, 92, 2677, 218, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3028, 0, 3, 95, 2686, 224, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3046, 0, 3, 98, 2695, 230, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3064, 0, 3, 101, 2704, 236, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3082, 0, 3, 104, 2713, 242, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3100, 0, 3, 107, 2722, 248, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3118, 0, 3, 110, 2731, 254, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3136, 0, 3, 113, 2740, 260, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3154, 0, 3, 116, 2749, 266, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3172, 0, 3, 119, 2758, 272, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3190, 0, 3, 122, 2767, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3208, 0, 3, 134, 2794, 304, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3238, 0, 3, 140, 2812, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3268, 0, 3, 146, 2830, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3298, 0, 3, 152, 2848, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3328, 0, 3, 158, 2866, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3358, 0, 3, 164, 2884, 354, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3388, 0, 3, 170, 2902, 364, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3418, 0, 3, 176, 2920, 374, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3448, 0, 3, 182, 2938, 384, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3478, 0, 3, 188, 2956, 394, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3508, 0, 3, 194, 2974, 404, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3538, 0, 3, 212, 3010, 434, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3568, 0, 3, 218, 3028, 444, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3598, 0, 3, 224, 3046, 454, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3628, 0, 3, 230, 3064, 464, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3658, 0, 3, 236, 3082, 474, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3688, 0, 3, 242, 3100, 484, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3718, 0, 3, 248, 3118, 494, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3748, 0, 3, 254, 3136, 504, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3778, 0, 3, 260, 3154, 514, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3808, 0, 3, 266, 3172, 524, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3838, 0, 3, 272, 3190, 534, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3868, 0, 3, 304, 3238, 559, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3913, 0, 3, 314, 3268, 574, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3958, 0, 3, 324, 3298, 589, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4003, 0, 3, 334, 3328, 604, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4048, 0, 3, 344, 3358, 619, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4093, 0, 3, 354, 3388, 634, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4138, 0, 3, 364, 3418, 649, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4183, 0, 3, 374, 3448, 664, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4228, 0, 3, 384, 3478, 679, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4273, 0, 3, 394, 3508, 694, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4318, 0, 3, 434, 3568, 724, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4363, 0, 3, 444, 3598, 739, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4408, 0, 3, 454, 3628, 754, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4453, 0, 3, 464, 3658, 769, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4498, 0, 3, 474, 3688, 784, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4543, 0, 3, 484, 3718, 799, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4588, 0, 3, 494, 3748, 814, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4633, 0, 3, 504, 3778, 829, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4678, 0, 3, 514, 3808, 844, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4723, 0, 3, 524, 3838, 859, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4768, 0, 3, 559, 3913, 916, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4831, 0, 3, 574, 3958, 937, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4894, 0, 3, 589, 4003, 958, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4957, 0, 3, 604, 4048, 979, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5020, 0, 3, 619, 4093, 1000, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5083, 0, 3, 634, 4138, 1021, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5146, 0, 3, 649, 4183, 1042, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5209, 0, 3, 664, 4228, 1063, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5272, 0, 3, 679, 4273, 1084, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5335, 0, 3, 724, 4363, 1147, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5398, 0, 3, 739, 4408, 1168, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5461, 0, 3, 754, 4453, 1189, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5524, 0, 3, 769, 4498, 1210, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5587, 0, 3, 784, 4543, 1231, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5650, 0, 3, 799, 4588, 1252, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5713, 0, 3, 814, 4633, 1273, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5776, 0, 3, 829, 4678, 1294, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5839, 0, 3, 844, 4723, 1315, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5902, 0, 3, 916, 4831, 1364, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5986, 0, 3, 937, 4894, 1392, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6070, 0, 3, 958, 4957, 1420, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6154, 0, 3, 979, 5020, 1448, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6238, 0, 3, 1000, 5083, 1476, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6322, 0, 3, 1021, 5146, 1504, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6406, 0, 3, 1042, 5209, 1532, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6490, 0, 3, 1063, 5272, 1560, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6574, 0, 3, 1147, 5398, 1616, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6658, 0, 3, 1168, 5461, 1644, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6742, 0, 3, 1189, 5524, 1672, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6826, 0, 3, 1210, 5587, 1700, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6910, 0, 3, 1231, 5650, 1728, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6994, 0, 3, 1252, 5713, 1756, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7078, 0, 3, 1273, 5776, 1784, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7162, 0, 3, 1294, 5839, 1812, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7246, 0, 3, 1364, 5986, 1912, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7354, 0, 3, 1392, 6070, 1948, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7462, 0, 3, 1420, 6154, 1984, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7570, 0, 3, 1448, 6238, 2020, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7678, 0, 3, 1476, 6322, 2056, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7786, 0, 3, 1504, 6406, 2092, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7894, 0, 3, 1532, 6490, 2128, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8002, 0, 3, 1616, 6658, 2236, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8110, 0, 3, 1644, 6742, 2272, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8218, 0, 3, 1672, 6826, 2308, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8326, 0, 3, 1700, 6910, 2344, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8434, 0, 3, 1728, 6994, 2380, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8542, 0, 3, 1756, 7078, 2416, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8650, 0, 3, 1784, 7162, 2452, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 8758, 3, 9, 10, 2491, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8764, 3, 10, 11, 2494, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8770, 3, 11, 12, 2497, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8776, 3, 12, 13, 2500, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8782, 3, 13, 14, 2503, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8788, 3, 14, 15, 2506, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8794, 3, 15, 16, 2509, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8800, 3, 16, 17, 2512, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8806, 3, 17, 18, 2515, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8812, 3, 18, 19, 2518, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8818, 3, 19, 20, 2521, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8824, 3, 25, 26, 2527, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8830, 3, 26, 27, 2530, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8836, 3, 27, 28, 2533, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8842, 3, 28, 29, 2536, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8848, 3, 29, 30, 2539, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8854, 3, 30, 31, 2542, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8860, 3, 31, 32, 2545, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8866, 3, 32, 33, 2548, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8872, 3, 33, 34, 2551, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8878, 3, 34, 35, 2554, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8884, 3, 35, 36, 2557, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 8890, 0, 3, 2488, 8758, 2569, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8908, 0, 3, 2491, 8764, 2578, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8926, 0, 3, 2494, 8770, 2587, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8944, 0, 3, 2497, 8776, 2596, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8962, 0, 3, 2500, 8782, 2605, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8980, 0, 3, 2503, 8788, 2614, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8998, 0, 3, 2506, 8794, 2623, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9016, 0, 3, 2509, 8800, 2632, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9034, 0, 3, 2512, 8806, 2641, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9052, 0, 3, 2515, 8812, 2650, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9070, 0, 3, 2518, 8818, 2659, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9088, 0, 3, 2524, 8824, 2677, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9106, 0, 3, 2527, 8830, 2686, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9124, 0, 3, 2530, 8836, 2695, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9142, 0, 3, 2533, 8842, 2704, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9160, 0, 3, 2536, 8848, 2713, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9178, 0, 3, 2539, 8854, 2722, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9196, 0, 3, 2542, 8860, 2731, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9214, 0, 3, 2545, 8866, 2740, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9232, 0, 3, 2548, 8872, 2749, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9250, 0, 3, 2551, 8878, 2758, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9268, 0, 3, 2554, 8884, 2767, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 9286, 0, 3, 2560, 8890, 128, 134, 2794,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9322, 0, 3, 2569, 8908, 134, 140, 2812,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9358, 0, 3, 2578, 8926, 140, 146, 2830,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9394, 0, 3, 2587, 8944, 146, 152, 2848,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9430, 0, 3, 2596, 8962, 152, 158, 2866,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9466, 0, 3, 2605, 8980, 158, 164, 2884,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9502, 0, 3, 2614, 8998, 164, 170, 2902,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9538, 0, 3, 2623, 9016, 170, 176, 2920,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9574, 0, 3, 2632, 9034, 176, 182, 2938,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9610, 0, 3, 2641, 9052, 182, 188, 2956,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9646, 0, 3, 2650, 9070, 188, 194, 2974,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9682, 0, 3, 2668, 9088, 206, 212, 3010,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9718, 0, 3, 2677, 9106, 212, 218, 3028,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9754, 0, 3, 2686, 9124, 218, 224, 3046,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9790, 0, 3, 2695, 9142, 224, 230, 3064,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9826, 0, 3, 2704, 9160, 230, 236, 3082,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9862, 0, 3, 2713, 9178, 236, 242, 3100,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9898, 0, 3, 2722, 9196, 242, 248, 3118,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9934, 0, 3, 2731, 9214, 248, 254, 3136,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9970, 0, 3, 2740, 9232, 254, 260, 3154,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10006, 0, 3, 2749, 9250, 260, 266, 3172,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10042, 0, 3, 2758, 9268, 266, 272, 3190,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10078, 0, 3, 2776, 9286, 284, 294, 3208,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10138, 0, 3, 2794, 9322, 294, 304, 3238,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10198, 0, 3, 2812, 9358, 304, 314, 3268,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10258, 0, 3, 2830, 9394, 314, 324, 3298,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10318, 0, 3, 2848, 9430, 324, 334, 3328,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10378, 0, 3, 2866, 9466, 334, 344, 3358,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10438, 0, 3, 2884, 9502, 344, 354, 3388,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10498, 0, 3, 2902, 9538, 354, 364, 3418,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10558, 0, 3, 2920, 9574, 364, 374, 3448,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10618, 0, 3, 2938, 9610, 374, 384, 3478,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10678, 0, 3, 2956, 9646, 384, 394, 3508,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10738, 0, 3, 2992, 9682, 414, 424, 3538,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10798, 0, 3, 3010, 9718, 424, 434, 3568,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10858, 0, 3, 3028, 9754, 434, 444, 3598,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10918, 0, 3, 3046, 9790, 444, 454, 3628,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10978, 0, 3, 3064, 9826, 454, 464, 3658,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11038, 0, 3, 3082, 9862, 464, 474, 3688,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11098, 0, 3, 3100, 9898, 474, 484, 3718,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11158, 0, 3, 3118, 9934, 484, 494, 3748,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11218, 0, 3, 3136, 9970, 494, 504, 3778,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11278, 0, 3, 3154, 10006, 504, 514,
                                                 3808, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11338, 0, 3, 3172, 10042, 514, 524,
                                                 3838, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11398, 0, 3, 9286, 9322, 3238, 10198,
                                                 544, 559, 3913, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11488, 0, 3, 9322, 9358, 3268, 10258,
                                                 559, 574, 3958, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11578, 0, 3, 9358, 9394, 3298, 10318,
                                                 574, 589, 4003, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11668, 0, 3, 9394, 9430, 3328, 10378,
                                                 589, 604, 4048, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11758, 0, 3, 9430, 9466, 3358, 10438,
                                                 604, 619, 4093, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11848, 0, 3, 9466, 9502, 3388, 10498,
                                                 619, 634, 4138, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11938, 0, 3, 9502, 9538, 3418, 10558,
                                                 634, 649, 4183, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12028, 0, 3, 9538, 9574, 3448, 10618,
                                                 649, 664, 4228, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12118, 0, 3, 9574, 9610, 3478, 10678,
                                                 664, 679, 4273, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12208, 0, 3, 9682, 9718, 3568, 10858,
                                                 709, 724, 4363, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12298, 0, 3, 9718, 9754, 3598, 10918,
                                                 724, 739, 4408, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12388, 0, 3, 9754, 9790, 3628, 10978,
                                                 739, 754, 4453, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12478, 0, 3, 9790, 9826, 3658, 11038,
                                                 754, 769, 4498, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12568, 0, 3, 9826, 9862, 3688, 11098,
                                                 769, 784, 4543, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12658, 0, 3, 9862, 9898, 3718, 11158,
                                                 784, 799, 4588, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12748, 0, 3, 9898, 9934, 3748, 11218,
                                                 799, 814, 4633, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12838, 0, 3, 9934, 9970, 3778, 11278,
                                                 814, 829, 4678, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12928, 0, 3, 9970, 10006, 3808, 11338,
                                                 829, 844, 4723, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13018, 0, 3, 10078, 10138, 3868, 11398,
                                                 874, 895, 4768, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13144, 0, 3, 10138, 10198, 3913, 11488,
                                                 895, 916, 4831, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13270, 0, 3, 10198, 10258, 3958, 11578,
                                                 916, 937, 4894, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13396, 0, 3, 10258, 10318, 4003, 11668,
                                                 937, 958, 4957, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13522, 0, 3, 10318, 10378, 4048, 11758,
                                                 958, 979, 5020, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13648, 0, 3, 10378, 10438, 4093, 11848,
                                                 979, 1000, 5083, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13774, 0, 3, 10438, 10498, 4138, 11938,
                                                 1000, 1021, 5146, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13900, 0, 3, 10498, 10558, 4183, 12028,
                                                 1021, 1042, 5209, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14026, 0, 3, 10558, 10618, 4228, 12118,
                                                 1042, 1063, 5272, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14152, 0, 3, 10738, 10798, 4318, 12208,
                                                 1105, 1126, 5335, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14278, 0, 3, 10798, 10858, 4363, 12298,
                                                 1126, 1147, 5398, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14404, 0, 3, 10858, 10918, 4408, 12388,
                                                 1147, 1168, 5461, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14530, 0, 3, 10918, 10978, 4453, 12478,
                                                 1168, 1189, 5524, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14656, 0, 3, 10978, 11038, 4498, 12568,
                                                 1189, 1210, 5587, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14782, 0, 3, 11038, 11098, 4543, 12658,
                                                 1210, 1231, 5650, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14908, 0, 3, 11098, 11158, 4588, 12748,
                                                 1231, 1252, 5713, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15034, 0, 3, 11158, 11218, 4633, 12838,
                                                 1252, 1273, 5776, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15160, 0, 3, 11218, 11278, 4678, 12928,
                                                 1273, 1294, 5839, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15286, 0, 3, 11398, 11488, 4831, 13270,
                                                 1336, 1364, 5986, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15454, 0, 3, 11488, 11578, 4894, 13396,
                                                 1364, 1392, 6070, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15622, 0, 3, 11578, 11668, 4957, 13522,
                                                 1392, 1420, 6154, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15790, 0, 3, 11668, 11758, 5020, 13648,
                                                 1420, 1448, 6238, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15958, 0, 3, 11758, 11848, 5083, 13774,
                                                 1448, 1476, 6322, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16126, 0, 3, 11848, 11938, 5146, 13900,
                                                 1476, 1504, 6406, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16294, 0, 3, 11938, 12028, 5209, 14026,
                                                 1504, 1532, 6490, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16462, 0, 3, 12208, 12298, 5398, 14404,
                                                 1588, 1616, 6658, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16630, 0, 3, 12298, 12388, 5461, 14530,
                                                 1616, 1644, 6742, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16798, 0, 3, 12388, 12478, 5524, 14656,
                                                 1644, 1672, 6826, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16966, 0, 3, 12478, 12568, 5587, 14782,
                                                 1672, 1700, 6910, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17134, 0, 3, 12568, 12658, 5650, 14908,
                                                 1700, 1728, 6994, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17302, 0, 3, 12658, 12748, 5713, 15034,
                                                 1728, 1756, 7078, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17470, 0, 3, 12748, 12838, 5776, 15160,
                                                 1756, 1784, 7162, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17638, 0, 3, 13018, 13144, 5902, 15286,
                                                 1840, 1876, 7246, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17854, 0, 3, 13144, 13270, 5986, 15454,
                                                 1876, 1912, 7354, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18070, 0, 3, 13270, 13396, 6070, 15622,
                                                 1912, 1948, 7462, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18286, 0, 3, 13396, 13522, 6154, 15790,
                                                 1948, 1984, 7570, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18502, 0, 3, 13522, 13648, 6238, 15958,
                                                 1984, 2020, 7678, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18718, 0, 3, 13648, 13774, 6322, 16126,
                                                 2020, 2056, 7786, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18934, 0, 3, 13774, 13900, 6406, 16294,
                                                 2056, 2092, 7894, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19150, 0, 3, 14152, 14278, 6574, 16462,
                                                 2164, 2200, 8002, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19366, 0, 3, 14278, 14404, 6658, 16630,
                                                 2200, 2236, 8110, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19582, 0, 3, 14404, 14530, 6742, 16798,
                                                 2236, 2272, 8218, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19798, 0, 3, 14530, 14656, 6826, 16966,
                                                 2272, 2308, 8326, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20014, 0, 3, 14656, 14782, 6910, 17134,
                                                 2308, 2344, 8434, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20230, 0, 3, 14782, 14908, 6994, 17302,
                                                 2344, 2380, 8542, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 20446, 0, 3, 14908, 15034, 7078, 17470,
                                                 2380, 2416, 8650, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20662, 3, 2488, 2491, 8764, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20672, 3, 2491, 2494, 8770, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20682, 3, 2494, 2497, 8776, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20692, 3, 2497, 2500, 8782, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20702, 3, 2500, 2503, 8788, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20712, 3, 2503, 2506, 8794, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20722, 3, 2506, 2509, 8800, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20732, 3, 2509, 2512, 8806, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20742, 3, 2512, 2515, 8812, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20752, 3, 2515, 2518, 8818, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20762, 3, 2524, 2527, 8830, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20772, 3, 2527, 2530, 8836, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20782, 3, 2530, 2533, 8842, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20792, 3, 2533, 2536, 8848, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20802, 3, 2536, 2539, 8854, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20812, 3, 2539, 2542, 8860, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20822, 3, 2542, 2545, 8866, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20832, 3, 2545, 2548, 8872, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20842, 3, 2548, 2551, 8878, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 20852, 3, 2551, 2554, 8884, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 20862, 0, 3, 8758, 20662, 8908, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20892, 0, 3, 8764, 20672, 8926, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20922, 0, 3, 8770, 20682, 8944, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20952, 0, 3, 8776, 20692, 8962, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 20982, 0, 3, 8782, 20702, 8980, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21012, 0, 3, 8788, 20712, 8998, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21042, 0, 3, 8794, 20722, 9016, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21072, 0, 3, 8800, 20732, 9034, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21102, 0, 3, 8806, 20742, 9052, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21132, 0, 3, 8812, 20752, 9070, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21162, 0, 3, 8824, 20762, 9106, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21192, 0, 3, 8830, 20772, 9124, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21222, 0, 3, 8836, 20782, 9142, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21252, 0, 3, 8842, 20792, 9160, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21282, 0, 3, 8848, 20802, 9178, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21312, 0, 3, 8854, 20812, 9196, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21342, 0, 3, 8860, 20822, 9214, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21372, 0, 3, 8866, 20832, 9232, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21402, 0, 3, 8872, 20842, 9250, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 21432, 0, 3, 8878, 20852, 9268, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 21462, 0, 3, 8890, 20862, 2776, 2794,
                                                 9322, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21522, 0, 3, 8908, 20892, 2794, 2812,
                                                 9358, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21582, 0, 3, 8926, 20922, 2812, 2830,
                                                 9394, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21642, 0, 3, 8944, 20952, 2830, 2848,
                                                 9430, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21702, 0, 3, 8962, 20982, 2848, 2866,
                                                 9466, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21762, 0, 3, 8980, 21012, 2866, 2884,
                                                 9502, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21822, 0, 3, 8998, 21042, 2884, 2902,
                                                 9538, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21882, 0, 3, 9016, 21072, 2902, 2920,
                                                 9574, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 21942, 0, 3, 9034, 21102, 2920, 2938,
                                                 9610, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22002, 0, 3, 9052, 21132, 2938, 2956,
                                                 9646, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22062, 0, 3, 9088, 21162, 2992, 3010,
                                                 9718, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22122, 0, 3, 9106, 21192, 3010, 3028,
                                                 9754, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22182, 0, 3, 9124, 21222, 3028, 3046,
                                                 9790, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22242, 0, 3, 9142, 21252, 3046, 3064,
                                                 9826, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22302, 0, 3, 9160, 21282, 3064, 3082,
                                                 9862, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22362, 0, 3, 9178, 21312, 3082, 3100,
                                                 9898, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22422, 0, 3, 9196, 21342, 3100, 3118,
                                                 9934, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22482, 0, 3, 9214, 21372, 3118, 3136,
                                                 9970, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22542, 0, 3, 9232, 21402, 3136, 3154,
                                                 10006, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 22602, 0, 3, 9250, 21432, 3154, 3172,
                                                 10042, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 22662, 0, 3, 9322, 21522, 3208, 3238,
                                                 10198, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 22762, 0, 3, 9358, 21582, 3238, 3268,
                                                 10258, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 22862, 0, 3, 9394, 21642, 3268, 3298,
                                                 10318, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 22962, 0, 3, 9430, 21702, 3298, 3328,
                                                 10378, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23062, 0, 3, 9466, 21762, 3328, 3358,
                                                 10438, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23162, 0, 3, 9502, 21822, 3358, 3388,
                                                 10498, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23262, 0, 3, 9538, 21882, 3388, 3418,
                                                 10558, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23362, 0, 3, 9574, 21942, 3418, 3448,
                                                 10618, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23462, 0, 3, 9610, 22002, 3448, 3478,
                                                 10678, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23562, 0, 3, 9718, 22122, 3538, 3568,
                                                 10858, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23662, 0, 3, 9754, 22182, 3568, 3598,
                                                 10918, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23762, 0, 3, 9790, 22242, 3598, 3628,
                                                 10978, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23862, 0, 3, 9826, 22302, 3628, 3658,
                                                 11038, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 23962, 0, 3, 9862, 22362, 3658, 3688,
                                                 11098, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24062, 0, 3, 9898, 22422, 3688, 3718,
                                                 11158, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24162, 0, 3, 9934, 22482, 3718, 3748,
                                                 11218, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24262, 0, 3, 9970, 22542, 3748, 3778,
                                                 11278, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24362, 0, 3, 10006, 22602, 3778, 3808,
                                                 11338, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 24462, 0, 3, 21462, 21522, 10198, 22762,
                                                 3868, 3913, 11488, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 24612, 0, 3, 21522, 21582, 10258, 22862,
                                                 3913, 3958, 11578, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 24762, 0, 3, 21582, 21642, 10318, 22962,
                                                 3958, 4003, 11668, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 24912, 0, 3, 21642, 21702, 10378, 23062,
                                                 4003, 4048, 11758, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25062, 0, 3, 21702, 21762, 10438, 23162,
                                                 4048, 4093, 11848, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25212, 0, 3, 21762, 21822, 10498, 23262,
                                                 4093, 4138, 11938, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25362, 0, 3, 21822, 21882, 10558, 23362,
                                                 4138, 4183, 12028, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25512, 0, 3, 21882, 21942, 10618, 23462,
                                                 4183, 4228, 12118, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25662, 0, 3, 22062, 22122, 10858, 23662,
                                                 4318, 4363, 12298, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25812, 0, 3, 22122, 22182, 10918, 23762,
                                                 4363, 4408, 12388, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 25962, 0, 3, 22182, 22242, 10978, 23862,
                                                 4408, 4453, 12478, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26112, 0, 3, 22242, 22302, 11038, 23962,
                                                 4453, 4498, 12568, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26262, 0, 3, 22302, 22362, 11098, 24062,
                                                 4498, 4543, 12658, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26412, 0, 3, 22362, 22422, 11158, 24162,
                                                 4543, 4588, 12748, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26562, 0, 3, 22422, 22482, 11218, 24262,
                                                 4588, 4633, 12838, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26712, 0, 3, 22482, 22542, 11278, 24362,
                                                 4633, 4678, 12928, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 26862, 0, 3, 22662, 22762, 11488, 24612,
                                                 4768, 4831, 13270, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27072, 0, 3, 22762, 22862, 11578, 24762,
                                                 4831, 4894, 13396, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27282, 0, 3, 22862, 22962, 11668, 24912,
                                                 4894, 4957, 13522, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27492, 0, 3, 22962, 23062, 11758, 25062,
                                                 4957, 5020, 13648, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27702, 0, 3, 23062, 23162, 11848, 25212,
                                                 5020, 5083, 13774, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 27912, 0, 3, 23162, 23262, 11938, 25362,
                                                 5083, 5146, 13900, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28122, 0, 3, 23262, 23362, 12028, 25512,
                                                 5146, 5209, 14026, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28332, 0, 3, 23562, 23662, 12298, 25812,
                                                 5335, 5398, 14404, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28542, 0, 3, 23662, 23762, 12388, 25962,
                                                 5398, 5461, 14530, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28752, 0, 3, 23762, 23862, 12478, 26112,
                                                 5461, 5524, 14656, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28962, 0, 3, 23862, 23962, 12568, 26262,
                                                 5524, 5587, 14782, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29172, 0, 3, 23962, 24062, 12658, 26412,
                                                 5587, 5650, 14908, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29382, 0, 3, 24062, 24162, 12748, 26562,
                                                 5650, 5713, 15034, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29592, 0, 3, 24162, 24262, 12838, 26712,
                                                 5713, 5776, 15160, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 29802, 0, 3, 24462, 24612, 13270, 27072,
                                                 5902, 5986, 15454, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 30082, 0, 3, 24612, 24762, 13396, 27282,
                                                 5986, 6070, 15622, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 30362, 0, 3, 24762, 24912, 13522, 27492,
                                                 6070, 6154, 15790, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 30642, 0, 3, 24912, 25062, 13648, 27702,
                                                 6154, 6238, 15958, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 30922, 0, 3, 25062, 25212, 13774, 27912,
                                                 6238, 6322, 16126, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31202, 0, 3, 25212, 25362, 13900, 28122,
                                                 6322, 6406, 16294, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31482, 0, 3, 25662, 25812, 14404, 28542,
                                                 6574, 6658, 16630, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31762, 0, 3, 25812, 25962, 14530, 28752,
                                                 6658, 6742, 16798, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 32042, 0, 3, 25962, 26112, 14656, 28962,
                                                 6742, 6826, 16966, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 32322, 0, 3, 26112, 26262, 14782, 29172,
                                                 6826, 6910, 17134, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 32602, 0, 3, 26262, 26412, 14908, 29382,
                                                 6910, 6994, 17302, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 32882, 0, 3, 26412, 26562, 15034, 29592,
                                                 6994, 7078, 17470, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 33162, 0, 3, 26862, 27072, 15454, 30082,
                                                 7246, 7354, 18070, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 33522, 0, 3, 27072, 27282, 15622, 30362,
                                                 7354, 7462, 18286, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 33882, 0, 3, 27282, 27492, 15790, 30642,
                                                 7462, 7570, 18502, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 34242, 0, 3, 27492, 27702, 15958, 30922,
                                                 7570, 7678, 18718, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 34602, 0, 3, 27702, 27912, 16126, 31202,
                                                 7678, 7786, 18934, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 34962, 0, 3, 28332, 28542, 16630, 31762,
                                                 8002, 8110, 19582, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 35322, 0, 3, 28542, 28752, 16798, 32042,
                                                 8110, 8218, 19798, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 35682, 0, 3, 28752, 28962, 16966, 32322,
                                                 8218, 8326, 20014, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 36042, 0, 3, 28962, 29172, 17134, 32602,
                                                 8326, 8434, 20230, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 36402, 0, 3, 29172, 29382, 17302, 32882,
                                                 8434, 8542, 20446, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36762, 3, 8758, 8764, 20672, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36777, 3, 8764, 8770, 20682, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36792, 3, 8770, 8776, 20692, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36807, 3, 8776, 8782, 20702, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36822, 3, 8782, 8788, 20712, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36837, 3, 8788, 8794, 20722, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36852, 3, 8794, 8800, 20732, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36867, 3, 8800, 8806, 20742, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36882, 3, 8806, 8812, 20752, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36897, 3, 8824, 8830, 20772, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36912, 3, 8830, 8836, 20782, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36927, 3, 8836, 8842, 20792, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36942, 3, 8842, 8848, 20802, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36957, 3, 8848, 8854, 20812, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36972, 3, 8854, 8860, 20822, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 36987, 3, 8860, 8866, 20832, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 37002, 3, 8866, 8872, 20842, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 37017, 3, 8872, 8878, 20852, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 37032, 0, 3, 20662, 36762, 20892, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37077, 0, 3, 20672, 36777, 20922, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37122, 0, 3, 20682, 36792, 20952, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37167, 0, 3, 20692, 36807, 20982, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37212, 0, 3, 20702, 36822, 21012, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37257, 0, 3, 20712, 36837, 21042, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37302, 0, 3, 20722, 36852, 21072, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37347, 0, 3, 20732, 36867, 21102, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37392, 0, 3, 20742, 36882, 21132, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37437, 0, 3, 20762, 36897, 21192, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37482, 0, 3, 20772, 36912, 21222, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37527, 0, 3, 20782, 36927, 21252, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37572, 0, 3, 20792, 36942, 21282, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37617, 0, 3, 20802, 36957, 21312, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37662, 0, 3, 20812, 36972, 21342, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37707, 0, 3, 20822, 36987, 21372, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37752, 0, 3, 20832, 37002, 21402, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 37797, 0, 3, 20842, 37017, 21432, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 37842, 0, 3, 20862, 37032, 9286, 9322,
                                                 21522, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 37932, 0, 3, 20892, 37077, 9322, 9358,
                                                 21582, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38022, 0, 3, 20922, 37122, 9358, 9394,
                                                 21642, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38112, 0, 3, 20952, 37167, 9394, 9430,
                                                 21702, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38202, 0, 3, 20982, 37212, 9430, 9466,
                                                 21762, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38292, 0, 3, 21012, 37257, 9466, 9502,
                                                 21822, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38382, 0, 3, 21042, 37302, 9502, 9538,
                                                 21882, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38472, 0, 3, 21072, 37347, 9538, 9574,
                                                 21942, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38562, 0, 3, 21102, 37392, 9574, 9610,
                                                 22002, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38652, 0, 3, 21162, 37437, 9682, 9718,
                                                 22122, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38742, 0, 3, 21192, 37482, 9718, 9754,
                                                 22182, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38832, 0, 3, 21222, 37527, 9754, 9790,
                                                 22242, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 38922, 0, 3, 21252, 37572, 9790, 9826,
                                                 22302, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 39012, 0, 3, 21282, 37617, 9826, 9862,
                                                 22362, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 39102, 0, 3, 21312, 37662, 9862, 9898,
                                                 22422, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 39192, 0, 3, 21342, 37707, 9898, 9934,
                                                 22482, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 39282, 0, 3, 21372, 37752, 9934, 9970,
                                                 22542, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 39372, 0, 3, 21402, 37797, 9970, 10006,
                                                 22602, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39462, 0, 3, 21462, 37842, 10078, 10138,
                                                 22662, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39612, 0, 3, 21522, 37932, 10138, 10198,
                                                 22762, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39762, 0, 3, 21582, 38022, 10198, 10258,
                                                 22862, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 39912, 0, 3, 21642, 38112, 10258, 10318,
                                                 22962, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40062, 0, 3, 21702, 38202, 10318, 10378,
                                                 23062, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40212, 0, 3, 21762, 38292, 10378, 10438,
                                                 23162, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40362, 0, 3, 21822, 38382, 10438, 10498,
                                                 23262, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40512, 0, 3, 21882, 38472, 10498, 10558,
                                                 23362, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40662, 0, 3, 21942, 38562, 10558, 10618,
                                                 23462, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40812, 0, 3, 22062, 38652, 10738, 10798,
                                                 23562, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 40962, 0, 3, 22122, 38742, 10798, 10858,
                                                 23662, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 41112, 0, 3, 22182, 38832, 10858, 10918,
                                                 23762, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 41262, 0, 3, 22242, 38922, 10918, 10978,
                                                 23862, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 41412, 0, 3, 22302, 39012, 10978, 11038,
                                                 23962, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 41562, 0, 3, 22362, 39102, 11038, 11098,
                                                 24062, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 41712, 0, 3, 22422, 39192, 11098, 11158,
                                                 24162, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 41862, 0, 3, 22482, 39282, 11158, 11218,
                                                 24262, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 42012, 0, 3, 22542, 39372, 11218, 11278,
                                                 24362, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42162, 0, 3, 37842, 37932, 22762, 39762,
                                                 11398, 11488, 24612, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42387, 0, 3, 37932, 38022, 22862, 39912,
                                                 11488, 11578, 24762, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42612, 0, 3, 38022, 38112, 22962, 40062,
                                                 11578, 11668, 24912, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 42837, 0, 3, 38112, 38202, 23062, 40212,
                                                 11668, 11758, 25062, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 43062, 0, 3, 38202, 38292, 23162, 40362,
                                                 11758, 11848, 25212, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 43287, 0, 3, 38292, 38382, 23262, 40512,
                                                 11848, 11938, 25362, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 43512, 0, 3, 38382, 38472, 23362, 40662,
                                                 11938, 12028, 25512, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 43737, 0, 3, 38652, 38742, 23662, 41112,
                                                 12208, 12298, 25812, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 43962, 0, 3, 38742, 38832, 23762, 41262,
                                                 12298, 12388, 25962, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44187, 0, 3, 38832, 38922, 23862, 41412,
                                                 12388, 12478, 26112, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44412, 0, 3, 38922, 39012, 23962, 41562,
                                                 12478, 12568, 26262, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44637, 0, 3, 39012, 39102, 24062, 41712,
                                                 12568, 12658, 26412, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44862, 0, 3, 39102, 39192, 24162, 41862,
                                                 12658, 12748, 26562, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 45087, 0, 3, 39192, 39282, 24262, 42012,
                                                 12748, 12838, 26712, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 45312, 0, 3, 39462, 39612, 24462, 42162,
                                                 13018, 13144, 26862, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 45627, 0, 3, 39612, 39762, 24612, 42387,
                                                 13144, 13270, 27072, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 45942, 0, 3, 39762, 39912, 24762, 42612,
                                                 13270, 13396, 27282, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 46257, 0, 3, 39912, 40062, 24912, 42837,
                                                 13396, 13522, 27492, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 46572, 0, 3, 40062, 40212, 25062, 43062,
                                                 13522, 13648, 27702, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 46887, 0, 3, 40212, 40362, 25212, 43287,
                                                 13648, 13774, 27912, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 47202, 0, 3, 40362, 40512, 25362, 43512,
                                                 13774, 13900, 28122, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 47517, 0, 3, 40812, 40962, 25662, 43737,
                                                 14152, 14278, 28332, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 47832, 0, 3, 40962, 41112, 25812, 43962,
                                                 14278, 14404, 28542, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 48147, 0, 3, 41112, 41262, 25962, 44187,
                                                 14404, 14530, 28752, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 48462, 0, 3, 41262, 41412, 26112, 44412,
                                                 14530, 14656, 28962, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 48777, 0, 3, 41412, 41562, 26262, 44637,
                                                 14656, 14782, 29172, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 49092, 0, 3, 41562, 41712, 26412, 44862,
                                                 14782, 14908, 29382, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 49407, 0, 3, 41712, 41862, 26562, 45087,
                                                 14908, 15034, 29592, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 49722, 0, 3, 42162, 42387, 27072, 45942,
                                                 15286, 15454, 30082, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 50142, 0, 3, 42387, 42612, 27282, 46257,
                                                 15454, 15622, 30362, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 50562, 0, 3, 42612, 42837, 27492, 46572,
                                                 15622, 15790, 30642, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 50982, 0, 3, 42837, 43062, 27702, 46887,
                                                 15790, 15958, 30922, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 51402, 0, 3, 43062, 43287, 27912, 47202,
                                                 15958, 16126, 31202, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 51822, 0, 3, 43737, 43962, 28542, 48147,
                                                 16462, 16630, 31762, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 52242, 0, 3, 43962, 44187, 28752, 48462,
                                                 16630, 16798, 32042, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 52662, 0, 3, 44187, 44412, 28962, 48777,
                                                 16798, 16966, 32322, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 53082, 0, 3, 44412, 44637, 29172, 49092,
                                                 16966, 17134, 32602, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 53502, 0, 3, 44637, 44862, 29382, 49407,
                                                 17134, 17302, 32882, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 53922, 0, 3, 45312, 45627, 29802, 49722,
                                                 17638, 17854, 33162, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 54462, 0, 3, 45627, 45942, 30082, 50142,
                                                 17854, 18070, 33522, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 55002, 0, 3, 45942, 46257, 30362, 50562,
                                                 18070, 18286, 33882, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 55542, 0, 3, 46257, 46572, 30642, 50982,
                                                 18286, 18502, 34242, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 56082, 0, 3, 46572, 46887, 30922, 51402,
                                                 18502, 18718, 34602, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 56622, 0, 3, 47517, 47832, 31482, 51822,
                                                 19150, 19366, 34962, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 57162, 0, 3, 47832, 48147, 31762, 52242,
                                                 19366, 19582, 35322, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 57702, 0, 3, 48147, 48462, 32042, 52662,
                                                 19582, 19798, 35682, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 58242, 0, 3, 48462, 48777, 32322, 53082,
                                                 19798, 20014, 36042, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 58782, 0, 3, 48777, 49092, 32602, 53502,
                                                 20014, 20230, 36402, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59322, 3, 20662, 20672, 36777, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59343, 3, 20672, 20682, 36792, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59364, 3, 20682, 20692, 36807, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59385, 3, 20692, 20702, 36822, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59406, 3, 20702, 20712, 36837, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59427, 3, 20712, 20722, 36852, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59448, 3, 20722, 20732, 36867, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59469, 3, 20732, 20742, 36882, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59490, 3, 20762, 20772, 36912, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59511, 3, 20772, 20782, 36927, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59532, 3, 20782, 20792, 36942, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59553, 3, 20792, 20802, 36957, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59574, 3, 20802, 20812, 36972, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59595, 3, 20812, 20822, 36987, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59616, 3, 20822, 20832, 37002, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 59637, 3, 20832, 20842, 37017, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 59658, 0, 3, 36762, 59322, 37077, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 59721, 0, 3, 36777, 59343, 37122, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 59784, 0, 3, 36792, 59364, 37167, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 59847, 0, 3, 36807, 59385, 37212, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 59910, 0, 3, 36822, 59406, 37257, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 59973, 0, 3, 36837, 59427, 37302, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60036, 0, 3, 36852, 59448, 37347, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60099, 0, 3, 36867, 59469, 37392, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60162, 0, 3, 36897, 59490, 37482, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60225, 0, 3, 36912, 59511, 37527, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60288, 0, 3, 36927, 59532, 37572, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60351, 0, 3, 36942, 59553, 37617, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60414, 0, 3, 36957, 59574, 37662, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60477, 0, 3, 36972, 59595, 37707, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60540, 0, 3, 36987, 59616, 37752, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 60603, 0, 3, 37002, 59637, 37797, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 60666, 0, 3, 37032, 59658, 21462, 21522,
                                                 37932, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 60792, 0, 3, 37077, 59721, 21522, 21582,
                                                 38022, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 60918, 0, 3, 37122, 59784, 21582, 21642,
                                                 38112, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 61044, 0, 3, 37167, 59847, 21642, 21702,
                                                 38202, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 61170, 0, 3, 37212, 59910, 21702, 21762,
                                                 38292, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 61296, 0, 3, 37257, 59973, 21762, 21822,
                                                 38382, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 61422, 0, 3, 37302, 60036, 21822, 21882,
                                                 38472, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 61548, 0, 3, 37347, 60099, 21882, 21942,
                                                 38562, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 61674, 0, 3, 37437, 60162, 22062, 22122,
                                                 38742, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 61800, 0, 3, 37482, 60225, 22122, 22182,
                                                 38832, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 61926, 0, 3, 37527, 60288, 22182, 22242,
                                                 38922, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 62052, 0, 3, 37572, 60351, 22242, 22302,
                                                 39012, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 62178, 0, 3, 37617, 60414, 22302, 22362,
                                                 39102, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 62304, 0, 3, 37662, 60477, 22362, 22422,
                                                 39192, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 62430, 0, 3, 37707, 60540, 22422, 22482,
                                                 39282, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 62556, 0, 3, 37752, 60603, 22482, 22542,
                                                 39372, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 62682, 0, 3, 37932, 60792, 22662, 22762,
                                                 39762, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 62892, 0, 3, 38022, 60918, 22762, 22862,
                                                 39912, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 63102, 0, 3, 38112, 61044, 22862, 22962,
                                                 40062, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 63312, 0, 3, 38202, 61170, 22962, 23062,
                                                 40212, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 63522, 0, 3, 38292, 61296, 23062, 23162,
                                                 40362, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 63732, 0, 3, 38382, 61422, 23162, 23262,
                                                 40512, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 63942, 0, 3, 38472, 61548, 23262, 23362,
                                                 40662, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 64152, 0, 3, 38742, 61800, 23562, 23662,
                                                 41112, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 64362, 0, 3, 38832, 61926, 23662, 23762,
                                                 41262, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 64572, 0, 3, 38922, 62052, 23762, 23862,
                                                 41412, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 64782, 0, 3, 39012, 62178, 23862, 23962,
                                                 41562, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 64992, 0, 3, 39102, 62304, 23962, 24062,
                                                 41712, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 65202, 0, 3, 39192, 62430, 24062, 24162,
                                                 41862, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 65412, 0, 3, 39282, 62556, 24162, 24262,
                                                 42012, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 65622, 0, 3, 60666, 60792, 39762, 62892,
                                                 24462, 24612, 42387, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 65937, 0, 3, 60792, 60918, 39912, 63102,
                                                 24612, 24762, 42612, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 66252, 0, 3, 60918, 61044, 40062, 63312,
                                                 24762, 24912, 42837, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 66567, 0, 3, 61044, 61170, 40212, 63522,
                                                 24912, 25062, 43062, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 66882, 0, 3, 61170, 61296, 40362, 63732,
                                                 25062, 25212, 43287, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 67197, 0, 3, 61296, 61422, 40512, 63942,
                                                 25212, 25362, 43512, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 67512, 0, 3, 61674, 61800, 41112, 64362,
                                                 25662, 25812, 43962, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 67827, 0, 3, 61800, 61926, 41262, 64572,
                                                 25812, 25962, 44187, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 68142, 0, 3, 61926, 62052, 41412, 64782,
                                                 25962, 26112, 44412, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 68457, 0, 3, 62052, 62178, 41562, 64992,
                                                 26112, 26262, 44637, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 68772, 0, 3, 62178, 62304, 41712, 65202,
                                                 26262, 26412, 44862, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 69087, 0, 3, 62304, 62430, 41862, 65412,
                                                 26412, 26562, 45087, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 69402, 0, 3, 62682, 62892, 42387, 65937,
                                                 26862, 27072, 45942, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 69843, 0, 3, 62892, 63102, 42612, 66252,
                                                 27072, 27282, 46257, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 70284, 0, 3, 63102, 63312, 42837, 66567,
                                                 27282, 27492, 46572, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 70725, 0, 3, 63312, 63522, 43062, 66882,
                                                 27492, 27702, 46887, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 71166, 0, 3, 63522, 63732, 43287, 67197,
                                                 27702, 27912, 47202, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 71607, 0, 3, 64152, 64362, 43962, 67827,
                                                 28332, 28542, 48147, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 72048, 0, 3, 64362, 64572, 44187, 68142,
                                                 28542, 28752, 48462, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 72489, 0, 3, 64572, 64782, 44412, 68457,
                                                 28752, 28962, 48777, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 72930, 0, 3, 64782, 64992, 44637, 68772,
                                                 28962, 29172, 49092, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 73371, 0, 3, 64992, 65202, 44862, 69087,
                                                 29172, 29382, 49407, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 73812, 0, 3, 65622, 65937, 45942, 69843,
                                                 29802, 30082, 50142, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 74400, 0, 3, 65937, 66252, 46257, 70284,
                                                 30082, 30362, 50562, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 74988, 0, 3, 66252, 66567, 46572, 70725,
                                                 30362, 30642, 50982, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 75576, 0, 3, 66567, 66882, 46887, 71166,
                                                 30642, 30922, 51402, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 76164, 0, 3, 67512, 67827, 48147, 72048,
                                                 31482, 31762, 52242, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 76752, 0, 3, 67827, 68142, 48462, 72489,
                                                 31762, 32042, 52662, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 77340, 0, 3, 68142, 68457, 48777, 72930,
                                                 32042, 32322, 53082, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 77928, 0, 3, 68457, 68772, 49092, 73371,
                                                 32322, 32602, 53502, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 78516, 0, 3, 69402, 69843, 50142, 74400,
                                                 33162, 33522, 55002, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 79272, 0, 3, 69843, 70284, 50562, 74988,
                                                 33522, 33882, 55542, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 80028, 0, 3, 70284, 70725, 50982, 75576,
                                                 33882, 34242, 56082, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 80784, 0, 3, 71607, 72048, 52242, 76752,
                                                 34962, 35322, 57702, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 81540, 0, 3, 72048, 72489, 52662, 77340,
                                                 35322, 35682, 58242, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 82296, 0, 3, 72489, 72930, 53082, 77928,
                                                 35682, 36042, 58782, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83052, 3, 36762, 36777, 59343, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83080, 3, 36777, 36792, 59364, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83108, 3, 36792, 36807, 59385, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83136, 3, 36807, 36822, 59406, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83164, 3, 36822, 36837, 59427, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83192, 3, 36837, 36852, 59448, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83220, 3, 36852, 36867, 59469, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83248, 3, 36897, 36912, 59511, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83276, 3, 36912, 36927, 59532, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83304, 3, 36927, 36942, 59553, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83332, 3, 36942, 36957, 59574, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83360, 3, 36957, 36972, 59595, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83388, 3, 36972, 36987, 59616, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 83416, 3, 36987, 37002, 59637, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 83444, 0, 3, 59322, 83052, 59721, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83528, 0, 3, 59343, 83080, 59784, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83612, 0, 3, 59364, 83108, 59847, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83696, 0, 3, 59385, 83136, 59910, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83780, 0, 3, 59406, 83164, 59973, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83864, 0, 3, 59427, 83192, 60036, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83948, 0, 3, 59448, 83220, 60099, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 84032, 0, 3, 59490, 83248, 60225, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 84116, 0, 3, 59511, 83276, 60288, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 84200, 0, 3, 59532, 83304, 60351, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 84284, 0, 3, 59553, 83332, 60414, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 84368, 0, 3, 59574, 83360, 60477, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 84452, 0, 3, 59595, 83388, 60540, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 84536, 0, 3, 59616, 83416, 60603, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 84620, 0, 3, 59658, 83444, 37842, 37932,
                                                 60792, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84788, 0, 3, 59721, 83528, 37932, 38022,
                                                 60918, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84956, 0, 3, 59784, 83612, 38022, 38112,
                                                 61044, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85124, 0, 3, 59847, 83696, 38112, 38202,
                                                 61170, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85292, 0, 3, 59910, 83780, 38202, 38292,
                                                 61296, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85460, 0, 3, 59973, 83864, 38292, 38382,
                                                 61422, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85628, 0, 3, 60036, 83948, 38382, 38472,
                                                 61548, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85796, 0, 3, 60162, 84032, 38652, 38742,
                                                 61800, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85964, 0, 3, 60225, 84116, 38742, 38832,
                                                 61926, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 86132, 0, 3, 60288, 84200, 38832, 38922,
                                                 62052, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 86300, 0, 3, 60351, 84284, 38922, 39012,
                                                 62178, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 86468, 0, 3, 60414, 84368, 39012, 39102,
                                                 62304, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 86636, 0, 3, 60477, 84452, 39102, 39192,
                                                 62430, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 86804, 0, 3, 60540, 84536, 39192, 39282,
                                                 62556, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 86972, 0, 3, 60666, 84620, 39462, 39612,
                                                 62682, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 87252, 0, 3, 60792, 84788, 39612, 39762,
                                                 62892, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 87532, 0, 3, 60918, 84956, 39762, 39912,
                                                 63102, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 87812, 0, 3, 61044, 85124, 39912, 40062,
                                                 63312, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 88092, 0, 3, 61170, 85292, 40062, 40212,
                                                 63522, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 88372, 0, 3, 61296, 85460, 40212, 40362,
                                                 63732, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 88652, 0, 3, 61422, 85628, 40362, 40512,
                                                 63942, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 88932, 0, 3, 61674, 85796, 40812, 40962,
                                                 64152, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 89212, 0, 3, 61800, 85964, 40962, 41112,
                                                 64362, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 89492, 0, 3, 61926, 86132, 41112, 41262,
                                                 64572, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 89772, 0, 3, 62052, 86300, 41262, 41412,
                                                 64782, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 90052, 0, 3, 62178, 86468, 41412, 41562,
                                                 64992, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 90332, 0, 3, 62304, 86636, 41562, 41712,
                                                 65202, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 90612, 0, 3, 62430, 86804, 41712, 41862,
                                                 65412, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 90892, 0, 3, 84620, 84788, 62892, 87532,
                                                 42162, 42387, 65937, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 91312, 0, 3, 84788, 84956, 63102, 87812,
                                                 42387, 42612, 66252, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 91732, 0, 3, 84956, 85124, 63312, 88092,
                                                 42612, 42837, 66567, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 92152, 0, 3, 85124, 85292, 63522, 88372,
                                                 42837, 43062, 66882, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 92572, 0, 3, 85292, 85460, 63732, 88652,
                                                 43062, 43287, 67197, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 92992, 0, 3, 85796, 85964, 64362, 89492,
                                                 43737, 43962, 67827, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 93412, 0, 3, 85964, 86132, 64572, 89772,
                                                 43962, 44187, 68142, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 93832, 0, 3, 86132, 86300, 64782, 90052,
                                                 44187, 44412, 68457, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 94252, 0, 3, 86300, 86468, 64992, 90332,
                                                 44412, 44637, 68772, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 94672, 0, 3, 86468, 86636, 65202, 90612,
                                                 44637, 44862, 69087, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 95092, 0, 3, 86972, 87252, 65622, 90892,
                                                 45312, 45627, 69402, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 95680, 0, 3, 87252, 87532, 65937, 91312,
                                                 45627, 45942, 69843, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 96268, 0, 3, 87532, 87812, 66252, 91732,
                                                 45942, 46257, 70284, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 96856, 0, 3, 87812, 88092, 66567, 92152,
                                                 46257, 46572, 70725, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 97444, 0, 3, 88092, 88372, 66882, 92572,
                                                 46572, 46887, 71166, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 98032, 0, 3, 88932, 89212, 67512, 92992,
                                                 47517, 47832, 71607, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 98620, 0, 3, 89212, 89492, 67827, 93412,
                                                 47832, 48147, 72048, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 99208, 0, 3, 89492, 89772, 68142, 93832,
                                                 48147, 48462, 72489, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 99796, 0, 3, 89772, 90052, 68457, 94252,
                                                 48462, 48777, 72930, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 100384, 0, 3, 90052, 90332, 68772,
                                                 94672, 48777, 49092, 73371, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 100972, 0, 3, 90892, 91312, 69843,
                                                 96268, 49722, 50142, 74400, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 101756, 0, 3, 91312, 91732, 70284,
                                                 96856, 50142, 50562, 74988, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 102540, 0, 3, 91732, 92152, 70725,
                                                 97444, 50562, 50982, 75576, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 103324, 0, 3, 92992, 93412, 72048,
                                                 99208, 51822, 52242, 76752, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 104108, 0, 3, 93412, 93832, 72489,
                                                 99796, 52242, 52662, 77340, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 104892, 0, 3, 93832, 94252, 72930,
                                                 100384, 52662, 53082, 77928, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 105676, 0, 3, 95092, 95680, 73812,
                                                 100972, 53922, 54462, 78516, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 106684, 0, 3, 95680, 96268, 74400,
                                                 101756, 54462, 55002, 79272, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 107692, 0, 3, 96268, 96856, 74988,
                                                 102540, 55002, 55542, 80028, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 108700, 0, 3, 98032, 98620, 76164,
                                                 103324, 56622, 57162, 80784, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 109708, 0, 3, 98620, 99208, 76752,
                                                 104108, 57162, 57702, 81540, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 110716, 0, 3, 99208, 99796, 77340,
                                                 104892, 57702, 58242, 82296, ncols, alpha, beta,
                                                 p);

            compute_prim_sk_electron_repulsion_0(buffer, 111724, 3, 59322, 59343, 83080, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 111760, 3, 59343, 59364, 83108, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 111796, 3, 59364, 59385, 83136, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 111832, 3, 59385, 59406, 83164, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 111868, 3, 59406, 59427, 83192, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 111904, 3, 59427, 59448, 83220, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 111940, 3, 59490, 59511, 83276, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 111976, 3, 59511, 59532, 83304, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112012, 3, 59532, 59553, 83332, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112048, 3, 59553, 59574, 83360, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112084, 3, 59574, 59595, 83388, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 112120, 3, 59595, 59616, 83416, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112156, 0, 3, 83052, 111724, 83528,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112264, 0, 3, 83080, 111760, 83612,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112372, 0, 3, 83108, 111796, 83696,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112480, 0, 3, 83136, 111832, 83780,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112588, 0, 3, 83164, 111868, 83864,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112696, 0, 3, 83192, 111904, 83948,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112804, 0, 3, 83248, 111940, 84116,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 112912, 0, 3, 83276, 111976, 84200,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 113020, 0, 3, 83304, 112012, 84284,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 113128, 0, 3, 83332, 112048, 84368,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 113236, 0, 3, 83360, 112084, 84452,
                                                 ncols, p);

            compute_prim_pk_electron_repulsion_0(buffer, 113344, 0, 3, 83388, 112120, 84536,
                                                 ncols, p);

            compute_prim_dk_electron_repulsion_0(buffer, 113452, 0, 3, 83444, 112156, 60666,
                                                 60792, 84788, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 113668, 0, 3, 83528, 112264, 60792,
                                                 60918, 84956, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 113884, 0, 3, 83612, 112372, 60918,
                                                 61044, 85124, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 114100, 0, 3, 83696, 112480, 61044,
                                                 61170, 85292, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 114316, 0, 3, 83780, 112588, 61170,
                                                 61296, 85460, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 114532, 0, 3, 83864, 112696, 61296,
                                                 61422, 85628, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 114748, 0, 3, 84032, 112804, 61674,
                                                 61800, 85964, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 114964, 0, 3, 84116, 112912, 61800,
                                                 61926, 86132, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 115180, 0, 3, 84200, 113020, 61926,
                                                 62052, 86300, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 115396, 0, 3, 84284, 113128, 62052,
                                                 62178, 86468, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 115612, 0, 3, 84368, 113236, 62178,
                                                 62304, 86636, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 115828, 0, 3, 84452, 113344, 62304,
                                                 62430, 86804, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 116044, 0, 3, 84788, 113668, 62682,
                                                 62892, 87532, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 116404, 0, 3, 84956, 113884, 62892,
                                                 63102, 87812, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 116764, 0, 3, 85124, 114100, 63102,
                                                 63312, 88092, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 117124, 0, 3, 85292, 114316, 63312,
                                                 63522, 88372, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 117484, 0, 3, 85460, 114532, 63522,
                                                 63732, 88652, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 117844, 0, 3, 85964, 114964, 64152,
                                                 64362, 89492, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 118204, 0, 3, 86132, 115180, 64362,
                                                 64572, 89772, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 118564, 0, 3, 86300, 115396, 64572,
                                                 64782, 90052, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 118924, 0, 3, 86468, 115612, 64782,
                                                 64992, 90332, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 119284, 0, 3, 86636, 115828, 64992,
                                                 65202, 90612, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 119644, 0, 3, 113452, 113668, 87532,
                                                 116404, 65622, 65937, 91312, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 120184, 0, 3, 113668, 113884, 87812,
                                                 116764, 65937, 66252, 91732, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 120724, 0, 3, 113884, 114100, 88092,
                                                 117124, 66252, 66567, 92152, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 121264, 0, 3, 114100, 114316, 88372,
                                                 117484, 66567, 66882, 92572, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 121804, 0, 3, 114748, 114964, 89492,
                                                 118204, 67512, 67827, 93412, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 122344, 0, 3, 114964, 115180, 89772,
                                                 118564, 67827, 68142, 93832, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 122884, 0, 3, 115180, 115396, 90052,
                                                 118924, 68142, 68457, 94252, ncols, alpha, beta,
                                                 p);

            compute_prim_gk_electron_repulsion_0(buffer, 123424, 0, 3, 115396, 115612, 90332,
                                                 119284, 68457, 68772, 94672, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 123964, 0, 3, 116044, 116404, 91312,
                                                 120184, 69402, 69843, 96268, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 124720, 0, 3, 116404, 116764, 91732,
                                                 120724, 69843, 70284, 96856, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 125476, 0, 3, 116764, 117124, 92152,
                                                 121264, 70284, 70725, 97444, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 126232, 0, 3, 117844, 118204, 93412,
                                                 122344, 71607, 72048, 99208, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 126988, 0, 3, 118204, 118564, 93832,
                                                 122884, 72048, 72489, 99796, ncols, alpha, beta,
                                                 p);

            compute_prim_hk_electron_repulsion_0(buffer, 127744, 0, 3, 118564, 118924, 94252,
                                                 123424, 72489, 72930, 100384, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 128500, 0, 3, 119644, 120184, 96268,
                                                 124720, 73812, 74400, 101756, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 129508, 0, 3, 120184, 120724, 96856,
                                                 125476, 74400, 74988, 102540, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 130516, 0, 3, 121804, 122344, 99208,
                                                 126988, 76164, 76752, 104108, ncols, alpha,
                                                 beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 131524, 0, 3, 122344, 122884, 99796,
                                                 127744, 76752, 77340, 104892, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 132532, 0, 3, 123964, 124720, 101756,
                                                 129508, 78516, 79272, 107692, ncols, alpha,
                                                 beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 133828, 0, 3, 126232, 126988, 104108,
                                                 131524, 80784, 81540, 110716, ncols, alpha,
                                                 beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135124, 3, 83052, 83080, 111760, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135169, 3, 83080, 83108, 111796, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135214, 3, 83108, 83136, 111832, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135259, 3, 83136, 83164, 111868, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135304, 3, 83164, 83192, 111904, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135349, 3, 83248, 83276, 111976, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135394, 3, 83276, 83304, 112012, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135439, 3, 83304, 83332, 112048, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135484, 3, 83332, 83360, 112084, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 135529, 3, 83360, 83388, 112120, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 135574, 0, 3, 111724, 135124, 112264,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 135709, 0, 3, 111760, 135169, 112372,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 135844, 0, 3, 111796, 135214, 112480,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 135979, 0, 3, 111832, 135259, 112588,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 136114, 0, 3, 111868, 135304, 112696,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 136249, 0, 3, 111940, 135349, 112912,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 136384, 0, 3, 111976, 135394, 113020,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 136519, 0, 3, 112012, 135439, 113128,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 136654, 0, 3, 112048, 135484, 113236,
                                                 ncols, p);

            compute_prim_pl_electron_repulsion_0(buffer, 136789, 0, 3, 112084, 135529, 113344,
                                                 ncols, p);

            compute_prim_dl_electron_repulsion_0(buffer, 136924, 0, 3, 112156, 135574, 84620,
                                                 84788, 113668, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 137194, 0, 3, 112264, 135709, 84788,
                                                 84956, 113884, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 137464, 0, 3, 112372, 135844, 84956,
                                                 85124, 114100, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 137734, 0, 3, 112480, 135979, 85124,
                                                 85292, 114316, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 138004, 0, 3, 112588, 136114, 85292,
                                                 85460, 114532, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 138274, 0, 3, 112804, 136249, 85796,
                                                 85964, 114964, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 138544, 0, 3, 112912, 136384, 85964,
                                                 86132, 115180, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 138814, 0, 3, 113020, 136519, 86132,
                                                 86300, 115396, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 139084, 0, 3, 113128, 136654, 86300,
                                                 86468, 115612, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 139354, 0, 3, 113236, 136789, 86468,
                                                 86636, 115828, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 139624, 0, 3, 113452, 136924, 86972,
                                                 87252, 116044, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 140074, 0, 3, 113668, 137194, 87252,
                                                 87532, 116404, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 140524, 0, 3, 113884, 137464, 87532,
                                                 87812, 116764, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 140974, 0, 3, 114100, 137734, 87812,
                                                 88092, 117124, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 141424, 0, 3, 114316, 138004, 88092,
                                                 88372, 117484, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 141874, 0, 3, 114748, 138274, 88932,
                                                 89212, 117844, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 142324, 0, 3, 114964, 138544, 89212,
                                                 89492, 118204, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 142774, 0, 3, 115180, 138814, 89492,
                                                 89772, 118564, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 143224, 0, 3, 115396, 139084, 89772,
                                                 90052, 118924, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 143674, 0, 3, 115612, 139354, 90052,
                                                 90332, 119284, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 144124, 0, 3, 136924, 137194, 116404,
                                                 140524, 90892, 91312, 120184, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 144799, 0, 3, 137194, 137464, 116764,
                                                 140974, 91312, 91732, 120724, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 145474, 0, 3, 137464, 137734, 117124,
                                                 141424, 91732, 92152, 121264, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 146149, 0, 3, 138274, 138544, 118204,
                                                 142774, 92992, 93412, 122344, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 146824, 0, 3, 138544, 138814, 118564,
                                                 143224, 93412, 93832, 122884, ncols, alpha,
                                                 beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 147499, 0, 3, 138814, 139084, 118924,
                                                 143674, 93832, 94252, 123424, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 148174, 0, 3, 139624, 140074, 119644,
                                                 144124, 95092, 95680, 123964, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 149119, 0, 3, 140074, 140524, 120184,
                                                 144799, 95680, 96268, 124720, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 150064, 0, 3, 140524, 140974, 120724,
                                                 145474, 96268, 96856, 125476, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 151009, 0, 3, 141874, 142324, 121804,
                                                 146149, 98032, 98620, 126232, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 151954, 0, 3, 142324, 142774, 122344,
                                                 146824, 98620, 99208, 126988, ncols, alpha,
                                                 beta, p);

            compute_prim_hl_electron_repulsion_0(buffer, 152899, 0, 3, 142774, 143224, 122884,
                                                 147499, 99208, 99796, 127744, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 153844, 0, 3, 144124, 144799, 124720,
                                                 150064, 100972, 101756, 129508, ncols, alpha,
                                                 beta, p);

            compute_prim_il_electron_repulsion_0(buffer, 155104, 0, 3, 146149, 146824, 126988,
                                                 152899, 103324, 104108, 131524, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 156364, 0, 3, 148174, 149119, 128500,
                                                 153844, 105676, 106684, 132532, ncols, alpha,
                                                 beta, p);

            compute_prim_kl_electron_repulsion_0(buffer, 157984, 0, 3, 151009, 151954, 130516,
                                                 155104, 108700, 109708, 133828, ncols, alpha,
                                                 beta, p);

            simdgeo::geom_i_x(buffer, 159604, 151009, 157984, 1, 45, ncols, alpha);

            simdgeo::geom_i_y(buffer, 160864, 151009, 157984, 1, 45, ncols, alpha);

            simdgeo::geom_i_z(buffer, 162124, 151009, 157984, 1, 45, ncols, alpha);

            simdgeo::geom_i_x(buffer, 163384, 148174, 156364, 1, 45, ncols, alpha);

            simdgeo::geom_i_y(buffer, 164644, 148174, 156364, 1, 45, ncols, alpha);

            simdgeo::geom_i_z(buffer, 165904, 148174, 156364, 1, 45, ncols, alpha);

            simdfunc::contract_primitives(buffer, 167164, 163384, 3780, ncols);

            simdfunc::contract_primitives(buffer, 170944, 159604, 3780, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 174724, 170944, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 174724, 17, nmax);

    simdtrf::transform_l_inner(buffer, 174724, 172204, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 221 * nvalues, nvalues, buffer, 174724, 17, nmax);

    simdtrf::transform_l_inner(buffer, 174724, 173464, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 442 * nvalues, nvalues, buffer, 174724, 17, nmax);

    simdtrf::transform_l_inner(buffer, 174724, 167164, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 663 * nvalues, nvalues, buffer, 174724, 17, nmax);

    simdtrf::transform_l_inner(buffer, 174724, 168424, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 884 * nvalues, nvalues, buffer, 174724, 17, nmax);

    simdtrf::transform_l_inner(buffer, 174724, 169684, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 1105 * nvalues, nvalues, buffer, 174724, 17, nmax);
}

}  // namespace simdt2ceri
