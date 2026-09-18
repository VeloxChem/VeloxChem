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


#include "SimdElectronRepulsionGeom10RsRecII.hpp"

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
#include "SimdGeometryI1.hpp"
#include "SimdTransformI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ii_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ii_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 81176, 76108, 4704, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 2008, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2011, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2014, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2017, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2020, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2023, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2026, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2029, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2032, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2035, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2038, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2041, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2044, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2047, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2050, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2053, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2056, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2059, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2062, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2065, 3, 33, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2068, 3, 9, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2077, 3, 10, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2086, 3, 11, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2095, 3, 12, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2104, 3, 13, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2113, 3, 14, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2122, 3, 15, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2131, 3, 16, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2140, 3, 17, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2149, 3, 18, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2158, 3, 23, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2167, 3, 24, 85, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2176, 3, 25, 88, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2185, 3, 26, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2194, 3, 27, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2203, 3, 28, 97, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2212, 3, 29, 100, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2221, 3, 30, 103, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2230, 3, 31, 106, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2239, 3, 32, 109, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2248, 0, 3, 40, 2068, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2266, 0, 3, 43, 2077, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2284, 0, 3, 46, 2086, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2302, 0, 3, 49, 2095, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2320, 0, 3, 52, 2104, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2338, 0, 3, 55, 2113, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2356, 0, 3, 58, 2122, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2374, 0, 3, 61, 2131, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2392, 0, 3, 64, 2140, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2410, 0, 3, 67, 2149, 172, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2428, 0, 3, 79, 2158, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2446, 0, 3, 82, 2167, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2464, 0, 3, 85, 2176, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2482, 0, 3, 88, 2185, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2500, 0, 3, 91, 2194, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2518, 0, 3, 94, 2203, 214, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2536, 0, 3, 97, 2212, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2554, 0, 3, 100, 2221, 226, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2572, 0, 3, 103, 2230, 232, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2590, 0, 3, 106, 2239, 238, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2608, 0, 3, 118, 2266, 264, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2638, 0, 3, 124, 2284, 274, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2668, 0, 3, 130, 2302, 284, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2698, 0, 3, 136, 2320, 294, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2728, 0, 3, 142, 2338, 304, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2758, 0, 3, 148, 2356, 314, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2788, 0, 3, 154, 2374, 324, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2818, 0, 3, 160, 2392, 334, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2848, 0, 3, 166, 2410, 344, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2878, 0, 3, 184, 2446, 374, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2908, 0, 3, 190, 2464, 384, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2938, 0, 3, 196, 2482, 394, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2968, 0, 3, 202, 2500, 404, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2998, 0, 3, 208, 2518, 414, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3028, 0, 3, 214, 2536, 424, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3058, 0, 3, 220, 2554, 434, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3088, 0, 3, 226, 2572, 444, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3118, 0, 3, 232, 2590, 454, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3148, 0, 3, 264, 2638, 479, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3193, 0, 3, 274, 2668, 494, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3238, 0, 3, 284, 2698, 509, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3283, 0, 3, 294, 2728, 524, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3328, 0, 3, 304, 2758, 539, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3373, 0, 3, 314, 2788, 554, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3418, 0, 3, 324, 2818, 569, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3463, 0, 3, 334, 2848, 584, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3508, 0, 3, 374, 2908, 614, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3553, 0, 3, 384, 2938, 629, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3598, 0, 3, 394, 2968, 644, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3643, 0, 3, 404, 2998, 659, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3688, 0, 3, 414, 3028, 674, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3733, 0, 3, 424, 3058, 689, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3778, 0, 3, 434, 3088, 704, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3823, 0, 3, 444, 3118, 719, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3868, 0, 3, 479, 3193, 776, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3931, 0, 3, 494, 3238, 797, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3994, 0, 3, 509, 3283, 818, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4057, 0, 3, 524, 3328, 839, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4120, 0, 3, 539, 3373, 860, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4183, 0, 3, 554, 3418, 881, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4246, 0, 3, 569, 3463, 902, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4309, 0, 3, 614, 3553, 965, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4372, 0, 3, 629, 3598, 986, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4435, 0, 3, 644, 3643, 1007, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4498, 0, 3, 659, 3688, 1028, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4561, 0, 3, 674, 3733, 1049, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4624, 0, 3, 689, 3778, 1070, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4687, 0, 3, 704, 3823, 1091, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4750, 0, 3, 776, 3931, 1140, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4834, 0, 3, 797, 3994, 1168, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4918, 0, 3, 818, 4057, 1196, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5002, 0, 3, 839, 4120, 1224, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5086, 0, 3, 860, 4183, 1252, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5170, 0, 3, 881, 4246, 1280, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5254, 0, 3, 965, 4372, 1336, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5338, 0, 3, 986, 4435, 1364, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5422, 0, 3, 1007, 4498, 1392, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5506, 0, 3, 1028, 4561, 1420, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5590, 0, 3, 1049, 4624, 1448, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5674, 0, 3, 1070, 4687, 1476, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5758, 0, 3, 1140, 4834, 1576, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5866, 0, 3, 1168, 4918, 1612, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 5974, 0, 3, 1196, 5002, 1648, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6082, 0, 3, 1224, 5086, 1684, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6190, 0, 3, 1252, 5170, 1720, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6298, 0, 3, 1336, 5338, 1828, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6406, 0, 3, 1364, 5422, 1864, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6514, 0, 3, 1392, 5506, 1900, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6622, 0, 3, 1420, 5590, 1936, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6730, 0, 3, 1448, 5674, 1972, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 6838, 3, 9, 10, 2011, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6844, 3, 10, 11, 2014, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6850, 3, 11, 12, 2017, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6856, 3, 12, 13, 2020, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6862, 3, 13, 14, 2023, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6868, 3, 14, 15, 2026, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6874, 3, 15, 16, 2029, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6880, 3, 16, 17, 2032, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6886, 3, 17, 18, 2035, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6892, 3, 23, 24, 2041, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6898, 3, 24, 25, 2044, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6904, 3, 25, 26, 2047, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6910, 3, 26, 27, 2050, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6916, 3, 27, 28, 2053, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6922, 3, 28, 29, 2056, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6928, 3, 29, 30, 2059, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6934, 3, 30, 31, 2062, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6940, 3, 31, 32, 2065, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 6946, 0, 3, 2008, 6838, 2077, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6964, 0, 3, 2011, 6844, 2086, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6982, 0, 3, 2014, 6850, 2095, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7000, 0, 3, 2017, 6856, 2104, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7018, 0, 3, 2020, 6862, 2113, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7036, 0, 3, 2023, 6868, 2122, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7054, 0, 3, 2026, 6874, 2131, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7072, 0, 3, 2029, 6880, 2140, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7090, 0, 3, 2032, 6886, 2149, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7108, 0, 3, 2038, 6892, 2167, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7126, 0, 3, 2041, 6898, 2176, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7144, 0, 3, 2044, 6904, 2185, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7162, 0, 3, 2047, 6910, 2194, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7180, 0, 3, 2050, 6916, 2203, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7198, 0, 3, 2053, 6922, 2212, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7216, 0, 3, 2056, 6928, 2221, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7234, 0, 3, 2059, 6934, 2230, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 7252, 0, 3, 2062, 6940, 2239, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 7270, 0, 3, 2068, 6946, 112, 118, 2266,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7306, 0, 3, 2077, 6964, 118, 124, 2284,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7342, 0, 3, 2086, 6982, 124, 130, 2302,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7378, 0, 3, 2095, 7000, 130, 136, 2320,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7414, 0, 3, 2104, 7018, 136, 142, 2338,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7450, 0, 3, 2113, 7036, 142, 148, 2356,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7486, 0, 3, 2122, 7054, 148, 154, 2374,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7522, 0, 3, 2131, 7072, 154, 160, 2392,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7558, 0, 3, 2140, 7090, 160, 166, 2410,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7594, 0, 3, 2158, 7108, 178, 184, 2446,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7630, 0, 3, 2167, 7126, 184, 190, 2464,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7666, 0, 3, 2176, 7144, 190, 196, 2482,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7702, 0, 3, 2185, 7162, 196, 202, 2500,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7738, 0, 3, 2194, 7180, 202, 208, 2518,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7774, 0, 3, 2203, 7198, 208, 214, 2536,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7810, 0, 3, 2212, 7216, 214, 220, 2554,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7846, 0, 3, 2221, 7234, 220, 226, 2572,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7882, 0, 3, 2230, 7252, 226, 232, 2590,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7918, 0, 3, 2248, 7270, 244, 254, 2608,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7978, 0, 3, 2266, 7306, 254, 264, 2638,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8038, 0, 3, 2284, 7342, 264, 274, 2668,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8098, 0, 3, 2302, 7378, 274, 284, 2698,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8158, 0, 3, 2320, 7414, 284, 294, 2728,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8218, 0, 3, 2338, 7450, 294, 304, 2758,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8278, 0, 3, 2356, 7486, 304, 314, 2788,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8338, 0, 3, 2374, 7522, 314, 324, 2818,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8398, 0, 3, 2392, 7558, 324, 334, 2848,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8458, 0, 3, 2428, 7594, 354, 364, 2878,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8518, 0, 3, 2446, 7630, 364, 374, 2908,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8578, 0, 3, 2464, 7666, 374, 384, 2938,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8638, 0, 3, 2482, 7702, 384, 394, 2968,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8698, 0, 3, 2500, 7738, 394, 404, 2998,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8758, 0, 3, 2518, 7774, 404, 414, 3028,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8818, 0, 3, 2536, 7810, 414, 424, 3058,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8878, 0, 3, 2554, 7846, 424, 434, 3088,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8938, 0, 3, 2572, 7882, 434, 444, 3118,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8998, 0, 3, 7270, 7306, 2638, 8038, 464,
                                                 479, 3193, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9088, 0, 3, 7306, 7342, 2668, 8098, 479,
                                                 494, 3238, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9178, 0, 3, 7342, 7378, 2698, 8158, 494,
                                                 509, 3283, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9268, 0, 3, 7378, 7414, 2728, 8218, 509,
                                                 524, 3328, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9358, 0, 3, 7414, 7450, 2758, 8278, 524,
                                                 539, 3373, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9448, 0, 3, 7450, 7486, 2788, 8338, 539,
                                                 554, 3418, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9538, 0, 3, 7486, 7522, 2818, 8398, 554,
                                                 569, 3463, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9628, 0, 3, 7594, 7630, 2908, 8578, 599,
                                                 614, 3553, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9718, 0, 3, 7630, 7666, 2938, 8638, 614,
                                                 629, 3598, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9808, 0, 3, 7666, 7702, 2968, 8698, 629,
                                                 644, 3643, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9898, 0, 3, 7702, 7738, 2998, 8758, 644,
                                                 659, 3688, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9988, 0, 3, 7738, 7774, 3028, 8818, 659,
                                                 674, 3733, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10078, 0, 3, 7774, 7810, 3058, 8878,
                                                 674, 689, 3778, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10168, 0, 3, 7810, 7846, 3088, 8938,
                                                 689, 704, 3823, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10258, 0, 3, 7918, 7978, 3148, 8998,
                                                 734, 755, 3868, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10384, 0, 3, 7978, 8038, 3193, 9088,
                                                 755, 776, 3931, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10510, 0, 3, 8038, 8098, 3238, 9178,
                                                 776, 797, 3994, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10636, 0, 3, 8098, 8158, 3283, 9268,
                                                 797, 818, 4057, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10762, 0, 3, 8158, 8218, 3328, 9358,
                                                 818, 839, 4120, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10888, 0, 3, 8218, 8278, 3373, 9448,
                                                 839, 860, 4183, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11014, 0, 3, 8278, 8338, 3418, 9538,
                                                 860, 881, 4246, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11140, 0, 3, 8458, 8518, 3508, 9628,
                                                 923, 944, 4309, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11266, 0, 3, 8518, 8578, 3553, 9718,
                                                 944, 965, 4372, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11392, 0, 3, 8578, 8638, 3598, 9808,
                                                 965, 986, 4435, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11518, 0, 3, 8638, 8698, 3643, 9898,
                                                 986, 1007, 4498, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11644, 0, 3, 8698, 8758, 3688, 9988,
                                                 1007, 1028, 4561, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11770, 0, 3, 8758, 8818, 3733, 10078,
                                                 1028, 1049, 4624, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11896, 0, 3, 8818, 8878, 3778, 10168,
                                                 1049, 1070, 4687, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12022, 0, 3, 8998, 9088, 3931, 10510,
                                                 1112, 1140, 4834, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12190, 0, 3, 9088, 9178, 3994, 10636,
                                                 1140, 1168, 4918, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12358, 0, 3, 9178, 9268, 4057, 10762,
                                                 1168, 1196, 5002, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12526, 0, 3, 9268, 9358, 4120, 10888,
                                                 1196, 1224, 5086, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12694, 0, 3, 9358, 9448, 4183, 11014,
                                                 1224, 1252, 5170, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12862, 0, 3, 9628, 9718, 4372, 11392,
                                                 1308, 1336, 5338, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13030, 0, 3, 9718, 9808, 4435, 11518,
                                                 1336, 1364, 5422, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13198, 0, 3, 9808, 9898, 4498, 11644,
                                                 1364, 1392, 5506, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13366, 0, 3, 9898, 9988, 4561, 11770,
                                                 1392, 1420, 5590, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13534, 0, 3, 9988, 10078, 4624, 11896,
                                                 1420, 1448, 5674, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13702, 0, 3, 10258, 10384, 4750, 12022,
                                                 1504, 1540, 5758, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 13918, 0, 3, 10384, 10510, 4834, 12190,
                                                 1540, 1576, 5866, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14134, 0, 3, 10510, 10636, 4918, 12358,
                                                 1576, 1612, 5974, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14350, 0, 3, 10636, 10762, 5002, 12526,
                                                 1612, 1648, 6082, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14566, 0, 3, 10762, 10888, 5086, 12694,
                                                 1648, 1684, 6190, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14782, 0, 3, 11140, 11266, 5254, 12862,
                                                 1756, 1792, 6298, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 14998, 0, 3, 11266, 11392, 5338, 13030,
                                                 1792, 1828, 6406, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15214, 0, 3, 11392, 11518, 5422, 13198,
                                                 1828, 1864, 6514, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15430, 0, 3, 11518, 11644, 5506, 13366,
                                                 1864, 1900, 6622, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15646, 0, 3, 11644, 11770, 5590, 13534,
                                                 1900, 1936, 6730, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15862, 3, 2008, 2011, 6844, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15872, 3, 2011, 2014, 6850, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15882, 3, 2014, 2017, 6856, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15892, 3, 2017, 2020, 6862, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15902, 3, 2020, 2023, 6868, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15912, 3, 2023, 2026, 6874, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15922, 3, 2026, 2029, 6880, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15932, 3, 2029, 2032, 6886, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15942, 3, 2038, 2041, 6898, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15952, 3, 2041, 2044, 6904, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15962, 3, 2044, 2047, 6910, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15972, 3, 2047, 2050, 6916, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15982, 3, 2050, 2053, 6922, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 15992, 3, 2053, 2056, 6928, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16002, 3, 2056, 2059, 6934, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 16012, 3, 2059, 2062, 6940, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 16022, 0, 3, 6838, 15862, 6964, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16052, 0, 3, 6844, 15872, 6982, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16082, 0, 3, 6850, 15882, 7000, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16112, 0, 3, 6856, 15892, 7018, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16142, 0, 3, 6862, 15902, 7036, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16172, 0, 3, 6868, 15912, 7054, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16202, 0, 3, 6874, 15922, 7072, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16232, 0, 3, 6880, 15932, 7090, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16262, 0, 3, 6892, 15942, 7126, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16292, 0, 3, 6898, 15952, 7144, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16322, 0, 3, 6904, 15962, 7162, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16352, 0, 3, 6910, 15972, 7180, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16382, 0, 3, 6916, 15982, 7198, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16412, 0, 3, 6922, 15992, 7216, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16442, 0, 3, 6928, 16002, 7234, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 16472, 0, 3, 6934, 16012, 7252, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 16502, 0, 3, 6946, 16022, 2248, 2266,
                                                 7306, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16562, 0, 3, 6964, 16052, 2266, 2284,
                                                 7342, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16622, 0, 3, 6982, 16082, 2284, 2302,
                                                 7378, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16682, 0, 3, 7000, 16112, 2302, 2320,
                                                 7414, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16742, 0, 3, 7018, 16142, 2320, 2338,
                                                 7450, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16802, 0, 3, 7036, 16172, 2338, 2356,
                                                 7486, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16862, 0, 3, 7054, 16202, 2356, 2374,
                                                 7522, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16922, 0, 3, 7072, 16232, 2374, 2392,
                                                 7558, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 16982, 0, 3, 7108, 16262, 2428, 2446,
                                                 7630, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17042, 0, 3, 7126, 16292, 2446, 2464,
                                                 7666, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17102, 0, 3, 7144, 16322, 2464, 2482,
                                                 7702, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17162, 0, 3, 7162, 16352, 2482, 2500,
                                                 7738, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17222, 0, 3, 7180, 16382, 2500, 2518,
                                                 7774, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17282, 0, 3, 7198, 16412, 2518, 2536,
                                                 7810, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17342, 0, 3, 7216, 16442, 2536, 2554,
                                                 7846, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 17402, 0, 3, 7234, 16472, 2554, 2572,
                                                 7882, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17462, 0, 3, 7306, 16562, 2608, 2638,
                                                 8038, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17562, 0, 3, 7342, 16622, 2638, 2668,
                                                 8098, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17662, 0, 3, 7378, 16682, 2668, 2698,
                                                 8158, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17762, 0, 3, 7414, 16742, 2698, 2728,
                                                 8218, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17862, 0, 3, 7450, 16802, 2728, 2758,
                                                 8278, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17962, 0, 3, 7486, 16862, 2758, 2788,
                                                 8338, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18062, 0, 3, 7522, 16922, 2788, 2818,
                                                 8398, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18162, 0, 3, 7630, 17042, 2878, 2908,
                                                 8578, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18262, 0, 3, 7666, 17102, 2908, 2938,
                                                 8638, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18362, 0, 3, 7702, 17162, 2938, 2968,
                                                 8698, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18462, 0, 3, 7738, 17222, 2968, 2998,
                                                 8758, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18562, 0, 3, 7774, 17282, 2998, 3028,
                                                 8818, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18662, 0, 3, 7810, 17342, 3028, 3058,
                                                 8878, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 18762, 0, 3, 7846, 17402, 3058, 3088,
                                                 8938, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18862, 0, 3, 16502, 16562, 8038, 17562,
                                                 3148, 3193, 9088, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19012, 0, 3, 16562, 16622, 8098, 17662,
                                                 3193, 3238, 9178, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19162, 0, 3, 16622, 16682, 8158, 17762,
                                                 3238, 3283, 9268, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19312, 0, 3, 16682, 16742, 8218, 17862,
                                                 3283, 3328, 9358, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19462, 0, 3, 16742, 16802, 8278, 17962,
                                                 3328, 3373, 9448, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19612, 0, 3, 16802, 16862, 8338, 18062,
                                                 3373, 3418, 9538, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19762, 0, 3, 16982, 17042, 8578, 18262,
                                                 3508, 3553, 9718, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19912, 0, 3, 17042, 17102, 8638, 18362,
                                                 3553, 3598, 9808, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20062, 0, 3, 17102, 17162, 8698, 18462,
                                                 3598, 3643, 9898, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20212, 0, 3, 17162, 17222, 8758, 18562,
                                                 3643, 3688, 9988, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20362, 0, 3, 17222, 17282, 8818, 18662,
                                                 3688, 3733, 10078, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 20512, 0, 3, 17282, 17342, 8878, 18762,
                                                 3733, 3778, 10168, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20662, 0, 3, 17462, 17562, 9088, 19012,
                                                 3868, 3931, 10510, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20872, 0, 3, 17562, 17662, 9178, 19162,
                                                 3931, 3994, 10636, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21082, 0, 3, 17662, 17762, 9268, 19312,
                                                 3994, 4057, 10762, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21292, 0, 3, 17762, 17862, 9358, 19462,
                                                 4057, 4120, 10888, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21502, 0, 3, 17862, 17962, 9448, 19612,
                                                 4120, 4183, 11014, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21712, 0, 3, 18162, 18262, 9718, 19912,
                                                 4309, 4372, 11392, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21922, 0, 3, 18262, 18362, 9808, 20062,
                                                 4372, 4435, 11518, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22132, 0, 3, 18362, 18462, 9898, 20212,
                                                 4435, 4498, 11644, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22342, 0, 3, 18462, 18562, 9988, 20362,
                                                 4498, 4561, 11770, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22552, 0, 3, 18562, 18662, 10078, 20512,
                                                 4561, 4624, 11896, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22762, 0, 3, 18862, 19012, 10510, 20872,
                                                 4750, 4834, 12190, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23042, 0, 3, 19012, 19162, 10636, 21082,
                                                 4834, 4918, 12358, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23322, 0, 3, 19162, 19312, 10762, 21292,
                                                 4918, 5002, 12526, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23602, 0, 3, 19312, 19462, 10888, 21502,
                                                 5002, 5086, 12694, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23882, 0, 3, 19762, 19912, 11392, 21922,
                                                 5254, 5338, 13030, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24162, 0, 3, 19912, 20062, 11518, 22132,
                                                 5338, 5422, 13198, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24442, 0, 3, 20062, 20212, 11644, 22342,
                                                 5422, 5506, 13366, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24722, 0, 3, 20212, 20362, 11770, 22552,
                                                 5506, 5590, 13534, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25002, 0, 3, 20662, 20872, 12190, 23042,
                                                 5758, 5866, 14134, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25362, 0, 3, 20872, 21082, 12358, 23322,
                                                 5866, 5974, 14350, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 25722, 0, 3, 21082, 21292, 12526, 23602,
                                                 5974, 6082, 14566, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 26082, 0, 3, 21712, 21922, 13030, 24162,
                                                 6298, 6406, 15214, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 26442, 0, 3, 21922, 22132, 13198, 24442,
                                                 6406, 6514, 15430, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 26802, 0, 3, 22132, 22342, 13366, 24722,
                                                 6514, 6622, 15646, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27162, 3, 6838, 6844, 15872, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27177, 3, 6844, 6850, 15882, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27192, 3, 6850, 6856, 15892, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27207, 3, 6856, 6862, 15902, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27222, 3, 6862, 6868, 15912, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27237, 3, 6868, 6874, 15922, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27252, 3, 6874, 6880, 15932, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27267, 3, 6892, 6898, 15952, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27282, 3, 6898, 6904, 15962, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27297, 3, 6904, 6910, 15972, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27312, 3, 6910, 6916, 15982, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27327, 3, 6916, 6922, 15992, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27342, 3, 6922, 6928, 16002, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 27357, 3, 6928, 6934, 16012, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 27372, 0, 3, 15862, 27162, 16052, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27417, 0, 3, 15872, 27177, 16082, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27462, 0, 3, 15882, 27192, 16112, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27507, 0, 3, 15892, 27207, 16142, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27552, 0, 3, 15902, 27222, 16172, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27597, 0, 3, 15912, 27237, 16202, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27642, 0, 3, 15922, 27252, 16232, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27687, 0, 3, 15942, 27267, 16292, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27732, 0, 3, 15952, 27282, 16322, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27777, 0, 3, 15962, 27297, 16352, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27822, 0, 3, 15972, 27312, 16382, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27867, 0, 3, 15982, 27327, 16412, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27912, 0, 3, 15992, 27342, 16442, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 27957, 0, 3, 16002, 27357, 16472, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 28002, 0, 3, 16022, 27372, 7270, 7306,
                                                 16562, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28092, 0, 3, 16052, 27417, 7306, 7342,
                                                 16622, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28182, 0, 3, 16082, 27462, 7342, 7378,
                                                 16682, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28272, 0, 3, 16112, 27507, 7378, 7414,
                                                 16742, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28362, 0, 3, 16142, 27552, 7414, 7450,
                                                 16802, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28452, 0, 3, 16172, 27597, 7450, 7486,
                                                 16862, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28542, 0, 3, 16202, 27642, 7486, 7522,
                                                 16922, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28632, 0, 3, 16262, 27687, 7594, 7630,
                                                 17042, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28722, 0, 3, 16292, 27732, 7630, 7666,
                                                 17102, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28812, 0, 3, 16322, 27777, 7666, 7702,
                                                 17162, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28902, 0, 3, 16352, 27822, 7702, 7738,
                                                 17222, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 28992, 0, 3, 16382, 27867, 7738, 7774,
                                                 17282, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 29082, 0, 3, 16412, 27912, 7774, 7810,
                                                 17342, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 29172, 0, 3, 16442, 27957, 7810, 7846,
                                                 17402, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29262, 0, 3, 16502, 28002, 7918, 7978,
                                                 17462, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29412, 0, 3, 16562, 28092, 7978, 8038,
                                                 17562, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29562, 0, 3, 16622, 28182, 8038, 8098,
                                                 17662, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29712, 0, 3, 16682, 28272, 8098, 8158,
                                                 17762, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29862, 0, 3, 16742, 28362, 8158, 8218,
                                                 17862, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 30012, 0, 3, 16802, 28452, 8218, 8278,
                                                 17962, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 30162, 0, 3, 16862, 28542, 8278, 8338,
                                                 18062, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 30312, 0, 3, 16982, 28632, 8458, 8518,
                                                 18162, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 30462, 0, 3, 17042, 28722, 8518, 8578,
                                                 18262, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 30612, 0, 3, 17102, 28812, 8578, 8638,
                                                 18362, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 30762, 0, 3, 17162, 28902, 8638, 8698,
                                                 18462, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 30912, 0, 3, 17222, 28992, 8698, 8758,
                                                 18562, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 31062, 0, 3, 17282, 29082, 8758, 8818,
                                                 18662, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 31212, 0, 3, 17342, 29172, 8818, 8878,
                                                 18762, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31362, 0, 3, 28002, 28092, 17562, 29562,
                                                 8998, 9088, 19012, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31587, 0, 3, 28092, 28182, 17662, 29712,
                                                 9088, 9178, 19162, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31812, 0, 3, 28182, 28272, 17762, 29862,
                                                 9178, 9268, 19312, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32037, 0, 3, 28272, 28362, 17862, 30012,
                                                 9268, 9358, 19462, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32262, 0, 3, 28362, 28452, 17962, 30162,
                                                 9358, 9448, 19612, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32487, 0, 3, 28632, 28722, 18262, 30612,
                                                 9628, 9718, 19912, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32712, 0, 3, 28722, 28812, 18362, 30762,
                                                 9718, 9808, 20062, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32937, 0, 3, 28812, 28902, 18462, 30912,
                                                 9808, 9898, 20212, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 33162, 0, 3, 28902, 28992, 18562, 31062,
                                                 9898, 9988, 20362, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 33387, 0, 3, 28992, 29082, 18662, 31212,
                                                 9988, 10078, 20512, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33612, 0, 3, 29262, 29412, 18862, 31362,
                                                 10258, 10384, 20662, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33927, 0, 3, 29412, 29562, 19012, 31587,
                                                 10384, 10510, 20872, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 34242, 0, 3, 29562, 29712, 19162, 31812,
                                                 10510, 10636, 21082, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 34557, 0, 3, 29712, 29862, 19312, 32037,
                                                 10636, 10762, 21292, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 34872, 0, 3, 29862, 30012, 19462, 32262,
                                                 10762, 10888, 21502, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 35187, 0, 3, 30312, 30462, 19762, 32487,
                                                 11140, 11266, 21712, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 35502, 0, 3, 30462, 30612, 19912, 32712,
                                                 11266, 11392, 21922, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 35817, 0, 3, 30612, 30762, 20062, 32937,
                                                 11392, 11518, 22132, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 36132, 0, 3, 30762, 30912, 20212, 33162,
                                                 11518, 11644, 22342, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 36447, 0, 3, 30912, 31062, 20362, 33387,
                                                 11644, 11770, 22552, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 36762, 0, 3, 31362, 31587, 20872, 34242,
                                                 12022, 12190, 23042, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 37182, 0, 3, 31587, 31812, 21082, 34557,
                                                 12190, 12358, 23322, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 37602, 0, 3, 31812, 32037, 21292, 34872,
                                                 12358, 12526, 23602, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 38022, 0, 3, 32487, 32712, 21922, 35817,
                                                 12862, 13030, 24162, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 38442, 0, 3, 32712, 32937, 22132, 36132,
                                                 13030, 13198, 24442, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 38862, 0, 3, 32937, 33162, 22342, 36447,
                                                 13198, 13366, 24722, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 39282, 0, 3, 33612, 33927, 22762, 36762,
                                                 13702, 13918, 25002, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 39822, 0, 3, 33927, 34242, 23042, 37182,
                                                 13918, 14134, 25362, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 40362, 0, 3, 34242, 34557, 23322, 37602,
                                                 14134, 14350, 25722, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 40902, 0, 3, 35187, 35502, 23882, 38022,
                                                 14782, 14998, 26082, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 41442, 0, 3, 35502, 35817, 24162, 38442,
                                                 14998, 15214, 26442, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 41982, 0, 3, 35817, 36132, 24442, 38862,
                                                 15214, 15430, 26802, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42522, 3, 15862, 15872, 27177, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42543, 3, 15872, 15882, 27192, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42564, 3, 15882, 15892, 27207, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42585, 3, 15892, 15902, 27222, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42606, 3, 15902, 15912, 27237, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42627, 3, 15912, 15922, 27252, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42648, 3, 15942, 15952, 27282, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42669, 3, 15952, 15962, 27297, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42690, 3, 15962, 15972, 27312, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42711, 3, 15972, 15982, 27327, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42732, 3, 15982, 15992, 27342, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 42753, 3, 15992, 16002, 27357, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 42774, 0, 3, 27162, 42522, 27417, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 42837, 0, 3, 27177, 42543, 27462, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 42900, 0, 3, 27192, 42564, 27507, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 42963, 0, 3, 27207, 42585, 27552, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43026, 0, 3, 27222, 42606, 27597, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43089, 0, 3, 27237, 42627, 27642, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43152, 0, 3, 27267, 42648, 27732, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43215, 0, 3, 27282, 42669, 27777, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43278, 0, 3, 27297, 42690, 27822, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43341, 0, 3, 27312, 42711, 27867, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43404, 0, 3, 27327, 42732, 27912, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 43467, 0, 3, 27342, 42753, 27957, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 43530, 0, 3, 27372, 42774, 16502, 16562,
                                                 28092, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43656, 0, 3, 27417, 42837, 16562, 16622,
                                                 28182, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43782, 0, 3, 27462, 42900, 16622, 16682,
                                                 28272, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43908, 0, 3, 27507, 42963, 16682, 16742,
                                                 28362, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44034, 0, 3, 27552, 43026, 16742, 16802,
                                                 28452, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44160, 0, 3, 27597, 43089, 16802, 16862,
                                                 28542, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44286, 0, 3, 27687, 43152, 16982, 17042,
                                                 28722, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44412, 0, 3, 27732, 43215, 17042, 17102,
                                                 28812, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44538, 0, 3, 27777, 43278, 17102, 17162,
                                                 28902, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44664, 0, 3, 27822, 43341, 17162, 17222,
                                                 28992, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44790, 0, 3, 27867, 43404, 17222, 17282,
                                                 29082, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 44916, 0, 3, 27912, 43467, 17282, 17342,
                                                 29172, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45042, 0, 3, 28092, 43656, 17462, 17562,
                                                 29562, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45252, 0, 3, 28182, 43782, 17562, 17662,
                                                 29712, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45462, 0, 3, 28272, 43908, 17662, 17762,
                                                 29862, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45672, 0, 3, 28362, 44034, 17762, 17862,
                                                 30012, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45882, 0, 3, 28452, 44160, 17862, 17962,
                                                 30162, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 46092, 0, 3, 28722, 44412, 18162, 18262,
                                                 30612, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 46302, 0, 3, 28812, 44538, 18262, 18362,
                                                 30762, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 46512, 0, 3, 28902, 44664, 18362, 18462,
                                                 30912, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 46722, 0, 3, 28992, 44790, 18462, 18562,
                                                 31062, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 46932, 0, 3, 29082, 44916, 18562, 18662,
                                                 31212, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47142, 0, 3, 43530, 43656, 29562, 45252,
                                                 18862, 19012, 31587, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47457, 0, 3, 43656, 43782, 29712, 45462,
                                                 19012, 19162, 31812, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47772, 0, 3, 43782, 43908, 29862, 45672,
                                                 19162, 19312, 32037, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 48087, 0, 3, 43908, 44034, 30012, 45882,
                                                 19312, 19462, 32262, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 48402, 0, 3, 44286, 44412, 30612, 46302,
                                                 19762, 19912, 32712, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 48717, 0, 3, 44412, 44538, 30762, 46512,
                                                 19912, 20062, 32937, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 49032, 0, 3, 44538, 44664, 30912, 46722,
                                                 20062, 20212, 33162, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 49347, 0, 3, 44664, 44790, 31062, 46932,
                                                 20212, 20362, 33387, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 49662, 0, 3, 45042, 45252, 31587, 47457,
                                                 20662, 20872, 34242, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 50103, 0, 3, 45252, 45462, 31812, 47772,
                                                 20872, 21082, 34557, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 50544, 0, 3, 45462, 45672, 32037, 48087,
                                                 21082, 21292, 34872, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 50985, 0, 3, 46092, 46302, 32712, 48717,
                                                 21712, 21922, 35817, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 51426, 0, 3, 46302, 46512, 32937, 49032,
                                                 21922, 22132, 36132, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 51867, 0, 3, 46512, 46722, 33162, 49347,
                                                 22132, 22342, 36447, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 52308, 0, 3, 47142, 47457, 34242, 50103,
                                                 22762, 23042, 37182, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 52896, 0, 3, 47457, 47772, 34557, 50544,
                                                 23042, 23322, 37602, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 53484, 0, 3, 48402, 48717, 35817, 51426,
                                                 23882, 24162, 38442, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 54072, 0, 3, 48717, 49032, 36132, 51867,
                                                 24162, 24442, 38862, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 54660, 0, 3, 49662, 50103, 37182, 52896,
                                                 25002, 25362, 40362, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 55416, 0, 3, 50985, 51426, 38442, 54072,
                                                 26082, 26442, 41982, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56172, 3, 27162, 27177, 42543, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56200, 3, 27177, 27192, 42564, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56228, 3, 27192, 27207, 42585, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56256, 3, 27207, 27222, 42606, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56284, 3, 27222, 27237, 42627, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56312, 3, 27267, 27282, 42669, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56340, 3, 27282, 27297, 42690, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56368, 3, 27297, 27312, 42711, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56396, 3, 27312, 27327, 42732, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56424, 3, 27327, 27342, 42753, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 56452, 0, 3, 42522, 56172, 42837, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56536, 0, 3, 42543, 56200, 42900, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56620, 0, 3, 42564, 56228, 42963, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56704, 0, 3, 42585, 56256, 43026, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56788, 0, 3, 42606, 56284, 43089, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56872, 0, 3, 42648, 56312, 43215, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56956, 0, 3, 42669, 56340, 43278, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57040, 0, 3, 42690, 56368, 43341, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57124, 0, 3, 42711, 56396, 43404, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57208, 0, 3, 42732, 56424, 43467, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 57292, 0, 3, 42774, 56452, 28002, 28092,
                                                 43656, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 57460, 0, 3, 42837, 56536, 28092, 28182,
                                                 43782, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 57628, 0, 3, 42900, 56620, 28182, 28272,
                                                 43908, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 57796, 0, 3, 42963, 56704, 28272, 28362,
                                                 44034, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 57964, 0, 3, 43026, 56788, 28362, 28452,
                                                 44160, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58132, 0, 3, 43152, 56872, 28632, 28722,
                                                 44412, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58300, 0, 3, 43215, 56956, 28722, 28812,
                                                 44538, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58468, 0, 3, 43278, 57040, 28812, 28902,
                                                 44664, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58636, 0, 3, 43341, 57124, 28902, 28992,
                                                 44790, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58804, 0, 3, 43404, 57208, 28992, 29082,
                                                 44916, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 58972, 0, 3, 43530, 57292, 29262, 29412,
                                                 45042, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 59252, 0, 3, 43656, 57460, 29412, 29562,
                                                 45252, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 59532, 0, 3, 43782, 57628, 29562, 29712,
                                                 45462, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 59812, 0, 3, 43908, 57796, 29712, 29862,
                                                 45672, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60092, 0, 3, 44034, 57964, 29862, 30012,
                                                 45882, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60372, 0, 3, 44286, 58132, 30312, 30462,
                                                 46092, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60652, 0, 3, 44412, 58300, 30462, 30612,
                                                 46302, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60932, 0, 3, 44538, 58468, 30612, 30762,
                                                 46512, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 61212, 0, 3, 44664, 58636, 30762, 30912,
                                                 46722, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 61492, 0, 3, 44790, 58804, 30912, 31062,
                                                 46932, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 61772, 0, 3, 57292, 57460, 45252, 59532,
                                                 31362, 31587, 47457, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 62192, 0, 3, 57460, 57628, 45462, 59812,
                                                 31587, 31812, 47772, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 62612, 0, 3, 57628, 57796, 45672, 60092,
                                                 31812, 32037, 48087, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 63032, 0, 3, 58132, 58300, 46302, 60932,
                                                 32487, 32712, 48717, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 63452, 0, 3, 58300, 58468, 46512, 61212,
                                                 32712, 32937, 49032, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 63872, 0, 3, 58468, 58636, 46722, 61492,
                                                 32937, 33162, 49347, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 64292, 0, 3, 58972, 59252, 47142, 61772,
                                                 33612, 33927, 49662, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 64880, 0, 3, 59252, 59532, 47457, 62192,
                                                 33927, 34242, 50103, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 65468, 0, 3, 59532, 59812, 47772, 62612,
                                                 34242, 34557, 50544, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 66056, 0, 3, 60372, 60652, 48402, 63032,
                                                 35187, 35502, 50985, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 66644, 0, 3, 60652, 60932, 48717, 63452,
                                                 35502, 35817, 51426, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 67232, 0, 3, 60932, 61212, 49032, 63872,
                                                 35817, 36132, 51867, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 67820, 0, 3, 61772, 62192, 50103, 65468,
                                                 36762, 37182, 52896, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 68604, 0, 3, 63032, 63452, 51426, 67232,
                                                 38022, 38442, 54072, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 69388, 0, 3, 64292, 64880, 52308, 67820,
                                                 39282, 39822, 54660, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 70396, 0, 3, 66056, 66644, 53484, 68604,
                                                 40902, 41442, 55416, ncols, alpha, beta, p);

            simdgeo::geom_i_x(buffer, 71404, 66056, 70396, 1, 28, ncols, alpha);

            simdgeo::geom_i_y(buffer, 72188, 66056, 70396, 1, 28, ncols, alpha);

            simdgeo::geom_i_z(buffer, 72972, 66056, 70396, 1, 28, ncols, alpha);

            simdgeo::geom_i_x(buffer, 73756, 64292, 69388, 1, 28, ncols, alpha);

            simdgeo::geom_i_y(buffer, 74540, 64292, 69388, 1, 28, ncols, alpha);

            simdgeo::geom_i_z(buffer, 75324, 64292, 69388, 1, 28, ncols, alpha);

            simdfunc::contract_primitives(buffer, 76108, 73756, 2352, ncols);

            simdfunc::contract_primitives(buffer, 78460, 71404, 2352, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 80812, 78460, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 80812, 13, nmax);

    simdtrf::transform_i_inner(buffer, 80812, 79244, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 169 * nvalues, nvalues, buffer, 80812, 13, nmax);

    simdtrf::transform_i_inner(buffer, 80812, 80028, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 338 * nvalues, nvalues, buffer, 80812, 13, nmax);

    simdtrf::transform_i_inner(buffer, 80812, 76108, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 507 * nvalues, nvalues, buffer, 80812, 13, nmax);

    simdtrf::transform_i_inner(buffer, 80812, 76892, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 676 * nvalues, nvalues, buffer, 80812, 13, nmax);

    simdtrf::transform_i_inner(buffer, 80812, 77676, 28, 1, nmax);

    simdtrf::transform_i_outer(values + 845 * nvalues, nvalues, buffer, 80812, 13, nmax);
}

}  // namespace simdt2ceri
