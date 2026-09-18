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


#include "SimdElectronRepulsionGeom10RsRecKH.hpp"

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
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGG.hpp"
#include "SimdElectronRepulsionVrrRecGH.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecHD.hpp"
#include "SimdElectronRepulsionVrrRecHF.hpp"
#include "SimdElectronRepulsionVrrRecHG.hpp"
#include "SimdElectronRepulsionVrrRecHH.hpp"
#include "SimdElectronRepulsionVrrRecHP.hpp"
#include "SimdElectronRepulsionVrrRecHS.hpp"
#include "SimdElectronRepulsionVrrRecID.hpp"
#include "SimdElectronRepulsionVrrRecIF.hpp"
#include "SimdElectronRepulsionVrrRecIG.hpp"
#include "SimdElectronRepulsionVrrRecIH.hpp"
#include "SimdElectronRepulsionVrrRecIP.hpp"
#include "SimdElectronRepulsionVrrRecIS.hpp"
#include "SimdElectronRepulsionVrrRecKD.hpp"
#include "SimdElectronRepulsionVrrRecKF.hpp"
#include "SimdElectronRepulsionVrrRecKG.hpp"
#include "SimdElectronRepulsionVrrRecKH.hpp"
#include "SimdElectronRepulsionVrrRecKP.hpp"
#include "SimdElectronRepulsionVrrRecKS.hpp"
#include "SimdElectronRepulsionVrrRecLD.hpp"
#include "SimdElectronRepulsionVrrRecLF.hpp"
#include "SimdElectronRepulsionVrrRecLG.hpp"
#include "SimdElectronRepulsionVrrRecLH.hpp"
#include "SimdElectronRepulsionVrrRecLP.hpp"
#include "SimdElectronRepulsionVrrRecLS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryK1.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_kh_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_kh_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 75818, 70886, 4536, nvalues);

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

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 22, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 82, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 85, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 88, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 91, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 94, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 97, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 100, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 103, 0, 33, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 7, 8, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 8, 9, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 9, 10, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 10, 11, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 11, 12, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 12, 13, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 13, 14, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 14, 15, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 15, 16, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 16, 17, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 17, 18, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 21, 22, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 178, 0, 22, 23, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 184, 0, 23, 24, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 190, 0, 24, 25, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 196, 0, 25, 26, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 202, 0, 26, 27, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 208, 0, 27, 28, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 214, 0, 28, 29, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 220, 0, 29, 30, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 226, 0, 30, 31, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 232, 0, 31, 32, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 238, 0, 34, 37, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 248, 0, 37, 40, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 258, 0, 40, 43, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 268, 0, 43, 46, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 278, 0, 46, 49, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 288, 0, 49, 52, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 298, 0, 52, 55, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 308, 0, 55, 58, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 318, 0, 58, 61, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 328, 0, 61, 64, 166, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 338, 0, 70, 73, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 348, 0, 73, 76, 184, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 358, 0, 76, 79, 190, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 368, 0, 79, 82, 196, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 378, 0, 82, 85, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 388, 0, 85, 88, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 398, 0, 88, 91, 214, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 408, 0, 91, 94, 220, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 418, 0, 94, 97, 226, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 428, 0, 97, 100, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 438, 0, 106, 112, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 453, 0, 112, 118, 258, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 468, 0, 118, 124, 268, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 483, 0, 124, 130, 278, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 498, 0, 130, 136, 288, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 513, 0, 136, 142, 298, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 528, 0, 142, 148, 308, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 543, 0, 148, 154, 318, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 558, 0, 154, 160, 328, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 573, 0, 172, 178, 348, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 588, 0, 178, 184, 358, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 603, 0, 184, 190, 368, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 618, 0, 190, 196, 378, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 633, 0, 196, 202, 388, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 648, 0, 202, 208, 398, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 663, 0, 208, 214, 408, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 678, 0, 214, 220, 418, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 693, 0, 220, 226, 428, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 708, 0, 238, 248, 453, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 729, 0, 248, 258, 468, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 750, 0, 258, 268, 483, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 771, 0, 268, 278, 498, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 792, 0, 278, 288, 513, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 813, 0, 288, 298, 528, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 834, 0, 298, 308, 543, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 855, 0, 308, 318, 558, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 876, 0, 338, 348, 588, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 897, 0, 348, 358, 603, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 918, 0, 358, 368, 618, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 939, 0, 368, 378, 633, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 960, 0, 378, 388, 648, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 981, 0, 388, 398, 663, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1002, 0, 398, 408, 678, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1023, 0, 408, 418, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1044, 0, 438, 453, 729, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1072, 0, 453, 468, 750, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1100, 0, 468, 483, 771, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1128, 0, 483, 498, 792, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1156, 0, 498, 513, 813, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1184, 0, 513, 528, 834, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1212, 0, 528, 543, 855, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1240, 0, 573, 588, 897, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1268, 0, 588, 603, 918, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1296, 0, 603, 618, 939, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1324, 0, 618, 633, 960, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1352, 0, 633, 648, 981, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1380, 0, 648, 663, 1002, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1408, 0, 663, 678, 1023, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1436, 0, 708, 729, 1072, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1472, 0, 729, 750, 1100, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1508, 0, 750, 771, 1128, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1544, 0, 771, 792, 1156, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1580, 0, 792, 813, 1184, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1616, 0, 813, 834, 1212, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1652, 0, 876, 897, 1268, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1688, 0, 897, 918, 1296, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1724, 0, 918, 939, 1324, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1760, 0, 939, 960, 1352, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1796, 0, 960, 981, 1380, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1832, 0, 981, 1002, 1408, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1868, 0, 1044, 1072, 1472, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1913, 0, 1072, 1100, 1508, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1958, 0, 1100, 1128, 1544, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2003, 0, 1128, 1156, 1580, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2048, 0, 1156, 1184, 1616, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2093, 0, 1240, 1268, 1688, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2138, 0, 1268, 1296, 1724, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2183, 0, 1296, 1324, 1760, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2228, 0, 1324, 1352, 1796, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2273, 0, 1352, 1380, 1832, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 2318, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2321, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2324, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2327, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2330, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2333, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2336, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2339, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2342, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2345, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2348, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2351, 3, 23, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2354, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2357, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2360, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2363, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2366, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2369, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2372, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2375, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2378, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2381, 3, 33, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2384, 3, 8, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2393, 3, 9, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2402, 3, 10, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2411, 3, 11, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2420, 3, 12, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2429, 3, 13, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2438, 3, 14, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2447, 3, 15, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2456, 3, 16, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2465, 3, 17, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2474, 3, 18, 67, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2483, 3, 22, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2492, 3, 23, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2501, 3, 24, 79, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2510, 3, 25, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2519, 3, 26, 85, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2528, 3, 27, 88, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2537, 3, 28, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2546, 3, 29, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2555, 3, 30, 97, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2564, 3, 31, 100, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2573, 3, 32, 103, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2582, 0, 3, 34, 2384, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2600, 0, 3, 37, 2393, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2618, 0, 3, 40, 2402, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2636, 0, 3, 43, 2411, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2654, 0, 3, 46, 2420, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2672, 0, 3, 49, 2429, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2690, 0, 3, 52, 2438, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2708, 0, 3, 55, 2447, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2726, 0, 3, 58, 2456, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2744, 0, 3, 61, 2465, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2762, 0, 3, 64, 2474, 166, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2780, 0, 3, 70, 2483, 172, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2798, 0, 3, 73, 2492, 178, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2816, 0, 3, 76, 2501, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2834, 0, 3, 79, 2510, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2852, 0, 3, 82, 2519, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2870, 0, 3, 85, 2528, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2888, 0, 3, 88, 2537, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2906, 0, 3, 91, 2546, 214, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2924, 0, 3, 94, 2555, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2942, 0, 3, 97, 2564, 226, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2960, 0, 3, 100, 2573, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2978, 0, 3, 112, 2618, 248, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3008, 0, 3, 118, 2636, 258, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3038, 0, 3, 124, 2654, 268, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3068, 0, 3, 130, 2672, 278, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3098, 0, 3, 136, 2690, 288, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3128, 0, 3, 142, 2708, 298, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3158, 0, 3, 148, 2726, 308, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3188, 0, 3, 154, 2744, 318, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3218, 0, 3, 160, 2762, 328, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3248, 0, 3, 178, 2816, 348, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3278, 0, 3, 184, 2834, 358, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3308, 0, 3, 190, 2852, 368, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3338, 0, 3, 196, 2870, 378, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3368, 0, 3, 202, 2888, 388, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3398, 0, 3, 208, 2906, 398, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3428, 0, 3, 214, 2924, 408, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3458, 0, 3, 220, 2942, 418, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3488, 0, 3, 226, 2960, 428, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3518, 0, 3, 238, 2978, 438, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3563, 0, 3, 248, 3008, 453, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3608, 0, 3, 258, 3038, 468, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3653, 0, 3, 268, 3068, 483, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3698, 0, 3, 278, 3098, 498, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3743, 0, 3, 288, 3128, 513, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3788, 0, 3, 298, 3158, 528, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3833, 0, 3, 308, 3188, 543, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3878, 0, 3, 318, 3218, 558, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3923, 0, 3, 338, 3248, 573, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3968, 0, 3, 348, 3278, 588, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4013, 0, 3, 358, 3308, 603, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4058, 0, 3, 368, 3338, 618, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4103, 0, 3, 378, 3368, 633, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4148, 0, 3, 388, 3398, 648, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4193, 0, 3, 398, 3428, 663, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4238, 0, 3, 408, 3458, 678, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4283, 0, 3, 418, 3488, 693, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4328, 0, 3, 453, 3608, 729, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4391, 0, 3, 468, 3653, 750, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4454, 0, 3, 483, 3698, 771, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4517, 0, 3, 498, 3743, 792, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4580, 0, 3, 513, 3788, 813, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4643, 0, 3, 528, 3833, 834, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4706, 0, 3, 543, 3878, 855, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4769, 0, 3, 588, 4013, 897, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4832, 0, 3, 603, 4058, 918, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4895, 0, 3, 618, 4103, 939, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4958, 0, 3, 633, 4148, 960, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5021, 0, 3, 648, 4193, 981, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5084, 0, 3, 663, 4238, 1002, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5147, 0, 3, 678, 4283, 1023, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5210, 0, 3, 708, 4328, 1044, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5294, 0, 3, 729, 4391, 1072, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5378, 0, 3, 750, 4454, 1100, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5462, 0, 3, 771, 4517, 1128, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5546, 0, 3, 792, 4580, 1156, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5630, 0, 3, 813, 4643, 1184, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5714, 0, 3, 834, 4706, 1212, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5798, 0, 3, 876, 4769, 1240, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5882, 0, 3, 897, 4832, 1268, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5966, 0, 3, 918, 4895, 1296, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6050, 0, 3, 939, 4958, 1324, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6134, 0, 3, 960, 5021, 1352, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6218, 0, 3, 981, 5084, 1380, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6302, 0, 3, 1002, 5147, 1408, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6386, 0, 3, 1072, 5378, 1472, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6494, 0, 3, 1100, 5462, 1508, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6602, 0, 3, 1128, 5546, 1544, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6710, 0, 3, 1156, 5630, 1580, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6818, 0, 3, 1184, 5714, 1616, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 6926, 0, 3, 1268, 5966, 1688, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7034, 0, 3, 1296, 6050, 1724, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7142, 0, 3, 1324, 6134, 1760, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7250, 0, 3, 1352, 6218, 1796, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7358, 0, 3, 1380, 6302, 1832, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7466, 0, 3, 1436, 6386, 1868, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7601, 0, 3, 1472, 6494, 1913, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7736, 0, 3, 1508, 6602, 1958, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 7871, 0, 3, 1544, 6710, 2003, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8006, 0, 3, 1580, 6818, 2048, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8141, 0, 3, 1652, 6926, 2093, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8276, 0, 3, 1688, 7034, 2138, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8411, 0, 3, 1724, 7142, 2183, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8546, 0, 3, 1760, 7250, 2228, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8681, 0, 3, 1796, 7358, 2273, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 8816, 3, 8, 9, 2321, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 8822, 3, 9, 10, 2324, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8828, 3, 10, 11, 2327, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8834, 3, 11, 12, 2330, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8840, 3, 12, 13, 2333, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8846, 3, 13, 14, 2336, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8852, 3, 14, 15, 2339, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8858, 3, 15, 16, 2342, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8864, 3, 16, 17, 2345, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8870, 3, 17, 18, 2348, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8876, 3, 22, 23, 2354, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8882, 3, 23, 24, 2357, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8888, 3, 24, 25, 2360, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8894, 3, 25, 26, 2363, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8900, 3, 26, 27, 2366, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8906, 3, 27, 28, 2369, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8912, 3, 28, 29, 2372, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8918, 3, 29, 30, 2375, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8924, 3, 30, 31, 2378, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 8930, 3, 31, 32, 2381, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 8936, 0, 3, 2318, 8816, 2393, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8954, 0, 3, 2321, 8822, 2402, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8972, 0, 3, 2324, 8828, 2411, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 8990, 0, 3, 2327, 8834, 2420, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9008, 0, 3, 2330, 8840, 2429, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9026, 0, 3, 2333, 8846, 2438, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9044, 0, 3, 2336, 8852, 2447, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9062, 0, 3, 2339, 8858, 2456, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9080, 0, 3, 2342, 8864, 2465, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9098, 0, 3, 2345, 8870, 2474, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9116, 0, 3, 2351, 8876, 2492, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9134, 0, 3, 2354, 8882, 2501, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9152, 0, 3, 2357, 8888, 2510, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9170, 0, 3, 2360, 8894, 2519, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9188, 0, 3, 2363, 8900, 2528, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9206, 0, 3, 2366, 8906, 2537, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9224, 0, 3, 2369, 8912, 2546, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9242, 0, 3, 2372, 8918, 2555, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9260, 0, 3, 2375, 8924, 2564, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9278, 0, 3, 2378, 8930, 2573, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 9296, 0, 3, 2393, 8954, 106, 112, 2618,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9332, 0, 3, 2402, 8972, 112, 118, 2636,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9368, 0, 3, 2411, 8990, 118, 124, 2654,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9404, 0, 3, 2420, 9008, 124, 130, 2672,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9440, 0, 3, 2429, 9026, 130, 136, 2690,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9476, 0, 3, 2438, 9044, 136, 142, 2708,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9512, 0, 3, 2447, 9062, 142, 148, 2726,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9548, 0, 3, 2456, 9080, 148, 154, 2744,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9584, 0, 3, 2465, 9098, 154, 160, 2762,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9620, 0, 3, 2492, 9134, 172, 178, 2816,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9656, 0, 3, 2501, 9152, 178, 184, 2834,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9692, 0, 3, 2510, 9170, 184, 190, 2852,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9728, 0, 3, 2519, 9188, 190, 196, 2870,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9764, 0, 3, 2528, 9206, 196, 202, 2888,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9800, 0, 3, 2537, 9224, 202, 208, 2906,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9836, 0, 3, 2546, 9242, 208, 214, 2924,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9872, 0, 3, 2555, 9260, 214, 220, 2942,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 9908, 0, 3, 2564, 9278, 220, 226, 2960,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 9944, 0, 3, 2618, 9332, 238, 248, 3008,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10004, 0, 3, 2636, 9368, 248, 258, 3038,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10064, 0, 3, 2654, 9404, 258, 268, 3068,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10124, 0, 3, 2672, 9440, 268, 278, 3098,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10184, 0, 3, 2690, 9476, 278, 288, 3128,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10244, 0, 3, 2708, 9512, 288, 298, 3158,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10304, 0, 3, 2726, 9548, 298, 308, 3188,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10364, 0, 3, 2744, 9584, 308, 318, 3218,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10424, 0, 3, 2816, 9656, 338, 348, 3278,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10484, 0, 3, 2834, 9692, 348, 358, 3308,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10544, 0, 3, 2852, 9728, 358, 368, 3338,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10604, 0, 3, 2870, 9764, 368, 378, 3368,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10664, 0, 3, 2888, 9800, 378, 388, 3398,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10724, 0, 3, 2906, 9836, 388, 398, 3428,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10784, 0, 3, 2924, 9872, 398, 408, 3458,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 10844, 0, 3, 2942, 9908, 408, 418, 3488,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10904, 0, 3, 9296, 9332, 3008, 10004,
                                                 438, 453, 3608, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 10994, 0, 3, 9332, 9368, 3038, 10064,
                                                 453, 468, 3653, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11084, 0, 3, 9368, 9404, 3068, 10124,
                                                 468, 483, 3698, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11174, 0, 3, 9404, 9440, 3098, 10184,
                                                 483, 498, 3743, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11264, 0, 3, 9440, 9476, 3128, 10244,
                                                 498, 513, 3788, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11354, 0, 3, 9476, 9512, 3158, 10304,
                                                 513, 528, 3833, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11444, 0, 3, 9512, 9548, 3188, 10364,
                                                 528, 543, 3878, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11534, 0, 3, 9620, 9656, 3278, 10484,
                                                 573, 588, 4013, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11624, 0, 3, 9656, 9692, 3308, 10544,
                                                 588, 603, 4058, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11714, 0, 3, 9692, 9728, 3338, 10604,
                                                 603, 618, 4103, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11804, 0, 3, 9728, 9764, 3368, 10664,
                                                 618, 633, 4148, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11894, 0, 3, 9764, 9800, 3398, 10724,
                                                 633, 648, 4193, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 11984, 0, 3, 9800, 9836, 3428, 10784,
                                                 648, 663, 4238, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12074, 0, 3, 9836, 9872, 3458, 10844,
                                                 663, 678, 4283, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12164, 0, 3, 9944, 10004, 3608, 10994,
                                                 708, 729, 4391, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12290, 0, 3, 10004, 10064, 3653, 11084,
                                                 729, 750, 4454, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12416, 0, 3, 10064, 10124, 3698, 11174,
                                                 750, 771, 4517, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12542, 0, 3, 10124, 10184, 3743, 11264,
                                                 771, 792, 4580, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12668, 0, 3, 10184, 10244, 3788, 11354,
                                                 792, 813, 4643, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12794, 0, 3, 10244, 10304, 3833, 11444,
                                                 813, 834, 4706, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 12920, 0, 3, 10424, 10484, 4013, 11624,
                                                 876, 897, 4832, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13046, 0, 3, 10484, 10544, 4058, 11714,
                                                 897, 918, 4895, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13172, 0, 3, 10544, 10604, 4103, 11804,
                                                 918, 939, 4958, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13298, 0, 3, 10604, 10664, 4148, 11894,
                                                 939, 960, 5021, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13424, 0, 3, 10664, 10724, 4193, 11984,
                                                 960, 981, 5084, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13550, 0, 3, 10724, 10784, 4238, 12074,
                                                 981, 1002, 5147, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13676, 0, 3, 10904, 10994, 4391, 12290,
                                                 1044, 1072, 5378, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13844, 0, 3, 10994, 11084, 4454, 12416,
                                                 1072, 1100, 5462, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14012, 0, 3, 11084, 11174, 4517, 12542,
                                                 1100, 1128, 5546, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14180, 0, 3, 11174, 11264, 4580, 12668,
                                                 1128, 1156, 5630, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14348, 0, 3, 11264, 11354, 4643, 12794,
                                                 1156, 1184, 5714, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14516, 0, 3, 11534, 11624, 4832, 13046,
                                                 1240, 1268, 5966, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14684, 0, 3, 11624, 11714, 4895, 13172,
                                                 1268, 1296, 6050, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 14852, 0, 3, 11714, 11804, 4958, 13298,
                                                 1296, 1324, 6134, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15020, 0, 3, 11804, 11894, 5021, 13424,
                                                 1324, 1352, 6218, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15188, 0, 3, 11894, 11984, 5084, 13550,
                                                 1352, 1380, 6302, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15356, 0, 3, 12164, 12290, 5378, 13844,
                                                 1436, 1472, 6494, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15572, 0, 3, 12290, 12416, 5462, 14012,
                                                 1472, 1508, 6602, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 15788, 0, 3, 12416, 12542, 5546, 14180,
                                                 1508, 1544, 6710, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16004, 0, 3, 12542, 12668, 5630, 14348,
                                                 1544, 1580, 6818, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16220, 0, 3, 12920, 13046, 5966, 14684,
                                                 1652, 1688, 7034, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16436, 0, 3, 13046, 13172, 6050, 14852,
                                                 1688, 1724, 7142, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16652, 0, 3, 13172, 13298, 6134, 15020,
                                                 1724, 1760, 7250, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 16868, 0, 3, 13298, 13424, 6218, 15188,
                                                 1760, 1796, 7358, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 17084, 0, 3, 13676, 13844, 6494, 15572,
                                                 1868, 1913, 7736, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 17354, 0, 3, 13844, 14012, 6602, 15788,
                                                 1913, 1958, 7871, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 17624, 0, 3, 14012, 14180, 6710, 16004,
                                                 1958, 2003, 8006, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 17894, 0, 3, 14516, 14684, 7034, 16436,
                                                 2093, 2138, 8411, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 18164, 0, 3, 14684, 14852, 7142, 16652,
                                                 2138, 2183, 8546, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 18434, 0, 3, 14852, 15020, 7250, 16868,
                                                 2183, 2228, 8681, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18704, 3, 2318, 2321, 8822, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18714, 3, 2321, 2324, 8828, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18724, 3, 2324, 2327, 8834, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18734, 3, 2327, 2330, 8840, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18744, 3, 2330, 2333, 8846, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18754, 3, 2333, 2336, 8852, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18764, 3, 2336, 2339, 8858, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18774, 3, 2339, 2342, 8864, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18784, 3, 2342, 2345, 8870, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18794, 3, 2351, 2354, 8882, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18804, 3, 2354, 2357, 8888, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18814, 3, 2357, 2360, 8894, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18824, 3, 2360, 2363, 8900, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18834, 3, 2363, 2366, 8906, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18844, 3, 2366, 2369, 8912, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18854, 3, 2369, 2372, 8918, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18864, 3, 2372, 2375, 8924, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 18874, 3, 2375, 2378, 8930, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 18884, 0, 3, 8816, 18704, 8954, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18914, 0, 3, 8822, 18714, 8972, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18944, 0, 3, 8828, 18724, 8990, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 18974, 0, 3, 8834, 18734, 9008, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19004, 0, 3, 8840, 18744, 9026, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19034, 0, 3, 8846, 18754, 9044, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19064, 0, 3, 8852, 18764, 9062, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19094, 0, 3, 8858, 18774, 9080, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19124, 0, 3, 8864, 18784, 9098, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19154, 0, 3, 8876, 18794, 9134, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19184, 0, 3, 8882, 18804, 9152, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19214, 0, 3, 8888, 18814, 9170, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19244, 0, 3, 8894, 18824, 9188, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19274, 0, 3, 8900, 18834, 9206, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19304, 0, 3, 8906, 18844, 9224, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19334, 0, 3, 8912, 18854, 9242, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19364, 0, 3, 8918, 18864, 9260, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 19394, 0, 3, 8924, 18874, 9278, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 19424, 0, 3, 8936, 18884, 2582, 2600,
                                                 9296, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19484, 0, 3, 8954, 18914, 2600, 2618,
                                                 9332, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19544, 0, 3, 8972, 18944, 2618, 2636,
                                                 9368, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19604, 0, 3, 8990, 18974, 2636, 2654,
                                                 9404, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19664, 0, 3, 9008, 19004, 2654, 2672,
                                                 9440, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19724, 0, 3, 9026, 19034, 2672, 2690,
                                                 9476, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19784, 0, 3, 9044, 19064, 2690, 2708,
                                                 9512, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19844, 0, 3, 9062, 19094, 2708, 2726,
                                                 9548, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19904, 0, 3, 9080, 19124, 2726, 2744,
                                                 9584, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 19964, 0, 3, 9116, 19154, 2780, 2798,
                                                 9620, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20024, 0, 3, 9134, 19184, 2798, 2816,
                                                 9656, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20084, 0, 3, 9152, 19214, 2816, 2834,
                                                 9692, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20144, 0, 3, 9170, 19244, 2834, 2852,
                                                 9728, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20204, 0, 3, 9188, 19274, 2852, 2870,
                                                 9764, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20264, 0, 3, 9206, 19304, 2870, 2888,
                                                 9800, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20324, 0, 3, 9224, 19334, 2888, 2906,
                                                 9836, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20384, 0, 3, 9242, 19364, 2906, 2924,
                                                 9872, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 20444, 0, 3, 9260, 19394, 2924, 2942,
                                                 9908, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20504, 0, 3, 9332, 19544, 2978, 3008,
                                                 10004, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20604, 0, 3, 9368, 19604, 3008, 3038,
                                                 10064, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20704, 0, 3, 9404, 19664, 3038, 3068,
                                                 10124, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20804, 0, 3, 9440, 19724, 3068, 3098,
                                                 10184, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 20904, 0, 3, 9476, 19784, 3098, 3128,
                                                 10244, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21004, 0, 3, 9512, 19844, 3128, 3158,
                                                 10304, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21104, 0, 3, 9548, 19904, 3158, 3188,
                                                 10364, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21204, 0, 3, 9656, 20084, 3248, 3278,
                                                 10484, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21304, 0, 3, 9692, 20144, 3278, 3308,
                                                 10544, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21404, 0, 3, 9728, 20204, 3308, 3338,
                                                 10604, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21504, 0, 3, 9764, 20264, 3338, 3368,
                                                 10664, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21604, 0, 3, 9800, 20324, 3368, 3398,
                                                 10724, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21704, 0, 3, 9836, 20384, 3398, 3428,
                                                 10784, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 21804, 0, 3, 9872, 20444, 3428, 3458,
                                                 10844, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 21904, 0, 3, 19424, 19484, 9944, 20504,
                                                 3518, 3563, 10904, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22054, 0, 3, 19484, 19544, 10004, 20604,
                                                 3563, 3608, 10994, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22204, 0, 3, 19544, 19604, 10064, 20704,
                                                 3608, 3653, 11084, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22354, 0, 3, 19604, 19664, 10124, 20804,
                                                 3653, 3698, 11174, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22504, 0, 3, 19664, 19724, 10184, 20904,
                                                 3698, 3743, 11264, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22654, 0, 3, 19724, 19784, 10244, 21004,
                                                 3743, 3788, 11354, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22804, 0, 3, 19784, 19844, 10304, 21104,
                                                 3788, 3833, 11444, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 22954, 0, 3, 19964, 20024, 10424, 21204,
                                                 3923, 3968, 11534, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23104, 0, 3, 20024, 20084, 10484, 21304,
                                                 3968, 4013, 11624, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23254, 0, 3, 20084, 20144, 10544, 21404,
                                                 4013, 4058, 11714, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23404, 0, 3, 20144, 20204, 10604, 21504,
                                                 4058, 4103, 11804, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23554, 0, 3, 20204, 20264, 10664, 21604,
                                                 4103, 4148, 11894, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23704, 0, 3, 20264, 20324, 10724, 21704,
                                                 4148, 4193, 11984, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 23854, 0, 3, 20324, 20384, 10784, 21804,
                                                 4193, 4238, 12074, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24004, 0, 3, 20504, 20604, 10994, 22204,
                                                 4328, 4391, 12290, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24214, 0, 3, 20604, 20704, 11084, 22354,
                                                 4391, 4454, 12416, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24424, 0, 3, 20704, 20804, 11174, 22504,
                                                 4454, 4517, 12542, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24634, 0, 3, 20804, 20904, 11264, 22654,
                                                 4517, 4580, 12668, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 24844, 0, 3, 20904, 21004, 11354, 22804,
                                                 4580, 4643, 12794, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25054, 0, 3, 21204, 21304, 11624, 23254,
                                                 4769, 4832, 13046, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25264, 0, 3, 21304, 21404, 11714, 23404,
                                                 4832, 4895, 13172, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25474, 0, 3, 21404, 21504, 11804, 23554,
                                                 4895, 4958, 13298, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25684, 0, 3, 21504, 21604, 11894, 23704,
                                                 4958, 5021, 13424, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 25894, 0, 3, 21604, 21704, 11984, 23854,
                                                 5021, 5084, 13550, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26104, 0, 3, 21904, 22054, 12164, 24004,
                                                 5210, 5294, 13676, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26384, 0, 3, 22054, 22204, 12290, 24214,
                                                 5294, 5378, 13844, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26664, 0, 3, 22204, 22354, 12416, 24424,
                                                 5378, 5462, 14012, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 26944, 0, 3, 22354, 22504, 12542, 24634,
                                                 5462, 5546, 14180, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27224, 0, 3, 22504, 22654, 12668, 24844,
                                                 5546, 5630, 14348, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27504, 0, 3, 22954, 23104, 12920, 25054,
                                                 5798, 5882, 14516, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 27784, 0, 3, 23104, 23254, 13046, 25264,
                                                 5882, 5966, 14684, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 28064, 0, 3, 23254, 23404, 13172, 25474,
                                                 5966, 6050, 14852, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 28344, 0, 3, 23404, 23554, 13298, 25684,
                                                 6050, 6134, 15020, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 28624, 0, 3, 23554, 23704, 13424, 25894,
                                                 6134, 6218, 15188, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 28904, 0, 3, 24004, 24214, 13844, 26664,
                                                 6386, 6494, 15572, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29264, 0, 3, 24214, 24424, 14012, 26944,
                                                 6494, 6602, 15788, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29624, 0, 3, 24424, 24634, 14180, 27224,
                                                 6602, 6710, 16004, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 29984, 0, 3, 25054, 25264, 14684, 28064,
                                                 6926, 7034, 16436, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 30344, 0, 3, 25264, 25474, 14852, 28344,
                                                 7034, 7142, 16652, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 30704, 0, 3, 25474, 25684, 15020, 28624,
                                                 7142, 7250, 16868, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 31064, 0, 3, 26104, 26384, 15356, 28904,
                                                 7466, 7601, 17084, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 31514, 0, 3, 26384, 26664, 15572, 29264,
                                                 7601, 7736, 17354, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 31964, 0, 3, 26664, 26944, 15788, 29624,
                                                 7736, 7871, 17624, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 32414, 0, 3, 27504, 27784, 16220, 29984,
                                                 8141, 8276, 17894, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 32864, 0, 3, 27784, 28064, 16436, 30344,
                                                 8276, 8411, 18164, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 33314, 0, 3, 28064, 28344, 16652, 30704,
                                                 8411, 8546, 18434, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33764, 3, 8816, 8822, 18714, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33779, 3, 8822, 8828, 18724, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33794, 3, 8828, 8834, 18734, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33809, 3, 8834, 8840, 18744, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33824, 3, 8840, 8846, 18754, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33839, 3, 8846, 8852, 18764, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33854, 3, 8852, 8858, 18774, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33869, 3, 8858, 8864, 18784, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33884, 3, 8876, 8882, 18804, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33899, 3, 8882, 8888, 18814, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33914, 3, 8888, 8894, 18824, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33929, 3, 8894, 8900, 18834, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33944, 3, 8900, 8906, 18844, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33959, 3, 8906, 8912, 18854, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33974, 3, 8912, 8918, 18864, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 33989, 3, 8918, 8924, 18874, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 34004, 0, 3, 18704, 33764, 18914, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34049, 0, 3, 18714, 33779, 18944, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34094, 0, 3, 18724, 33794, 18974, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34139, 0, 3, 18734, 33809, 19004, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34184, 0, 3, 18744, 33824, 19034, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34229, 0, 3, 18754, 33839, 19064, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34274, 0, 3, 18764, 33854, 19094, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34319, 0, 3, 18774, 33869, 19124, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34364, 0, 3, 18794, 33884, 19184, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34409, 0, 3, 18804, 33899, 19214, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34454, 0, 3, 18814, 33914, 19244, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34499, 0, 3, 18824, 33929, 19274, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34544, 0, 3, 18834, 33944, 19304, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34589, 0, 3, 18844, 33959, 19334, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34634, 0, 3, 18854, 33974, 19364, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 34679, 0, 3, 18864, 33989, 19394, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 34724, 0, 3, 18914, 34049, 9296, 9332,
                                                 19544, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34814, 0, 3, 18944, 34094, 9332, 9368,
                                                 19604, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34904, 0, 3, 18974, 34139, 9368, 9404,
                                                 19664, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 34994, 0, 3, 19004, 34184, 9404, 9440,
                                                 19724, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35084, 0, 3, 19034, 34229, 9440, 9476,
                                                 19784, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35174, 0, 3, 19064, 34274, 9476, 9512,
                                                 19844, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35264, 0, 3, 19094, 34319, 9512, 9548,
                                                 19904, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35354, 0, 3, 19184, 34409, 9620, 9656,
                                                 20084, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35444, 0, 3, 19214, 34454, 9656, 9692,
                                                 20144, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35534, 0, 3, 19244, 34499, 9692, 9728,
                                                 20204, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35624, 0, 3, 19274, 34544, 9728, 9764,
                                                 20264, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35714, 0, 3, 19304, 34589, 9764, 9800,
                                                 20324, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35804, 0, 3, 19334, 34634, 9800, 9836,
                                                 20384, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 35894, 0, 3, 19364, 34679, 9836, 9872,
                                                 20444, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 35984, 0, 3, 19544, 34814, 9944, 10004,
                                                 20604, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36134, 0, 3, 19604, 34904, 10004, 10064,
                                                 20704, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36284, 0, 3, 19664, 34994, 10064, 10124,
                                                 20804, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36434, 0, 3, 19724, 35084, 10124, 10184,
                                                 20904, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36584, 0, 3, 19784, 35174, 10184, 10244,
                                                 21004, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36734, 0, 3, 19844, 35264, 10244, 10304,
                                                 21104, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 36884, 0, 3, 20084, 35444, 10424, 10484,
                                                 21304, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37034, 0, 3, 20144, 35534, 10484, 10544,
                                                 21404, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37184, 0, 3, 20204, 35624, 10544, 10604,
                                                 21504, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37334, 0, 3, 20264, 35714, 10604, 10664,
                                                 21604, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37484, 0, 3, 20324, 35804, 10664, 10724,
                                                 21704, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 37634, 0, 3, 20384, 35894, 10724, 10784,
                                                 21804, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 37784, 0, 3, 34724, 34814, 20604, 36134,
                                                 10904, 10994, 22204, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38009, 0, 3, 34814, 34904, 20704, 36284,
                                                 10994, 11084, 22354, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38234, 0, 3, 34904, 34994, 20804, 36434,
                                                 11084, 11174, 22504, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38459, 0, 3, 34994, 35084, 20904, 36584,
                                                 11174, 11264, 22654, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38684, 0, 3, 35084, 35174, 21004, 36734,
                                                 11264, 11354, 22804, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 38909, 0, 3, 35354, 35444, 21304, 37034,
                                                 11534, 11624, 23254, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39134, 0, 3, 35444, 35534, 21404, 37184,
                                                 11624, 11714, 23404, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39359, 0, 3, 35534, 35624, 21504, 37334,
                                                 11714, 11804, 23554, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39584, 0, 3, 35624, 35714, 21604, 37484,
                                                 11804, 11894, 23704, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 39809, 0, 3, 35714, 35804, 21704, 37634,
                                                 11894, 11984, 23854, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 40034, 0, 3, 35984, 36134, 22204, 38009,
                                                 12164, 12290, 24214, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 40349, 0, 3, 36134, 36284, 22354, 38234,
                                                 12290, 12416, 24424, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 40664, 0, 3, 36284, 36434, 22504, 38459,
                                                 12416, 12542, 24634, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 40979, 0, 3, 36434, 36584, 22654, 38684,
                                                 12542, 12668, 24844, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41294, 0, 3, 36884, 37034, 23254, 39134,
                                                 12920, 13046, 25264, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41609, 0, 3, 37034, 37184, 23404, 39359,
                                                 13046, 13172, 25474, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 41924, 0, 3, 37184, 37334, 23554, 39584,
                                                 13172, 13298, 25684, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 42239, 0, 3, 37334, 37484, 23704, 39809,
                                                 13298, 13424, 25894, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 42554, 0, 3, 37784, 38009, 24214, 40349,
                                                 13676, 13844, 26664, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 42974, 0, 3, 38009, 38234, 24424, 40664,
                                                 13844, 14012, 26944, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 43394, 0, 3, 38234, 38459, 24634, 40979,
                                                 14012, 14180, 27224, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 43814, 0, 3, 38909, 39134, 25264, 41609,
                                                 14516, 14684, 28064, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 44234, 0, 3, 39134, 39359, 25474, 41924,
                                                 14684, 14852, 28344, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 44654, 0, 3, 39359, 39584, 25684, 42239,
                                                 14852, 15020, 28624, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 45074, 0, 3, 40034, 40349, 26664, 42974,
                                                 15356, 15572, 29264, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 45614, 0, 3, 40349, 40664, 26944, 43394,
                                                 15572, 15788, 29624, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 46154, 0, 3, 41294, 41609, 28064, 44234,
                                                 16220, 16436, 30344, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 46694, 0, 3, 41609, 41924, 28344, 44654,
                                                 16436, 16652, 30704, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 47234, 0, 3, 42554, 42974, 29264, 45614,
                                                 17084, 17354, 31964, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 47909, 0, 3, 43814, 44234, 30344, 46694,
                                                 17894, 18164, 33314, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48584, 3, 18704, 18714, 33779, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48605, 3, 18714, 18724, 33794, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48626, 3, 18724, 18734, 33809, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48647, 3, 18734, 18744, 33824, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48668, 3, 18744, 18754, 33839, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48689, 3, 18754, 18764, 33854, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48710, 3, 18764, 18774, 33869, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48731, 3, 18794, 18804, 33899, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48752, 3, 18804, 18814, 33914, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48773, 3, 18814, 18824, 33929, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48794, 3, 18824, 18834, 33944, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48815, 3, 18834, 18844, 33959, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48836, 3, 18844, 18854, 33974, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 48857, 3, 18854, 18864, 33989, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 48878, 0, 3, 33764, 48584, 34049, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 48941, 0, 3, 33779, 48605, 34094, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49004, 0, 3, 33794, 48626, 34139, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49067, 0, 3, 33809, 48647, 34184, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49130, 0, 3, 33824, 48668, 34229, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49193, 0, 3, 33839, 48689, 34274, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49256, 0, 3, 33854, 48710, 34319, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49319, 0, 3, 33884, 48731, 34409, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49382, 0, 3, 33899, 48752, 34454, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49445, 0, 3, 33914, 48773, 34499, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49508, 0, 3, 33929, 48794, 34544, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49571, 0, 3, 33944, 48815, 34589, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49634, 0, 3, 33959, 48836, 34634, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 49697, 0, 3, 33974, 48857, 34679, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 49760, 0, 3, 34004, 48878, 19424, 19484,
                                                 34724, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 49886, 0, 3, 34049, 48941, 19484, 19544,
                                                 34814, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50012, 0, 3, 34094, 49004, 19544, 19604,
                                                 34904, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50138, 0, 3, 34139, 49067, 19604, 19664,
                                                 34994, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50264, 0, 3, 34184, 49130, 19664, 19724,
                                                 35084, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50390, 0, 3, 34229, 49193, 19724, 19784,
                                                 35174, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50516, 0, 3, 34274, 49256, 19784, 19844,
                                                 35264, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50642, 0, 3, 34364, 49319, 19964, 20024,
                                                 35354, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50768, 0, 3, 34409, 49382, 20024, 20084,
                                                 35444, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 50894, 0, 3, 34454, 49445, 20084, 20144,
                                                 35534, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51020, 0, 3, 34499, 49508, 20144, 20204,
                                                 35624, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51146, 0, 3, 34544, 49571, 20204, 20264,
                                                 35714, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51272, 0, 3, 34589, 49634, 20264, 20324,
                                                 35804, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 51398, 0, 3, 34634, 49697, 20324, 20384,
                                                 35894, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 51524, 0, 3, 34814, 50012, 20504, 20604,
                                                 36134, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 51734, 0, 3, 34904, 50138, 20604, 20704,
                                                 36284, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 51944, 0, 3, 34994, 50264, 20704, 20804,
                                                 36434, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52154, 0, 3, 35084, 50390, 20804, 20904,
                                                 36584, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52364, 0, 3, 35174, 50516, 20904, 21004,
                                                 36734, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52574, 0, 3, 35444, 50894, 21204, 21304,
                                                 37034, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52784, 0, 3, 35534, 51020, 21304, 21404,
                                                 37184, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 52994, 0, 3, 35624, 51146, 21404, 21504,
                                                 37334, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 53204, 0, 3, 35714, 51272, 21504, 21604,
                                                 37484, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 53414, 0, 3, 35804, 51398, 21604, 21704,
                                                 37634, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 53624, 0, 3, 49760, 49886, 35984, 51524,
                                                 21904, 22054, 37784, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 53939, 0, 3, 49886, 50012, 36134, 51734,
                                                 22054, 22204, 38009, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 54254, 0, 3, 50012, 50138, 36284, 51944,
                                                 22204, 22354, 38234, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 54569, 0, 3, 50138, 50264, 36434, 52154,
                                                 22354, 22504, 38459, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 54884, 0, 3, 50264, 50390, 36584, 52364,
                                                 22504, 22654, 38684, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 55199, 0, 3, 50642, 50768, 36884, 52574,
                                                 22954, 23104, 38909, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 55514, 0, 3, 50768, 50894, 37034, 52784,
                                                 23104, 23254, 39134, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 55829, 0, 3, 50894, 51020, 37184, 52994,
                                                 23254, 23404, 39359, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 56144, 0, 3, 51020, 51146, 37334, 53204,
                                                 23404, 23554, 39584, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 56459, 0, 3, 51146, 51272, 37484, 53414,
                                                 23554, 23704, 39809, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 56774, 0, 3, 51524, 51734, 38009, 54254,
                                                 24004, 24214, 40349, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 57215, 0, 3, 51734, 51944, 38234, 54569,
                                                 24214, 24424, 40664, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 57656, 0, 3, 51944, 52154, 38459, 54884,
                                                 24424, 24634, 40979, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 58097, 0, 3, 52574, 52784, 39134, 55829,
                                                 25054, 25264, 41609, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 58538, 0, 3, 52784, 52994, 39359, 56144,
                                                 25264, 25474, 41924, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 58979, 0, 3, 52994, 53204, 39584, 56459,
                                                 25474, 25684, 42239, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 59420, 0, 3, 53624, 53939, 40034, 56774,
                                                 26104, 26384, 42554, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 60008, 0, 3, 53939, 54254, 40349, 57215,
                                                 26384, 26664, 42974, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 60596, 0, 3, 54254, 54569, 40664, 57656,
                                                 26664, 26944, 43394, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 61184, 0, 3, 55199, 55514, 41294, 58097,
                                                 27504, 27784, 43814, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 61772, 0, 3, 55514, 55829, 41609, 58538,
                                                 27784, 28064, 44234, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 62360, 0, 3, 55829, 56144, 41924, 58979,
                                                 28064, 28344, 44654, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 62948, 0, 3, 56774, 57215, 42974, 60596,
                                                 28904, 29264, 45614, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 63704, 0, 3, 58097, 58538, 44234, 62360,
                                                 29984, 30344, 46694, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 64460, 0, 3, 59420, 60008, 45074, 62948,
                                                 31064, 31514, 47234, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 65405, 0, 3, 61184, 61772, 46154, 63704,
                                                 32414, 32864, 47909, ncols, alpha, beta, p);

            simdgeo::geom_k_x(buffer, 66350, 61184, 65405, 1, 21, ncols, alpha);

            simdgeo::geom_k_y(buffer, 67106, 61184, 65405, 1, 21, ncols, alpha);

            simdgeo::geom_k_z(buffer, 67862, 61184, 65405, 1, 21, ncols, alpha);

            simdgeo::geom_k_x(buffer, 68618, 59420, 64460, 1, 21, ncols, alpha);

            simdgeo::geom_k_y(buffer, 69374, 59420, 64460, 1, 21, ncols, alpha);

            simdgeo::geom_k_z(buffer, 70130, 59420, 64460, 1, 21, ncols, alpha);

            simdfunc::contract_primitives(buffer, 70886, 68618, 2268, ncols);

            simdfunc::contract_primitives(buffer, 73154, 66350, 2268, ncols);
        }
    }

    simdtrf::transform_h_inner(buffer, 75422, 73154, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 75422, 11, nmax);

    simdtrf::transform_h_inner(buffer, 75422, 73910, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 165 * nvalues, nvalues, buffer, 75422, 11, nmax);

    simdtrf::transform_h_inner(buffer, 75422, 74666, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 330 * nvalues, nvalues, buffer, 75422, 11, nmax);

    simdtrf::transform_h_inner(buffer, 75422, 70886, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 495 * nvalues, nvalues, buffer, 75422, 11, nmax);

    simdtrf::transform_h_inner(buffer, 75422, 71642, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 660 * nvalues, nvalues, buffer, 75422, 11, nmax);

    simdtrf::transform_h_inner(buffer, 75422, 72398, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 825 * nvalues, nvalues, buffer, 75422, 11, nmax);
}

}  // namespace simdt2ceri
