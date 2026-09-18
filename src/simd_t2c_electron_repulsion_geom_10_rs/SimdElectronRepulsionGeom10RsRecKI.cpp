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


#include "SimdElectronRepulsionGeom10RsRecKI.hpp"

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
#include "SimdGeometryK1.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ki_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_ki_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 118828, 112312, 6048, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 14, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 22, 14, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 38, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 41, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 44, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 47, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 50, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 53, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 56, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 59, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 62, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 65, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 68, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 71, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 74, 0, 21, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 77, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 80, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 83, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 86, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 89, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 92, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 95, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 98, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 101, 0, 33, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 104, 0, 34, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 107, 0, 35, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 110, 0, 36, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 113, 0, 37, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 116, 0, 7, 8, 38, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 122, 0, 8, 9, 41, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 128, 0, 9, 10, 44, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 134, 0, 10, 11, 47, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 140, 0, 11, 12, 50, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 146, 0, 12, 13, 53, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 152, 0, 13, 14, 56, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 158, 0, 14, 15, 59, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 164, 0, 15, 16, 62, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 170, 0, 16, 17, 65, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 176, 0, 17, 18, 68, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 182, 0, 18, 19, 71, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 188, 0, 19, 20, 74, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 194, 0, 23, 24, 77, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 200, 0, 24, 25, 80, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 206, 0, 25, 26, 83, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 212, 0, 26, 27, 86, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 218, 0, 27, 28, 89, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 224, 0, 28, 29, 92, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 230, 0, 29, 30, 95, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 236, 0, 30, 31, 98, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 242, 0, 31, 32, 101, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 248, 0, 32, 33, 104, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 254, 0, 33, 34, 107, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 260, 0, 34, 35, 110, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 266, 0, 35, 36, 113, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 272, 0, 38, 41, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 282, 0, 41, 44, 134, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 292, 0, 44, 47, 140, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 302, 0, 47, 50, 146, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 312, 0, 50, 53, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 322, 0, 53, 56, 158, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 332, 0, 56, 59, 164, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 342, 0, 59, 62, 170, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 352, 0, 62, 65, 176, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 362, 0, 65, 68, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 372, 0, 68, 71, 188, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 382, 0, 77, 80, 206, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 392, 0, 80, 83, 212, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 402, 0, 83, 86, 218, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 412, 0, 86, 89, 224, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 422, 0, 89, 92, 230, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 432, 0, 92, 95, 236, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 442, 0, 95, 98, 242, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 452, 0, 98, 101, 248, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 462, 0, 101, 104, 254, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 472, 0, 104, 107, 260, ncols, alpha,
                                                 beta, p);

            compute_prim_fs_electron_repulsion_0(buffer, 482, 0, 107, 110, 266, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 492, 0, 116, 122, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 507, 0, 122, 128, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 522, 0, 128, 134, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 537, 0, 134, 140, 302, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 552, 0, 140, 146, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 567, 0, 146, 152, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 582, 0, 152, 158, 332, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 597, 0, 158, 164, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 612, 0, 164, 170, 352, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 627, 0, 170, 176, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 642, 0, 176, 182, 372, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 657, 0, 194, 200, 382, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 672, 0, 200, 206, 392, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 687, 0, 206, 212, 402, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 702, 0, 212, 218, 412, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 717, 0, 218, 224, 422, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 732, 0, 224, 230, 432, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 747, 0, 230, 236, 442, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 762, 0, 236, 242, 452, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 777, 0, 242, 248, 462, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 792, 0, 248, 254, 472, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 807, 0, 254, 260, 482, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 822, 0, 272, 282, 522, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 843, 0, 282, 292, 537, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 864, 0, 292, 302, 552, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 885, 0, 302, 312, 567, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 906, 0, 312, 322, 582, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 927, 0, 322, 332, 597, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 948, 0, 332, 342, 612, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 969, 0, 342, 352, 627, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 990, 0, 352, 362, 642, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1011, 0, 382, 392, 687, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1032, 0, 392, 402, 702, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1053, 0, 402, 412, 717, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1074, 0, 412, 422, 732, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1095, 0, 422, 432, 747, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1116, 0, 432, 442, 762, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1137, 0, 442, 452, 777, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1158, 0, 452, 462, 792, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 1179, 0, 462, 472, 807, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1200, 0, 492, 507, 822, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1228, 0, 507, 522, 843, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1256, 0, 522, 537, 864, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1284, 0, 537, 552, 885, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1312, 0, 552, 567, 906, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1340, 0, 567, 582, 927, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1368, 0, 582, 597, 948, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1396, 0, 597, 612, 969, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1424, 0, 612, 627, 990, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1452, 0, 657, 672, 1011, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1480, 0, 672, 687, 1032, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1508, 0, 687, 702, 1053, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1536, 0, 702, 717, 1074, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1564, 0, 717, 732, 1095, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1592, 0, 732, 747, 1116, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1620, 0, 747, 762, 1137, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1648, 0, 762, 777, 1158, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 1676, 0, 777, 792, 1179, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1704, 0, 822, 843, 1256, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1740, 0, 843, 864, 1284, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1776, 0, 864, 885, 1312, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1812, 0, 885, 906, 1340, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1848, 0, 906, 927, 1368, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1884, 0, 927, 948, 1396, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1920, 0, 948, 969, 1424, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1956, 0, 1011, 1032, 1508, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1992, 0, 1032, 1053, 1536, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2028, 0, 1053, 1074, 1564, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2064, 0, 1074, 1095, 1592, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2100, 0, 1095, 1116, 1620, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2136, 0, 1116, 1137, 1648, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 2172, 0, 1137, 1158, 1676, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2208, 0, 1200, 1228, 1704, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2253, 0, 1228, 1256, 1740, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2298, 0, 1256, 1284, 1776, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2343, 0, 1284, 1312, 1812, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2388, 0, 1312, 1340, 1848, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2433, 0, 1340, 1368, 1884, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2478, 0, 1368, 1396, 1920, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2523, 0, 1452, 1480, 1956, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2568, 0, 1480, 1508, 1992, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2613, 0, 1508, 1536, 2028, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2658, 0, 1536, 1564, 2064, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2703, 0, 1564, 1592, 2100, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2748, 0, 1592, 1620, 2136, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 2793, 0, 1620, 1648, 2172, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 2838, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2841, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2844, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2847, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2850, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2853, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2856, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2859, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2862, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2865, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2868, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2871, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2874, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2877, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2880, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2883, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2886, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2889, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2892, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2895, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2898, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2901, 3, 35, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2904, 3, 36, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 2907, 3, 37, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 2910, 3, 9, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2919, 3, 10, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2928, 3, 11, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2937, 3, 12, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2946, 3, 13, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2955, 3, 14, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2964, 3, 15, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2973, 3, 16, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2982, 3, 17, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 2991, 3, 18, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3000, 3, 19, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3009, 3, 20, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3018, 3, 25, 80, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3027, 3, 26, 83, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3036, 3, 27, 86, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3045, 3, 28, 89, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3054, 3, 29, 92, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3063, 3, 30, 95, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3072, 3, 31, 98, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3081, 3, 32, 101, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3090, 3, 33, 104, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3099, 3, 34, 107, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3108, 3, 35, 110, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 3117, 3, 36, 113, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3126, 0, 3, 41, 2919, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3144, 0, 3, 44, 2928, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3162, 0, 3, 47, 2937, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3180, 0, 3, 50, 2946, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3198, 0, 3, 53, 2955, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3216, 0, 3, 56, 2964, 158, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3234, 0, 3, 59, 2973, 164, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3252, 0, 3, 62, 2982, 170, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3270, 0, 3, 65, 2991, 176, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3288, 0, 3, 68, 3000, 182, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3306, 0, 3, 71, 3009, 188, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3324, 0, 3, 80, 3027, 206, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3342, 0, 3, 83, 3036, 212, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3360, 0, 3, 86, 3045, 218, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3378, 0, 3, 89, 3054, 224, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3396, 0, 3, 92, 3063, 230, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3414, 0, 3, 95, 3072, 236, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3432, 0, 3, 98, 3081, 242, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3450, 0, 3, 101, 3090, 248, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3468, 0, 3, 104, 3099, 254, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3486, 0, 3, 107, 3108, 260, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 3504, 0, 3, 110, 3117, 266, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3522, 0, 3, 128, 3144, 282, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3552, 0, 3, 134, 3162, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3582, 0, 3, 140, 3180, 302, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3612, 0, 3, 146, 3198, 312, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3642, 0, 3, 152, 3216, 322, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3672, 0, 3, 158, 3234, 332, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3702, 0, 3, 164, 3252, 342, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3732, 0, 3, 170, 3270, 352, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3762, 0, 3, 176, 3288, 362, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3792, 0, 3, 182, 3306, 372, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3822, 0, 3, 206, 3342, 392, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3852, 0, 3, 212, 3360, 402, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3882, 0, 3, 218, 3378, 412, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3912, 0, 3, 224, 3396, 422, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3942, 0, 3, 230, 3414, 432, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 3972, 0, 3, 236, 3432, 442, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4002, 0, 3, 242, 3450, 452, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4032, 0, 3, 248, 3468, 462, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4062, 0, 3, 254, 3486, 472, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 4092, 0, 3, 260, 3504, 482, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4122, 0, 3, 282, 3552, 522, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4167, 0, 3, 292, 3582, 537, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4212, 0, 3, 302, 3612, 552, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4257, 0, 3, 312, 3642, 567, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4302, 0, 3, 322, 3672, 582, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4347, 0, 3, 332, 3702, 597, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4392, 0, 3, 342, 3732, 612, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4437, 0, 3, 352, 3762, 627, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4482, 0, 3, 362, 3792, 642, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4527, 0, 3, 392, 3852, 687, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4572, 0, 3, 402, 3882, 702, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4617, 0, 3, 412, 3912, 717, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4662, 0, 3, 422, 3942, 732, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4707, 0, 3, 432, 3972, 747, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4752, 0, 3, 442, 4002, 762, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4797, 0, 3, 452, 4032, 777, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4842, 0, 3, 462, 4062, 792, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 4887, 0, 3, 472, 4092, 807, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4932, 0, 3, 522, 4167, 843, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4995, 0, 3, 537, 4212, 864, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5058, 0, 3, 552, 4257, 885, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5121, 0, 3, 567, 4302, 906, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5184, 0, 3, 582, 4347, 927, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5247, 0, 3, 597, 4392, 948, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5310, 0, 3, 612, 4437, 969, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5373, 0, 3, 627, 4482, 990, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5436, 0, 3, 687, 4572, 1032, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5499, 0, 3, 702, 4617, 1053, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5562, 0, 3, 717, 4662, 1074, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5625, 0, 3, 732, 4707, 1095, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5688, 0, 3, 747, 4752, 1116, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5751, 0, 3, 762, 4797, 1137, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5814, 0, 3, 777, 4842, 1158, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 5877, 0, 3, 792, 4887, 1179, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5940, 0, 3, 843, 4995, 1256, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6024, 0, 3, 864, 5058, 1284, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6108, 0, 3, 885, 5121, 1312, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6192, 0, 3, 906, 5184, 1340, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6276, 0, 3, 927, 5247, 1368, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6360, 0, 3, 948, 5310, 1396, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6444, 0, 3, 969, 5373, 1424, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 6528, 0, 3, 1032, 5499, 1508, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6612, 0, 3, 1053, 5562, 1536, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6696, 0, 3, 1074, 5625, 1564, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6780, 0, 3, 1095, 5688, 1592, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6864, 0, 3, 1116, 5751, 1620, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 6948, 0, 3, 1137, 5814, 1648, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 7032, 0, 3, 1158, 5877, 1676, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7116, 0, 3, 1256, 6024, 1740, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7224, 0, 3, 1284, 6108, 1776, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7332, 0, 3, 1312, 6192, 1812, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7440, 0, 3, 1340, 6276, 1848, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7548, 0, 3, 1368, 6360, 1884, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7656, 0, 3, 1396, 6444, 1920, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7764, 0, 3, 1508, 6612, 1992, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7872, 0, 3, 1536, 6696, 2028, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 7980, 0, 3, 1564, 6780, 2064, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8088, 0, 3, 1592, 6864, 2100, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8196, 0, 3, 1620, 6948, 2136, ncols,
                                                 p);

            compute_prim_kp_electron_repulsion_0(buffer, 8304, 0, 3, 1648, 7032, 2172, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8412, 0, 3, 1740, 7224, 2298, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8547, 0, 3, 1776, 7332, 2343, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8682, 0, 3, 1812, 7440, 2388, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8817, 0, 3, 1848, 7548, 2433, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 8952, 0, 3, 1884, 7656, 2478, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9087, 0, 3, 1992, 7872, 2613, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9222, 0, 3, 2028, 7980, 2658, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9357, 0, 3, 2064, 8088, 2703, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9492, 0, 3, 2100, 8196, 2748, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 9627, 0, 3, 2136, 8304, 2793, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 9762, 3, 9, 10, 2841, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9768, 3, 10, 11, 2844, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9774, 3, 11, 12, 2847, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9780, 3, 12, 13, 2850, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9786, 3, 13, 14, 2853, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9792, 3, 14, 15, 2856, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9798, 3, 15, 16, 2859, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9804, 3, 16, 17, 2862, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9810, 3, 17, 18, 2865, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9816, 3, 18, 19, 2868, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9822, 3, 19, 20, 2871, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9828, 3, 25, 26, 2877, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9834, 3, 26, 27, 2880, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9840, 3, 27, 28, 2883, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9846, 3, 28, 29, 2886, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9852, 3, 29, 30, 2889, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9858, 3, 30, 31, 2892, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9864, 3, 31, 32, 2895, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9870, 3, 32, 33, 2898, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9876, 3, 33, 34, 2901, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9882, 3, 34, 35, 2904, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 9888, 3, 35, 36, 2907, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 9894, 0, 3, 2838, 9762, 2919, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9912, 0, 3, 2841, 9768, 2928, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9930, 0, 3, 2844, 9774, 2937, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9948, 0, 3, 2847, 9780, 2946, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9966, 0, 3, 2850, 9786, 2955, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 9984, 0, 3, 2853, 9792, 2964, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10002, 0, 3, 2856, 9798, 2973, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10020, 0, 3, 2859, 9804, 2982, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10038, 0, 3, 2862, 9810, 2991, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10056, 0, 3, 2865, 9816, 3000, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10074, 0, 3, 2868, 9822, 3009, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10092, 0, 3, 2874, 9828, 3027, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10110, 0, 3, 2877, 9834, 3036, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10128, 0, 3, 2880, 9840, 3045, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10146, 0, 3, 2883, 9846, 3054, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10164, 0, 3, 2886, 9852, 3063, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10182, 0, 3, 2889, 9858, 3072, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10200, 0, 3, 2892, 9864, 3081, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10218, 0, 3, 2895, 9870, 3090, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10236, 0, 3, 2898, 9876, 3099, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10254, 0, 3, 2901, 9882, 3108, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 10272, 0, 3, 2904, 9888, 3117, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 10290, 0, 3, 2910, 9894, 116, 122, 3126,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10326, 0, 3, 2919, 9912, 122, 128, 3144,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10362, 0, 3, 2928, 9930, 128, 134, 3162,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10398, 0, 3, 2937, 9948, 134, 140, 3180,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10434, 0, 3, 2946, 9966, 140, 146, 3198,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10470, 0, 3, 2955, 9984, 146, 152, 3216,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10506, 0, 3, 2964, 10002, 152, 158,
                                                 3234, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10542, 0, 3, 2973, 10020, 158, 164,
                                                 3252, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10578, 0, 3, 2982, 10038, 164, 170,
                                                 3270, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10614, 0, 3, 2991, 10056, 170, 176,
                                                 3288, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10650, 0, 3, 3000, 10074, 176, 182,
                                                 3306, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10686, 0, 3, 3018, 10092, 194, 200,
                                                 3324, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10722, 0, 3, 3027, 10110, 200, 206,
                                                 3342, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10758, 0, 3, 3036, 10128, 206, 212,
                                                 3360, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10794, 0, 3, 3045, 10146, 212, 218,
                                                 3378, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10830, 0, 3, 3054, 10164, 218, 224,
                                                 3396, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10866, 0, 3, 3063, 10182, 224, 230,
                                                 3414, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10902, 0, 3, 3072, 10200, 230, 236,
                                                 3432, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10938, 0, 3, 3081, 10218, 236, 242,
                                                 3450, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 10974, 0, 3, 3090, 10236, 242, 248,
                                                 3468, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 11010, 0, 3, 3099, 10254, 248, 254,
                                                 3486, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 11046, 0, 3, 3108, 10272, 254, 260,
                                                 3504, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11082, 0, 3, 3144, 10362, 272, 282,
                                                 3552, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11142, 0, 3, 3162, 10398, 282, 292,
                                                 3582, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11202, 0, 3, 3180, 10434, 292, 302,
                                                 3612, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11262, 0, 3, 3198, 10470, 302, 312,
                                                 3642, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11322, 0, 3, 3216, 10506, 312, 322,
                                                 3672, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11382, 0, 3, 3234, 10542, 322, 332,
                                                 3702, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11442, 0, 3, 3252, 10578, 332, 342,
                                                 3732, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11502, 0, 3, 3270, 10614, 342, 352,
                                                 3762, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11562, 0, 3, 3288, 10650, 352, 362,
                                                 3792, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11622, 0, 3, 3342, 10758, 382, 392,
                                                 3852, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11682, 0, 3, 3360, 10794, 392, 402,
                                                 3882, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11742, 0, 3, 3378, 10830, 402, 412,
                                                 3912, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11802, 0, 3, 3396, 10866, 412, 422,
                                                 3942, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11862, 0, 3, 3414, 10902, 422, 432,
                                                 3972, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11922, 0, 3, 3432, 10938, 432, 442,
                                                 4002, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 11982, 0, 3, 3450, 10974, 442, 452,
                                                 4032, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 12042, 0, 3, 3468, 11010, 452, 462,
                                                 4062, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 12102, 0, 3, 3486, 11046, 462, 472,
                                                 4092, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12162, 0, 3, 10290, 10326, 3522, 11082,
                                                 492, 507, 4122, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12252, 0, 3, 10326, 10362, 3552, 11142,
                                                 507, 522, 4167, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12342, 0, 3, 10362, 10398, 3582, 11202,
                                                 522, 537, 4212, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12432, 0, 3, 10398, 10434, 3612, 11262,
                                                 537, 552, 4257, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12522, 0, 3, 10434, 10470, 3642, 11322,
                                                 552, 567, 4302, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12612, 0, 3, 10470, 10506, 3672, 11382,
                                                 567, 582, 4347, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12702, 0, 3, 10506, 10542, 3702, 11442,
                                                 582, 597, 4392, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12792, 0, 3, 10542, 10578, 3732, 11502,
                                                 597, 612, 4437, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12882, 0, 3, 10578, 10614, 3762, 11562,
                                                 612, 627, 4482, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 12972, 0, 3, 10686, 10722, 3822, 11622,
                                                 657, 672, 4527, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13062, 0, 3, 10722, 10758, 3852, 11682,
                                                 672, 687, 4572, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13152, 0, 3, 10758, 10794, 3882, 11742,
                                                 687, 702, 4617, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13242, 0, 3, 10794, 10830, 3912, 11802,
                                                 702, 717, 4662, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13332, 0, 3, 10830, 10866, 3942, 11862,
                                                 717, 732, 4707, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13422, 0, 3, 10866, 10902, 3972, 11922,
                                                 732, 747, 4752, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13512, 0, 3, 10902, 10938, 4002, 11982,
                                                 747, 762, 4797, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13602, 0, 3, 10938, 10974, 4032, 12042,
                                                 762, 777, 4842, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 13692, 0, 3, 10974, 11010, 4062, 12102,
                                                 777, 792, 4887, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13782, 0, 3, 11082, 11142, 4167, 12342,
                                                 822, 843, 4995, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 13908, 0, 3, 11142, 11202, 4212, 12432,
                                                 843, 864, 5058, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14034, 0, 3, 11202, 11262, 4257, 12522,
                                                 864, 885, 5121, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14160, 0, 3, 11262, 11322, 4302, 12612,
                                                 885, 906, 5184, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14286, 0, 3, 11322, 11382, 4347, 12702,
                                                 906, 927, 5247, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14412, 0, 3, 11382, 11442, 4392, 12792,
                                                 927, 948, 5310, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14538, 0, 3, 11442, 11502, 4437, 12882,
                                                 948, 969, 5373, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14664, 0, 3, 11622, 11682, 4572, 13152,
                                                 1011, 1032, 5499, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14790, 0, 3, 11682, 11742, 4617, 13242,
                                                 1032, 1053, 5562, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 14916, 0, 3, 11742, 11802, 4662, 13332,
                                                 1053, 1074, 5625, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15042, 0, 3, 11802, 11862, 4707, 13422,
                                                 1074, 1095, 5688, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15168, 0, 3, 11862, 11922, 4752, 13512,
                                                 1095, 1116, 5751, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15294, 0, 3, 11922, 11982, 4797, 13602,
                                                 1116, 1137, 5814, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 15420, 0, 3, 11982, 12042, 4842, 13692,
                                                 1137, 1158, 5877, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15546, 0, 3, 12162, 12252, 4932, 13782,
                                                 1200, 1228, 5940, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15714, 0, 3, 12252, 12342, 4995, 13908,
                                                 1228, 1256, 6024, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 15882, 0, 3, 12342, 12432, 5058, 14034,
                                                 1256, 1284, 6108, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16050, 0, 3, 12432, 12522, 5121, 14160,
                                                 1284, 1312, 6192, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16218, 0, 3, 12522, 12612, 5184, 14286,
                                                 1312, 1340, 6276, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16386, 0, 3, 12612, 12702, 5247, 14412,
                                                 1340, 1368, 6360, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16554, 0, 3, 12702, 12792, 5310, 14538,
                                                 1368, 1396, 6444, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16722, 0, 3, 12972, 13062, 5436, 14664,
                                                 1452, 1480, 6528, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 16890, 0, 3, 13062, 13152, 5499, 14790,
                                                 1480, 1508, 6612, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17058, 0, 3, 13152, 13242, 5562, 14916,
                                                 1508, 1536, 6696, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17226, 0, 3, 13242, 13332, 5625, 15042,
                                                 1536, 1564, 6780, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17394, 0, 3, 13332, 13422, 5688, 15168,
                                                 1564, 1592, 6864, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17562, 0, 3, 13422, 13512, 5751, 15294,
                                                 1592, 1620, 6948, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 17730, 0, 3, 13512, 13602, 5814, 15420,
                                                 1620, 1648, 7032, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 17898, 0, 3, 13782, 13908, 6024, 15882,
                                                 1704, 1740, 7224, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18114, 0, 3, 13908, 14034, 6108, 16050,
                                                 1740, 1776, 7332, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18330, 0, 3, 14034, 14160, 6192, 16218,
                                                 1776, 1812, 7440, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18546, 0, 3, 14160, 14286, 6276, 16386,
                                                 1812, 1848, 7548, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18762, 0, 3, 14286, 14412, 6360, 16554,
                                                 1848, 1884, 7656, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 18978, 0, 3, 14664, 14790, 6612, 17058,
                                                 1956, 1992, 7872, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19194, 0, 3, 14790, 14916, 6696, 17226,
                                                 1992, 2028, 7980, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19410, 0, 3, 14916, 15042, 6780, 17394,
                                                 2028, 2064, 8088, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19626, 0, 3, 15042, 15168, 6864, 17562,
                                                 2064, 2100, 8196, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 19842, 0, 3, 15168, 15294, 6948, 17730,
                                                 2100, 2136, 8304, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 20058, 0, 3, 15546, 15714, 7116, 17898,
                                                 2208, 2253, 8412, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 20328, 0, 3, 15714, 15882, 7224, 18114,
                                                 2253, 2298, 8547, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 20598, 0, 3, 15882, 16050, 7332, 18330,
                                                 2298, 2343, 8682, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 20868, 0, 3, 16050, 16218, 7440, 18546,
                                                 2343, 2388, 8817, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 21138, 0, 3, 16218, 16386, 7548, 18762,
                                                 2388, 2433, 8952, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 21408, 0, 3, 16722, 16890, 7764, 18978,
                                                 2523, 2568, 9087, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 21678, 0, 3, 16890, 17058, 7872, 19194,
                                                 2568, 2613, 9222, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 21948, 0, 3, 17058, 17226, 7980, 19410,
                                                 2613, 2658, 9357, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 22218, 0, 3, 17226, 17394, 8088, 19626,
                                                 2658, 2703, 9492, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 22488, 0, 3, 17394, 17562, 8196, 19842,
                                                 2703, 2748, 9627, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22758, 3, 2838, 2841, 9768, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22768, 3, 2841, 2844, 9774, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22778, 3, 2844, 2847, 9780, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22788, 3, 2847, 2850, 9786, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22798, 3, 2850, 2853, 9792, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22808, 3, 2853, 2856, 9798, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22818, 3, 2856, 2859, 9804, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22828, 3, 2859, 2862, 9810, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22838, 3, 2862, 2865, 9816, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22848, 3, 2865, 2868, 9822, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22858, 3, 2874, 2877, 9834, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22868, 3, 2877, 2880, 9840, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22878, 3, 2880, 2883, 9846, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22888, 3, 2883, 2886, 9852, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22898, 3, 2886, 2889, 9858, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22908, 3, 2889, 2892, 9864, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22918, 3, 2892, 2895, 9870, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22928, 3, 2895, 2898, 9876, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22938, 3, 2898, 2901, 9882, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 22948, 3, 2901, 2904, 9888, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 22958, 0, 3, 9762, 22758, 9912, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 22988, 0, 3, 9768, 22768, 9930, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23018, 0, 3, 9774, 22778, 9948, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23048, 0, 3, 9780, 22788, 9966, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23078, 0, 3, 9786, 22798, 9984, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23108, 0, 3, 9792, 22808, 10002, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23138, 0, 3, 9798, 22818, 10020, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23168, 0, 3, 9804, 22828, 10038, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23198, 0, 3, 9810, 22838, 10056, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23228, 0, 3, 9816, 22848, 10074, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23258, 0, 3, 9828, 22858, 10110, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23288, 0, 3, 9834, 22868, 10128, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23318, 0, 3, 9840, 22878, 10146, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23348, 0, 3, 9846, 22888, 10164, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23378, 0, 3, 9852, 22898, 10182, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23408, 0, 3, 9858, 22908, 10200, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23438, 0, 3, 9864, 22918, 10218, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23468, 0, 3, 9870, 22928, 10236, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23498, 0, 3, 9876, 22938, 10254, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 23528, 0, 3, 9882, 22948, 10272, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 23558, 0, 3, 9912, 22988, 3126, 3144,
                                                 10362, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23618, 0, 3, 9930, 23018, 3144, 3162,
                                                 10398, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23678, 0, 3, 9948, 23048, 3162, 3180,
                                                 10434, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23738, 0, 3, 9966, 23078, 3180, 3198,
                                                 10470, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23798, 0, 3, 9984, 23108, 3198, 3216,
                                                 10506, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23858, 0, 3, 10002, 23138, 3216, 3234,
                                                 10542, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23918, 0, 3, 10020, 23168, 3234, 3252,
                                                 10578, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 23978, 0, 3, 10038, 23198, 3252, 3270,
                                                 10614, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24038, 0, 3, 10056, 23228, 3270, 3288,
                                                 10650, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24098, 0, 3, 10110, 23288, 3324, 3342,
                                                 10758, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24158, 0, 3, 10128, 23318, 3342, 3360,
                                                 10794, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24218, 0, 3, 10146, 23348, 3360, 3378,
                                                 10830, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24278, 0, 3, 10164, 23378, 3378, 3396,
                                                 10866, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24338, 0, 3, 10182, 23408, 3396, 3414,
                                                 10902, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24398, 0, 3, 10200, 23438, 3414, 3432,
                                                 10938, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24458, 0, 3, 10218, 23468, 3432, 3450,
                                                 10974, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24518, 0, 3, 10236, 23498, 3450, 3468,
                                                 11010, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 24578, 0, 3, 10254, 23528, 3468, 3486,
                                                 11046, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24638, 0, 3, 10362, 23618, 3522, 3552,
                                                 11142, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24738, 0, 3, 10398, 23678, 3552, 3582,
                                                 11202, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24838, 0, 3, 10434, 23738, 3582, 3612,
                                                 11262, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 24938, 0, 3, 10470, 23798, 3612, 3642,
                                                 11322, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25038, 0, 3, 10506, 23858, 3642, 3672,
                                                 11382, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25138, 0, 3, 10542, 23918, 3672, 3702,
                                                 11442, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25238, 0, 3, 10578, 23978, 3702, 3732,
                                                 11502, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25338, 0, 3, 10614, 24038, 3732, 3762,
                                                 11562, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25438, 0, 3, 10758, 24158, 3822, 3852,
                                                 11682, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25538, 0, 3, 10794, 24218, 3852, 3882,
                                                 11742, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25638, 0, 3, 10830, 24278, 3882, 3912,
                                                 11802, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25738, 0, 3, 10866, 24338, 3912, 3942,
                                                 11862, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25838, 0, 3, 10902, 24398, 3942, 3972,
                                                 11922, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 25938, 0, 3, 10938, 24458, 3972, 4002,
                                                 11982, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 26038, 0, 3, 10974, 24518, 4002, 4032,
                                                 12042, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 26138, 0, 3, 11010, 24578, 4032, 4062,
                                                 12102, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26238, 0, 3, 23558, 23618, 11142, 24738,
                                                 4122, 4167, 12342, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26388, 0, 3, 23618, 23678, 11202, 24838,
                                                 4167, 4212, 12432, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26538, 0, 3, 23678, 23738, 11262, 24938,
                                                 4212, 4257, 12522, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26688, 0, 3, 23738, 23798, 11322, 25038,
                                                 4257, 4302, 12612, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26838, 0, 3, 23798, 23858, 11382, 25138,
                                                 4302, 4347, 12702, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 26988, 0, 3, 23858, 23918, 11442, 25238,
                                                 4347, 4392, 12792, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 27138, 0, 3, 23918, 23978, 11502, 25338,
                                                 4392, 4437, 12882, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 27288, 0, 3, 24098, 24158, 11682, 25538,
                                                 4527, 4572, 13152, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 27438, 0, 3, 24158, 24218, 11742, 25638,
                                                 4572, 4617, 13242, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 27588, 0, 3, 24218, 24278, 11802, 25738,
                                                 4617, 4662, 13332, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 27738, 0, 3, 24278, 24338, 11862, 25838,
                                                 4662, 4707, 13422, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 27888, 0, 3, 24338, 24398, 11922, 25938,
                                                 4707, 4752, 13512, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 28038, 0, 3, 24398, 24458, 11982, 26038,
                                                 4752, 4797, 13602, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 28188, 0, 3, 24458, 24518, 12042, 26138,
                                                 4797, 4842, 13692, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28338, 0, 3, 24638, 24738, 12342, 26388,
                                                 4932, 4995, 13908, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28548, 0, 3, 24738, 24838, 12432, 26538,
                                                 4995, 5058, 14034, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28758, 0, 3, 24838, 24938, 12522, 26688,
                                                 5058, 5121, 14160, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 28968, 0, 3, 24938, 25038, 12612, 26838,
                                                 5121, 5184, 14286, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29178, 0, 3, 25038, 25138, 12702, 26988,
                                                 5184, 5247, 14412, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29388, 0, 3, 25138, 25238, 12792, 27138,
                                                 5247, 5310, 14538, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29598, 0, 3, 25438, 25538, 13152, 27438,
                                                 5436, 5499, 14790, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 29808, 0, 3, 25538, 25638, 13242, 27588,
                                                 5499, 5562, 14916, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 30018, 0, 3, 25638, 25738, 13332, 27738,
                                                 5562, 5625, 15042, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 30228, 0, 3, 25738, 25838, 13422, 27888,
                                                 5625, 5688, 15168, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 30438, 0, 3, 25838, 25938, 13512, 28038,
                                                 5688, 5751, 15294, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 30648, 0, 3, 25938, 26038, 13602, 28188,
                                                 5751, 5814, 15420, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 30858, 0, 3, 26238, 26388, 13908, 28548,
                                                 5940, 6024, 15882, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31138, 0, 3, 26388, 26538, 14034, 28758,
                                                 6024, 6108, 16050, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31418, 0, 3, 26538, 26688, 14160, 28968,
                                                 6108, 6192, 16218, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31698, 0, 3, 26688, 26838, 14286, 29178,
                                                 6192, 6276, 16386, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 31978, 0, 3, 26838, 26988, 14412, 29388,
                                                 6276, 6360, 16554, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 32258, 0, 3, 27288, 27438, 14790, 29808,
                                                 6528, 6612, 17058, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 32538, 0, 3, 27438, 27588, 14916, 30018,
                                                 6612, 6696, 17226, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 32818, 0, 3, 27588, 27738, 15042, 30228,
                                                 6696, 6780, 17394, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 33098, 0, 3, 27738, 27888, 15168, 30438,
                                                 6780, 6864, 17562, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 33378, 0, 3, 27888, 28038, 15294, 30648,
                                                 6864, 6948, 17730, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 33658, 0, 3, 28338, 28548, 15882, 31138,
                                                 7116, 7224, 18114, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 34018, 0, 3, 28548, 28758, 16050, 31418,
                                                 7224, 7332, 18330, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 34378, 0, 3, 28758, 28968, 16218, 31698,
                                                 7332, 7440, 18546, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 34738, 0, 3, 28968, 29178, 16386, 31978,
                                                 7440, 7548, 18762, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 35098, 0, 3, 29598, 29808, 17058, 32538,
                                                 7764, 7872, 19194, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 35458, 0, 3, 29808, 30018, 17226, 32818,
                                                 7872, 7980, 19410, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 35818, 0, 3, 30018, 30228, 17394, 33098,
                                                 7980, 8088, 19626, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 36178, 0, 3, 30228, 30438, 17562, 33378,
                                                 8088, 8196, 19842, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 36538, 0, 3, 30858, 31138, 18114, 34018,
                                                 8412, 8547, 20598, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 36988, 0, 3, 31138, 31418, 18330, 34378,
                                                 8547, 8682, 20868, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 37438, 0, 3, 31418, 31698, 18546, 34738,
                                                 8682, 8817, 21138, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 37888, 0, 3, 32258, 32538, 19194, 35458,
                                                 9087, 9222, 21948, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 38338, 0, 3, 32538, 32818, 19410, 35818,
                                                 9222, 9357, 22218, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 38788, 0, 3, 32818, 33098, 19626, 36178,
                                                 9357, 9492, 22488, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39238, 3, 9762, 9768, 22768, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39253, 3, 9768, 9774, 22778, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39268, 3, 9774, 9780, 22788, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39283, 3, 9780, 9786, 22798, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39298, 3, 9786, 9792, 22808, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39313, 3, 9792, 9798, 22818, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39328, 3, 9798, 9804, 22828, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39343, 3, 9804, 9810, 22838, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39358, 3, 9810, 9816, 22848, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39373, 3, 9828, 9834, 22868, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39388, 3, 9834, 9840, 22878, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39403, 3, 9840, 9846, 22888, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39418, 3, 9846, 9852, 22898, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39433, 3, 9852, 9858, 22908, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39448, 3, 9858, 9864, 22918, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39463, 3, 9864, 9870, 22928, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39478, 3, 9870, 9876, 22938, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 39493, 3, 9876, 9882, 22948, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 39508, 0, 3, 22758, 39238, 22988, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39553, 0, 3, 22768, 39253, 23018, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39598, 0, 3, 22778, 39268, 23048, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39643, 0, 3, 22788, 39283, 23078, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39688, 0, 3, 22798, 39298, 23108, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39733, 0, 3, 22808, 39313, 23138, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39778, 0, 3, 22818, 39328, 23168, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39823, 0, 3, 22828, 39343, 23198, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39868, 0, 3, 22838, 39358, 23228, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39913, 0, 3, 22858, 39373, 23288, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 39958, 0, 3, 22868, 39388, 23318, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 40003, 0, 3, 22878, 39403, 23348, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 40048, 0, 3, 22888, 39418, 23378, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 40093, 0, 3, 22898, 39433, 23408, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 40138, 0, 3, 22908, 39448, 23438, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 40183, 0, 3, 22918, 39463, 23468, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 40228, 0, 3, 22928, 39478, 23498, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 40273, 0, 3, 22938, 39493, 23528, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 40318, 0, 3, 22958, 39508, 10290, 10326,
                                                 23558, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 40408, 0, 3, 22988, 39553, 10326, 10362,
                                                 23618, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 40498, 0, 3, 23018, 39598, 10362, 10398,
                                                 23678, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 40588, 0, 3, 23048, 39643, 10398, 10434,
                                                 23738, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 40678, 0, 3, 23078, 39688, 10434, 10470,
                                                 23798, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 40768, 0, 3, 23108, 39733, 10470, 10506,
                                                 23858, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 40858, 0, 3, 23138, 39778, 10506, 10542,
                                                 23918, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 40948, 0, 3, 23168, 39823, 10542, 10578,
                                                 23978, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41038, 0, 3, 23198, 39868, 10578, 10614,
                                                 24038, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41128, 0, 3, 23258, 39913, 10686, 10722,
                                                 24098, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41218, 0, 3, 23288, 39958, 10722, 10758,
                                                 24158, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41308, 0, 3, 23318, 40003, 10758, 10794,
                                                 24218, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41398, 0, 3, 23348, 40048, 10794, 10830,
                                                 24278, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41488, 0, 3, 23378, 40093, 10830, 10866,
                                                 24338, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41578, 0, 3, 23408, 40138, 10866, 10902,
                                                 24398, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41668, 0, 3, 23438, 40183, 10902, 10938,
                                                 24458, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41758, 0, 3, 23468, 40228, 10938, 10974,
                                                 24518, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 41848, 0, 3, 23498, 40273, 10974, 11010,
                                                 24578, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 41938, 0, 3, 23618, 40498, 11082, 11142,
                                                 24738, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 42088, 0, 3, 23678, 40588, 11142, 11202,
                                                 24838, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 42238, 0, 3, 23738, 40678, 11202, 11262,
                                                 24938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 42388, 0, 3, 23798, 40768, 11262, 11322,
                                                 25038, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 42538, 0, 3, 23858, 40858, 11322, 11382,
                                                 25138, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 42688, 0, 3, 23918, 40948, 11382, 11442,
                                                 25238, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 42838, 0, 3, 23978, 41038, 11442, 11502,
                                                 25338, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 42988, 0, 3, 24158, 41308, 11622, 11682,
                                                 25538, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 43138, 0, 3, 24218, 41398, 11682, 11742,
                                                 25638, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 43288, 0, 3, 24278, 41488, 11742, 11802,
                                                 25738, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 43438, 0, 3, 24338, 41578, 11802, 11862,
                                                 25838, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 43588, 0, 3, 24398, 41668, 11862, 11922,
                                                 25938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 43738, 0, 3, 24458, 41758, 11922, 11982,
                                                 26038, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 43888, 0, 3, 24518, 41848, 11982, 12042,
                                                 26138, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44038, 0, 3, 40318, 40408, 24638, 41938,
                                                 12162, 12252, 26238, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44263, 0, 3, 40408, 40498, 24738, 42088,
                                                 12252, 12342, 26388, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44488, 0, 3, 40498, 40588, 24838, 42238,
                                                 12342, 12432, 26538, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44713, 0, 3, 40588, 40678, 24938, 42388,
                                                 12432, 12522, 26688, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 44938, 0, 3, 40678, 40768, 25038, 42538,
                                                 12522, 12612, 26838, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 45163, 0, 3, 40768, 40858, 25138, 42688,
                                                 12612, 12702, 26988, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 45388, 0, 3, 40858, 40948, 25238, 42838,
                                                 12702, 12792, 27138, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 45613, 0, 3, 41128, 41218, 25438, 42988,
                                                 12972, 13062, 27288, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 45838, 0, 3, 41218, 41308, 25538, 43138,
                                                 13062, 13152, 27438, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 46063, 0, 3, 41308, 41398, 25638, 43288,
                                                 13152, 13242, 27588, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 46288, 0, 3, 41398, 41488, 25738, 43438,
                                                 13242, 13332, 27738, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 46513, 0, 3, 41488, 41578, 25838, 43588,
                                                 13332, 13422, 27888, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 46738, 0, 3, 41578, 41668, 25938, 43738,
                                                 13422, 13512, 28038, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 46963, 0, 3, 41668, 41758, 26038, 43888,
                                                 13512, 13602, 28188, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 47188, 0, 3, 41938, 42088, 26388, 44488,
                                                 13782, 13908, 28548, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 47503, 0, 3, 42088, 42238, 26538, 44713,
                                                 13908, 14034, 28758, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 47818, 0, 3, 42238, 42388, 26688, 44938,
                                                 14034, 14160, 28968, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 48133, 0, 3, 42388, 42538, 26838, 45163,
                                                 14160, 14286, 29178, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 48448, 0, 3, 42538, 42688, 26988, 45388,
                                                 14286, 14412, 29388, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 48763, 0, 3, 42988, 43138, 27438, 46063,
                                                 14664, 14790, 29808, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 49078, 0, 3, 43138, 43288, 27588, 46288,
                                                 14790, 14916, 30018, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 49393, 0, 3, 43288, 43438, 27738, 46513,
                                                 14916, 15042, 30228, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 49708, 0, 3, 43438, 43588, 27888, 46738,
                                                 15042, 15168, 30438, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 50023, 0, 3, 43588, 43738, 28038, 46963,
                                                 15168, 15294, 30648, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 50338, 0, 3, 44038, 44263, 28338, 47188,
                                                 15546, 15714, 30858, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 50758, 0, 3, 44263, 44488, 28548, 47503,
                                                 15714, 15882, 31138, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 51178, 0, 3, 44488, 44713, 28758, 47818,
                                                 15882, 16050, 31418, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 51598, 0, 3, 44713, 44938, 28968, 48133,
                                                 16050, 16218, 31698, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 52018, 0, 3, 44938, 45163, 29178, 48448,
                                                 16218, 16386, 31978, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 52438, 0, 3, 45613, 45838, 29598, 48763,
                                                 16722, 16890, 32258, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 52858, 0, 3, 45838, 46063, 29808, 49078,
                                                 16890, 17058, 32538, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 53278, 0, 3, 46063, 46288, 30018, 49393,
                                                 17058, 17226, 32818, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 53698, 0, 3, 46288, 46513, 30228, 49708,
                                                 17226, 17394, 33098, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 54118, 0, 3, 46513, 46738, 30438, 50023,
                                                 17394, 17562, 33378, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 54538, 0, 3, 47188, 47503, 31138, 51178,
                                                 17898, 18114, 34018, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 55078, 0, 3, 47503, 47818, 31418, 51598,
                                                 18114, 18330, 34378, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 55618, 0, 3, 47818, 48133, 31698, 52018,
                                                 18330, 18546, 34738, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 56158, 0, 3, 48763, 49078, 32538, 53278,
                                                 18978, 19194, 35458, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 56698, 0, 3, 49078, 49393, 32818, 53698,
                                                 19194, 19410, 35818, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 57238, 0, 3, 49393, 49708, 33098, 54118,
                                                 19410, 19626, 36178, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 57778, 0, 3, 50338, 50758, 33658, 54538,
                                                 20058, 20328, 36538, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 58453, 0, 3, 50758, 51178, 34018, 55078,
                                                 20328, 20598, 36988, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 59128, 0, 3, 51178, 51598, 34378, 55618,
                                                 20598, 20868, 37438, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 59803, 0, 3, 52438, 52858, 35098, 56158,
                                                 21408, 21678, 37888, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 60478, 0, 3, 52858, 53278, 35458, 56698,
                                                 21678, 21948, 38338, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 61153, 0, 3, 53278, 53698, 35818, 57238,
                                                 21948, 22218, 38788, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61828, 3, 22758, 22768, 39253, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61849, 3, 22768, 22778, 39268, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61870, 3, 22778, 22788, 39283, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61891, 3, 22788, 22798, 39298, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61912, 3, 22798, 22808, 39313, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61933, 3, 22808, 22818, 39328, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61954, 3, 22818, 22828, 39343, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61975, 3, 22828, 22838, 39358, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 61996, 3, 22858, 22868, 39388, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 62017, 3, 22868, 22878, 39403, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 62038, 3, 22878, 22888, 39418, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 62059, 3, 22888, 22898, 39433, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 62080, 3, 22898, 22908, 39448, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 62101, 3, 22908, 22918, 39463, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 62122, 3, 22918, 22928, 39478, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 62143, 3, 22928, 22938, 39493, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 62164, 0, 3, 39238, 61828, 39553, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62227, 0, 3, 39253, 61849, 39598, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62290, 0, 3, 39268, 61870, 39643, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62353, 0, 3, 39283, 61891, 39688, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62416, 0, 3, 39298, 61912, 39733, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62479, 0, 3, 39313, 61933, 39778, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62542, 0, 3, 39328, 61954, 39823, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62605, 0, 3, 39343, 61975, 39868, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62668, 0, 3, 39373, 61996, 39958, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62731, 0, 3, 39388, 62017, 40003, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62794, 0, 3, 39403, 62038, 40048, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62857, 0, 3, 39418, 62059, 40093, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62920, 0, 3, 39433, 62080, 40138, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 62983, 0, 3, 39448, 62101, 40183, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 63046, 0, 3, 39463, 62122, 40228, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 63109, 0, 3, 39478, 62143, 40273, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 63172, 0, 3, 39553, 62227, 23558, 23618,
                                                 40498, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 63298, 0, 3, 39598, 62290, 23618, 23678,
                                                 40588, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 63424, 0, 3, 39643, 62353, 23678, 23738,
                                                 40678, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 63550, 0, 3, 39688, 62416, 23738, 23798,
                                                 40768, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 63676, 0, 3, 39733, 62479, 23798, 23858,
                                                 40858, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 63802, 0, 3, 39778, 62542, 23858, 23918,
                                                 40948, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 63928, 0, 3, 39823, 62605, 23918, 23978,
                                                 41038, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 64054, 0, 3, 39958, 62731, 24098, 24158,
                                                 41308, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 64180, 0, 3, 40003, 62794, 24158, 24218,
                                                 41398, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 64306, 0, 3, 40048, 62857, 24218, 24278,
                                                 41488, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 64432, 0, 3, 40093, 62920, 24278, 24338,
                                                 41578, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 64558, 0, 3, 40138, 62983, 24338, 24398,
                                                 41668, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 64684, 0, 3, 40183, 63046, 24398, 24458,
                                                 41758, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 64810, 0, 3, 40228, 63109, 24458, 24518,
                                                 41848, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 64936, 0, 3, 40498, 63298, 24638, 24738,
                                                 42088, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 65146, 0, 3, 40588, 63424, 24738, 24838,
                                                 42238, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 65356, 0, 3, 40678, 63550, 24838, 24938,
                                                 42388, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 65566, 0, 3, 40768, 63676, 24938, 25038,
                                                 42538, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 65776, 0, 3, 40858, 63802, 25038, 25138,
                                                 42688, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 65986, 0, 3, 40948, 63928, 25138, 25238,
                                                 42838, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 66196, 0, 3, 41308, 64180, 25438, 25538,
                                                 43138, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 66406, 0, 3, 41398, 64306, 25538, 25638,
                                                 43288, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 66616, 0, 3, 41488, 64432, 25638, 25738,
                                                 43438, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 66826, 0, 3, 41578, 64558, 25738, 25838,
                                                 43588, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 67036, 0, 3, 41668, 64684, 25838, 25938,
                                                 43738, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 67246, 0, 3, 41758, 64810, 25938, 26038,
                                                 43888, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 67456, 0, 3, 63172, 63298, 42088, 65146,
                                                 26238, 26388, 44488, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 67771, 0, 3, 63298, 63424, 42238, 65356,
                                                 26388, 26538, 44713, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 68086, 0, 3, 63424, 63550, 42388, 65566,
                                                 26538, 26688, 44938, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 68401, 0, 3, 63550, 63676, 42538, 65776,
                                                 26688, 26838, 45163, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 68716, 0, 3, 63676, 63802, 42688, 65986,
                                                 26838, 26988, 45388, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 69031, 0, 3, 64054, 64180, 43138, 66406,
                                                 27288, 27438, 46063, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 69346, 0, 3, 64180, 64306, 43288, 66616,
                                                 27438, 27588, 46288, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 69661, 0, 3, 64306, 64432, 43438, 66826,
                                                 27588, 27738, 46513, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 69976, 0, 3, 64432, 64558, 43588, 67036,
                                                 27738, 27888, 46738, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 70291, 0, 3, 64558, 64684, 43738, 67246,
                                                 27888, 28038, 46963, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 70606, 0, 3, 64936, 65146, 44488, 67771,
                                                 28338, 28548, 47503, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 71047, 0, 3, 65146, 65356, 44713, 68086,
                                                 28548, 28758, 47818, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 71488, 0, 3, 65356, 65566, 44938, 68401,
                                                 28758, 28968, 48133, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 71929, 0, 3, 65566, 65776, 45163, 68716,
                                                 28968, 29178, 48448, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 72370, 0, 3, 66196, 66406, 46063, 69346,
                                                 29598, 29808, 49078, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 72811, 0, 3, 66406, 66616, 46288, 69661,
                                                 29808, 30018, 49393, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 73252, 0, 3, 66616, 66826, 46513, 69976,
                                                 30018, 30228, 49708, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 73693, 0, 3, 66826, 67036, 46738, 70291,
                                                 30228, 30438, 50023, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 74134, 0, 3, 67456, 67771, 47503, 71047,
                                                 30858, 31138, 51178, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 74722, 0, 3, 67771, 68086, 47818, 71488,
                                                 31138, 31418, 51598, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 75310, 0, 3, 68086, 68401, 48133, 71929,
                                                 31418, 31698, 52018, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 75898, 0, 3, 69031, 69346, 49078, 72811,
                                                 32258, 32538, 53278, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 76486, 0, 3, 69346, 69661, 49393, 73252,
                                                 32538, 32818, 53698, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 77074, 0, 3, 69661, 69976, 49708, 73693,
                                                 32818, 33098, 54118, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 77662, 0, 3, 70606, 71047, 51178, 74722,
                                                 33658, 34018, 55078, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 78418, 0, 3, 71047, 71488, 51598, 75310,
                                                 34018, 34378, 55618, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 79174, 0, 3, 72370, 72811, 53278, 76486,
                                                 35098, 35458, 56698, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 79930, 0, 3, 72811, 73252, 53698, 77074,
                                                 35458, 35818, 57238, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 80686, 0, 3, 74134, 74722, 55078, 78418,
                                                 36538, 36988, 59128, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 81631, 0, 3, 75898, 76486, 56698, 79930,
                                                 37888, 38338, 61153, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82576, 3, 39238, 39253, 61849, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82604, 3, 39253, 39268, 61870, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82632, 3, 39268, 39283, 61891, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82660, 3, 39283, 39298, 61912, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82688, 3, 39298, 39313, 61933, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82716, 3, 39313, 39328, 61954, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82744, 3, 39328, 39343, 61975, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82772, 3, 39373, 39388, 62017, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82800, 3, 39388, 39403, 62038, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82828, 3, 39403, 39418, 62059, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82856, 3, 39418, 39433, 62080, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82884, 3, 39433, 39448, 62101, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82912, 3, 39448, 39463, 62122, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 82940, 3, 39463, 39478, 62143, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 82968, 0, 3, 61828, 82576, 62227, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83052, 0, 3, 61849, 82604, 62290, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83136, 0, 3, 61870, 82632, 62353, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83220, 0, 3, 61891, 82660, 62416, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83304, 0, 3, 61912, 82688, 62479, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83388, 0, 3, 61933, 82716, 62542, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83472, 0, 3, 61954, 82744, 62605, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83556, 0, 3, 61996, 82772, 62731, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83640, 0, 3, 62017, 82800, 62794, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83724, 0, 3, 62038, 82828, 62857, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83808, 0, 3, 62059, 82856, 62920, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83892, 0, 3, 62080, 82884, 62983, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 83976, 0, 3, 62101, 82912, 63046, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 84060, 0, 3, 62122, 82940, 63109, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 84144, 0, 3, 62164, 82968, 40318, 40408,
                                                 63172, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84312, 0, 3, 62227, 83052, 40408, 40498,
                                                 63298, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84480, 0, 3, 62290, 83136, 40498, 40588,
                                                 63424, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84648, 0, 3, 62353, 83220, 40588, 40678,
                                                 63550, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84816, 0, 3, 62416, 83304, 40678, 40768,
                                                 63676, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 84984, 0, 3, 62479, 83388, 40768, 40858,
                                                 63802, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85152, 0, 3, 62542, 83472, 40858, 40948,
                                                 63928, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85320, 0, 3, 62668, 83556, 41128, 41218,
                                                 64054, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85488, 0, 3, 62731, 83640, 41218, 41308,
                                                 64180, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85656, 0, 3, 62794, 83724, 41308, 41398,
                                                 64306, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85824, 0, 3, 62857, 83808, 41398, 41488,
                                                 64432, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 85992, 0, 3, 62920, 83892, 41488, 41578,
                                                 64558, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 86160, 0, 3, 62983, 83976, 41578, 41668,
                                                 64684, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 86328, 0, 3, 63046, 84060, 41668, 41758,
                                                 64810, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 86496, 0, 3, 63298, 84480, 41938, 42088,
                                                 65146, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 86776, 0, 3, 63424, 84648, 42088, 42238,
                                                 65356, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 87056, 0, 3, 63550, 84816, 42238, 42388,
                                                 65566, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 87336, 0, 3, 63676, 84984, 42388, 42538,
                                                 65776, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 87616, 0, 3, 63802, 85152, 42538, 42688,
                                                 65986, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 87896, 0, 3, 64180, 85656, 42988, 43138,
                                                 66406, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 88176, 0, 3, 64306, 85824, 43138, 43288,
                                                 66616, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 88456, 0, 3, 64432, 85992, 43288, 43438,
                                                 66826, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 88736, 0, 3, 64558, 86160, 43438, 43588,
                                                 67036, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 89016, 0, 3, 64684, 86328, 43588, 43738,
                                                 67246, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 89296, 0, 3, 84144, 84312, 64936, 86496,
                                                 44038, 44263, 67456, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 89716, 0, 3, 84312, 84480, 65146, 86776,
                                                 44263, 44488, 67771, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 90136, 0, 3, 84480, 84648, 65356, 87056,
                                                 44488, 44713, 68086, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 90556, 0, 3, 84648, 84816, 65566, 87336,
                                                 44713, 44938, 68401, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 90976, 0, 3, 84816, 84984, 65776, 87616,
                                                 44938, 45163, 68716, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 91396, 0, 3, 85320, 85488, 66196, 87896,
                                                 45613, 45838, 69031, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 91816, 0, 3, 85488, 85656, 66406, 88176,
                                                 45838, 46063, 69346, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 92236, 0, 3, 85656, 85824, 66616, 88456,
                                                 46063, 46288, 69661, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 92656, 0, 3, 85824, 85992, 66826, 88736,
                                                 46288, 46513, 69976, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 93076, 0, 3, 85992, 86160, 67036, 89016,
                                                 46513, 46738, 70291, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 93496, 0, 3, 86496, 86776, 67771, 90136,
                                                 47188, 47503, 71047, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 94084, 0, 3, 86776, 87056, 68086, 90556,
                                                 47503, 47818, 71488, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 94672, 0, 3, 87056, 87336, 68401, 90976,
                                                 47818, 48133, 71929, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 95260, 0, 3, 87896, 88176, 69346, 92236,
                                                 48763, 49078, 72811, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 95848, 0, 3, 88176, 88456, 69661, 92656,
                                                 49078, 49393, 73252, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 96436, 0, 3, 88456, 88736, 69976, 93076,
                                                 49393, 49708, 73693, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 97024, 0, 3, 89296, 89716, 70606, 93496,
                                                 50338, 50758, 74134, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 97808, 0, 3, 89716, 90136, 71047, 94084,
                                                 50758, 51178, 74722, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 98592, 0, 3, 90136, 90556, 71488, 94672,
                                                 51178, 51598, 75310, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 99376, 0, 3, 91396, 91816, 72370, 95260,
                                                 52438, 52858, 75898, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 100160, 0, 3, 91816, 92236, 72811,
                                                 95848, 52858, 53278, 76486, ncols, alpha, beta,
                                                 p);

            compute_prim_ii_electron_repulsion_0(buffer, 100944, 0, 3, 92236, 92656, 73252,
                                                 96436, 53278, 53698, 77074, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 101728, 0, 3, 93496, 94084, 74722,
                                                 98592, 54538, 55078, 78418, ncols, alpha, beta,
                                                 p);

            compute_prim_ki_electron_repulsion_0(buffer, 102736, 0, 3, 95260, 95848, 76486,
                                                 100944, 56158, 56698, 79930, ncols, alpha, beta,
                                                 p);

            compute_prim_li_electron_repulsion_0(buffer, 103744, 0, 3, 97024, 97808, 77662,
                                                 101728, 57778, 58453, 80686, ncols, alpha, beta,
                                                 p);

            compute_prim_li_electron_repulsion_0(buffer, 105004, 0, 3, 99376, 100160, 79174,
                                                 102736, 59803, 60478, 81631, ncols, alpha, beta,
                                                 p);

            simdgeo::geom_k_x(buffer, 106264, 99376, 105004, 1, 28, ncols, alpha);

            simdgeo::geom_k_y(buffer, 107272, 99376, 105004, 1, 28, ncols, alpha);

            simdgeo::geom_k_z(buffer, 108280, 99376, 105004, 1, 28, ncols, alpha);

            simdgeo::geom_k_x(buffer, 109288, 97024, 103744, 1, 28, ncols, alpha);

            simdgeo::geom_k_y(buffer, 110296, 97024, 103744, 1, 28, ncols, alpha);

            simdgeo::geom_k_z(buffer, 111304, 97024, 103744, 1, 28, ncols, alpha);

            simdfunc::contract_primitives(buffer, 112312, 109288, 3024, ncols);

            simdfunc::contract_primitives(buffer, 115336, 106264, 3024, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 118360, 115336, 36, 1, nmax);

    simdtrf::transform_k_outer(values, nvalues, buffer, 118360, 13, nmax);

    simdtrf::transform_i_inner(buffer, 118360, 116344, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 195 * nvalues, nvalues, buffer, 118360, 13, nmax);

    simdtrf::transform_i_inner(buffer, 118360, 117352, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 390 * nvalues, nvalues, buffer, 118360, 13, nmax);

    simdtrf::transform_i_inner(buffer, 118360, 112312, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 585 * nvalues, nvalues, buffer, 118360, 13, nmax);

    simdtrf::transform_i_inner(buffer, 118360, 113320, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 780 * nvalues, nvalues, buffer, 118360, 13, nmax);

    simdtrf::transform_i_inner(buffer, 118360, 114328, 36, 1, nmax);

    simdtrf::transform_k_outer(values + 975 * nvalues, nvalues, buffer, 118360, 13, nmax);
}

}  // namespace simdt2ceri
