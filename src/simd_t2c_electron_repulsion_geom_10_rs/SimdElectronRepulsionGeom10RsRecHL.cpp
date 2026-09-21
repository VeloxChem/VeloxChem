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


#include "SimdElectronRepulsionGeom10RsRecHL.hpp"

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
#include "SimdGeometryH1.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_hl_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_hl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 117709, 111682, 5670, nvalues);

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

            compute_prim_sp_electron_repulsion_0(buffer, 1704, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1707, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1710, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1713, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1716, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1719, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1722, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1725, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1728, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1731, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1734, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1737, 3, 21, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1740, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1743, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1746, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1749, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1752, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1755, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1758, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1761, 3, 33, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1764, 3, 34, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1767, 3, 35, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1770, 3, 36, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1773, 3, 37, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1776, 3, 9, 41, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1785, 3, 10, 44, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1794, 3, 11, 47, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1803, 3, 12, 50, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1812, 3, 13, 53, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1821, 3, 14, 56, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1830, 3, 15, 59, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1839, 3, 16, 62, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1848, 3, 17, 65, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1857, 3, 18, 68, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1866, 3, 19, 71, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1875, 3, 20, 74, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1884, 3, 25, 80, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1893, 3, 26, 83, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1902, 3, 27, 86, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1911, 3, 28, 89, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1920, 3, 29, 92, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1929, 3, 30, 95, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1938, 3, 31, 98, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1947, 3, 32, 101, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1956, 3, 33, 104, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1965, 3, 34, 107, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1974, 3, 35, 110, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1983, 3, 36, 113, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1992, 0, 3, 41, 1785, 128, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2010, 0, 3, 44, 1794, 134, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2028, 0, 3, 47, 1803, 140, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2046, 0, 3, 50, 1812, 146, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2064, 0, 3, 53, 1821, 152, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2082, 0, 3, 56, 1830, 158, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2100, 0, 3, 59, 1839, 164, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2118, 0, 3, 62, 1848, 170, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2136, 0, 3, 65, 1857, 176, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2154, 0, 3, 68, 1866, 182, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2172, 0, 3, 71, 1875, 188, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2190, 0, 3, 80, 1893, 206, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2208, 0, 3, 83, 1902, 212, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2226, 0, 3, 86, 1911, 218, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2244, 0, 3, 89, 1920, 224, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2262, 0, 3, 92, 1929, 230, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2280, 0, 3, 95, 1938, 236, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2298, 0, 3, 98, 1947, 242, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2316, 0, 3, 101, 1956, 248, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2334, 0, 3, 104, 1965, 254, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2352, 0, 3, 107, 1974, 260, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 2370, 0, 3, 110, 1983, 266, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2388, 0, 3, 128, 2010, 282, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2418, 0, 3, 134, 2028, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2448, 0, 3, 140, 2046, 302, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2478, 0, 3, 146, 2064, 312, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2508, 0, 3, 152, 2082, 322, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2538, 0, 3, 158, 2100, 332, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2568, 0, 3, 164, 2118, 342, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2598, 0, 3, 170, 2136, 352, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2628, 0, 3, 176, 2154, 362, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2658, 0, 3, 182, 2172, 372, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2688, 0, 3, 206, 2208, 392, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2718, 0, 3, 212, 2226, 402, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2748, 0, 3, 218, 2244, 412, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2778, 0, 3, 224, 2262, 422, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2808, 0, 3, 230, 2280, 432, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2838, 0, 3, 236, 2298, 442, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2868, 0, 3, 242, 2316, 452, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2898, 0, 3, 248, 2334, 462, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2928, 0, 3, 254, 2352, 472, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2958, 0, 3, 260, 2370, 482, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2988, 0, 3, 282, 2418, 522, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3033, 0, 3, 292, 2448, 537, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3078, 0, 3, 302, 2478, 552, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3123, 0, 3, 312, 2508, 567, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3168, 0, 3, 322, 2538, 582, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3213, 0, 3, 332, 2568, 597, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3258, 0, 3, 342, 2598, 612, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3303, 0, 3, 352, 2628, 627, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3348, 0, 3, 362, 2658, 642, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3393, 0, 3, 392, 2718, 687, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3438, 0, 3, 402, 2748, 702, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3483, 0, 3, 412, 2778, 717, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3528, 0, 3, 422, 2808, 732, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3573, 0, 3, 432, 2838, 747, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3618, 0, 3, 442, 2868, 762, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3663, 0, 3, 452, 2898, 777, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3708, 0, 3, 462, 2928, 792, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 3753, 0, 3, 472, 2958, 807, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3798, 0, 3, 522, 3033, 843, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3861, 0, 3, 537, 3078, 864, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3924, 0, 3, 552, 3123, 885, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3987, 0, 3, 567, 3168, 906, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4050, 0, 3, 582, 3213, 927, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4113, 0, 3, 597, 3258, 948, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4176, 0, 3, 612, 3303, 969, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4239, 0, 3, 627, 3348, 990, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4302, 0, 3, 687, 3438, 1032, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4365, 0, 3, 702, 3483, 1053, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4428, 0, 3, 717, 3528, 1074, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4491, 0, 3, 732, 3573, 1095, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4554, 0, 3, 747, 3618, 1116, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4617, 0, 3, 762, 3663, 1137, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4680, 0, 3, 777, 3708, 1158, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 4743, 0, 3, 792, 3753, 1179, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4806, 0, 3, 843, 3861, 1256, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4890, 0, 3, 864, 3924, 1284, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 4974, 0, 3, 885, 3987, 1312, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5058, 0, 3, 906, 4050, 1340, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5142, 0, 3, 927, 4113, 1368, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5226, 0, 3, 948, 4176, 1396, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5310, 0, 3, 969, 4239, 1424, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 5394, 0, 3, 1032, 4365, 1508, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5478, 0, 3, 1053, 4428, 1536, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5562, 0, 3, 1074, 4491, 1564, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5646, 0, 3, 1095, 4554, 1592, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5730, 0, 3, 1116, 4617, 1620, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5814, 0, 3, 1137, 4680, 1648, ncols,
                                                 p);

            compute_prim_ip_electron_repulsion_0(buffer, 5898, 0, 3, 1158, 4743, 1676, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 5982, 3, 9, 10, 1707, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5988, 3, 10, 11, 1710, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5994, 3, 11, 12, 1713, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6000, 3, 12, 13, 1716, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6006, 3, 13, 14, 1719, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6012, 3, 14, 15, 1722, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6018, 3, 15, 16, 1725, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6024, 3, 16, 17, 1728, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6030, 3, 17, 18, 1731, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6036, 3, 18, 19, 1734, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6042, 3, 19, 20, 1737, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6048, 3, 25, 26, 1743, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6054, 3, 26, 27, 1746, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6060, 3, 27, 28, 1749, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6066, 3, 28, 29, 1752, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6072, 3, 29, 30, 1755, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6078, 3, 30, 31, 1758, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6084, 3, 31, 32, 1761, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6090, 3, 32, 33, 1764, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6096, 3, 33, 34, 1767, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6102, 3, 34, 35, 1770, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 6108, 3, 35, 36, 1773, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 6114, 0, 3, 1704, 5982, 1785, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6132, 0, 3, 1707, 5988, 1794, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6150, 0, 3, 1710, 5994, 1803, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6168, 0, 3, 1713, 6000, 1812, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6186, 0, 3, 1716, 6006, 1821, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6204, 0, 3, 1719, 6012, 1830, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6222, 0, 3, 1722, 6018, 1839, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6240, 0, 3, 1725, 6024, 1848, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6258, 0, 3, 1728, 6030, 1857, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6276, 0, 3, 1731, 6036, 1866, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6294, 0, 3, 1734, 6042, 1875, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6312, 0, 3, 1740, 6048, 1893, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6330, 0, 3, 1743, 6054, 1902, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6348, 0, 3, 1746, 6060, 1911, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6366, 0, 3, 1749, 6066, 1920, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6384, 0, 3, 1752, 6072, 1929, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6402, 0, 3, 1755, 6078, 1938, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6420, 0, 3, 1758, 6084, 1947, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6438, 0, 3, 1761, 6090, 1956, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6456, 0, 3, 1764, 6096, 1965, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6474, 0, 3, 1767, 6102, 1974, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6492, 0, 3, 1770, 6108, 1983, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 6510, 0, 3, 1776, 6114, 116, 122, 1992,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6546, 0, 3, 1785, 6132, 122, 128, 2010,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6582, 0, 3, 1794, 6150, 128, 134, 2028,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6618, 0, 3, 1803, 6168, 134, 140, 2046,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6654, 0, 3, 1812, 6186, 140, 146, 2064,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6690, 0, 3, 1821, 6204, 146, 152, 2082,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6726, 0, 3, 1830, 6222, 152, 158, 2100,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6762, 0, 3, 1839, 6240, 158, 164, 2118,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6798, 0, 3, 1848, 6258, 164, 170, 2136,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6834, 0, 3, 1857, 6276, 170, 176, 2154,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6870, 0, 3, 1866, 6294, 176, 182, 2172,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6906, 0, 3, 1884, 6312, 194, 200, 2190,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6942, 0, 3, 1893, 6330, 200, 206, 2208,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6978, 0, 3, 1902, 6348, 206, 212, 2226,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7014, 0, 3, 1911, 6366, 212, 218, 2244,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7050, 0, 3, 1920, 6384, 218, 224, 2262,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7086, 0, 3, 1929, 6402, 224, 230, 2280,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7122, 0, 3, 1938, 6420, 230, 236, 2298,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7158, 0, 3, 1947, 6438, 236, 242, 2316,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7194, 0, 3, 1956, 6456, 242, 248, 2334,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7230, 0, 3, 1965, 6474, 248, 254, 2352,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 7266, 0, 3, 1974, 6492, 254, 260, 2370,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7302, 0, 3, 2010, 6582, 272, 282, 2418,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7362, 0, 3, 2028, 6618, 282, 292, 2448,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7422, 0, 3, 2046, 6654, 292, 302, 2478,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7482, 0, 3, 2064, 6690, 302, 312, 2508,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7542, 0, 3, 2082, 6726, 312, 322, 2538,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7602, 0, 3, 2100, 6762, 322, 332, 2568,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7662, 0, 3, 2118, 6798, 332, 342, 2598,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7722, 0, 3, 2136, 6834, 342, 352, 2628,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7782, 0, 3, 2154, 6870, 352, 362, 2658,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7842, 0, 3, 2208, 6978, 382, 392, 2718,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7902, 0, 3, 2226, 7014, 392, 402, 2748,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 7962, 0, 3, 2244, 7050, 402, 412, 2778,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8022, 0, 3, 2262, 7086, 412, 422, 2808,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8082, 0, 3, 2280, 7122, 422, 432, 2838,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8142, 0, 3, 2298, 7158, 432, 442, 2868,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8202, 0, 3, 2316, 7194, 442, 452, 2898,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8262, 0, 3, 2334, 7230, 452, 462, 2928,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 8322, 0, 3, 2352, 7266, 462, 472, 2958,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8382, 0, 3, 6510, 6546, 2388, 7302, 492,
                                                 507, 2988, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8472, 0, 3, 6546, 6582, 2418, 7362, 507,
                                                 522, 3033, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8562, 0, 3, 6582, 6618, 2448, 7422, 522,
                                                 537, 3078, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8652, 0, 3, 6618, 6654, 2478, 7482, 537,
                                                 552, 3123, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8742, 0, 3, 6654, 6690, 2508, 7542, 552,
                                                 567, 3168, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8832, 0, 3, 6690, 6726, 2538, 7602, 567,
                                                 582, 3213, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 8922, 0, 3, 6726, 6762, 2568, 7662, 582,
                                                 597, 3258, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9012, 0, 3, 6762, 6798, 2598, 7722, 597,
                                                 612, 3303, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9102, 0, 3, 6798, 6834, 2628, 7782, 612,
                                                 627, 3348, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9192, 0, 3, 6906, 6942, 2688, 7842, 657,
                                                 672, 3393, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9282, 0, 3, 6942, 6978, 2718, 7902, 672,
                                                 687, 3438, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9372, 0, 3, 6978, 7014, 2748, 7962, 687,
                                                 702, 3483, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9462, 0, 3, 7014, 7050, 2778, 8022, 702,
                                                 717, 3528, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9552, 0, 3, 7050, 7086, 2808, 8082, 717,
                                                 732, 3573, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9642, 0, 3, 7086, 7122, 2838, 8142, 732,
                                                 747, 3618, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9732, 0, 3, 7122, 7158, 2868, 8202, 747,
                                                 762, 3663, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9822, 0, 3, 7158, 7194, 2898, 8262, 762,
                                                 777, 3708, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 9912, 0, 3, 7194, 7230, 2928, 8322, 777,
                                                 792, 3753, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10002, 0, 3, 7302, 7362, 3033, 8562,
                                                 822, 843, 3861, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10128, 0, 3, 7362, 7422, 3078, 8652,
                                                 843, 864, 3924, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10254, 0, 3, 7422, 7482, 3123, 8742,
                                                 864, 885, 3987, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10380, 0, 3, 7482, 7542, 3168, 8832,
                                                 885, 906, 4050, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10506, 0, 3, 7542, 7602, 3213, 8922,
                                                 906, 927, 4113, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10632, 0, 3, 7602, 7662, 3258, 9012,
                                                 927, 948, 4176, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10758, 0, 3, 7662, 7722, 3303, 9102,
                                                 948, 969, 4239, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 10884, 0, 3, 7842, 7902, 3438, 9372,
                                                 1011, 1032, 4365, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11010, 0, 3, 7902, 7962, 3483, 9462,
                                                 1032, 1053, 4428, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11136, 0, 3, 7962, 8022, 3528, 9552,
                                                 1053, 1074, 4491, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11262, 0, 3, 8022, 8082, 3573, 9642,
                                                 1074, 1095, 4554, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11388, 0, 3, 8082, 8142, 3618, 9732,
                                                 1095, 1116, 4617, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11514, 0, 3, 8142, 8202, 3663, 9822,
                                                 1116, 1137, 4680, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 11640, 0, 3, 8202, 8262, 3708, 9912,
                                                 1137, 1158, 4743, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11766, 0, 3, 8382, 8472, 3798, 10002,
                                                 1200, 1228, 4806, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 11934, 0, 3, 8472, 8562, 3861, 10128,
                                                 1228, 1256, 4890, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12102, 0, 3, 8562, 8652, 3924, 10254,
                                                 1256, 1284, 4974, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12270, 0, 3, 8652, 8742, 3987, 10380,
                                                 1284, 1312, 5058, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12438, 0, 3, 8742, 8832, 4050, 10506,
                                                 1312, 1340, 5142, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12606, 0, 3, 8832, 8922, 4113, 10632,
                                                 1340, 1368, 5226, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12774, 0, 3, 8922, 9012, 4176, 10758,
                                                 1368, 1396, 5310, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 12942, 0, 3, 9192, 9282, 4302, 10884,
                                                 1452, 1480, 5394, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13110, 0, 3, 9282, 9372, 4365, 11010,
                                                 1480, 1508, 5478, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13278, 0, 3, 9372, 9462, 4428, 11136,
                                                 1508, 1536, 5562, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13446, 0, 3, 9462, 9552, 4491, 11262,
                                                 1536, 1564, 5646, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13614, 0, 3, 9552, 9642, 4554, 11388,
                                                 1564, 1592, 5730, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13782, 0, 3, 9642, 9732, 4617, 11514,
                                                 1592, 1620, 5814, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 13950, 0, 3, 9732, 9822, 4680, 11640,
                                                 1620, 1648, 5898, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14118, 3, 1704, 1707, 5988, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14128, 3, 1707, 1710, 5994, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14138, 3, 1710, 1713, 6000, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14148, 3, 1713, 1716, 6006, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14158, 3, 1716, 1719, 6012, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14168, 3, 1719, 1722, 6018, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14178, 3, 1722, 1725, 6024, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14188, 3, 1725, 1728, 6030, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14198, 3, 1728, 1731, 6036, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14208, 3, 1731, 1734, 6042, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14218, 3, 1740, 1743, 6054, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14228, 3, 1743, 1746, 6060, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14238, 3, 1746, 1749, 6066, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14248, 3, 1749, 1752, 6072, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14258, 3, 1752, 1755, 6078, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14268, 3, 1755, 1758, 6084, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14278, 3, 1758, 1761, 6090, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14288, 3, 1761, 1764, 6096, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14298, 3, 1764, 1767, 6102, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 14308, 3, 1767, 1770, 6108, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 14318, 0, 3, 5982, 14118, 6132, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14348, 0, 3, 5988, 14128, 6150, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14378, 0, 3, 5994, 14138, 6168, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14408, 0, 3, 6000, 14148, 6186, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14438, 0, 3, 6006, 14158, 6204, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14468, 0, 3, 6012, 14168, 6222, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14498, 0, 3, 6018, 14178, 6240, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14528, 0, 3, 6024, 14188, 6258, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14558, 0, 3, 6030, 14198, 6276, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14588, 0, 3, 6036, 14208, 6294, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14618, 0, 3, 6048, 14218, 6330, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14648, 0, 3, 6054, 14228, 6348, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14678, 0, 3, 6060, 14238, 6366, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14708, 0, 3, 6066, 14248, 6384, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14738, 0, 3, 6072, 14258, 6402, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14768, 0, 3, 6078, 14268, 6420, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14798, 0, 3, 6084, 14278, 6438, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14828, 0, 3, 6090, 14288, 6456, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14858, 0, 3, 6096, 14298, 6474, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 14888, 0, 3, 6102, 14308, 6492, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 14918, 0, 3, 6132, 14348, 1992, 2010,
                                                 6582, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 14978, 0, 3, 6150, 14378, 2010, 2028,
                                                 6618, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15038, 0, 3, 6168, 14408, 2028, 2046,
                                                 6654, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15098, 0, 3, 6186, 14438, 2046, 2064,
                                                 6690, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15158, 0, 3, 6204, 14468, 2064, 2082,
                                                 6726, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15218, 0, 3, 6222, 14498, 2082, 2100,
                                                 6762, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15278, 0, 3, 6240, 14528, 2100, 2118,
                                                 6798, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15338, 0, 3, 6258, 14558, 2118, 2136,
                                                 6834, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15398, 0, 3, 6276, 14588, 2136, 2154,
                                                 6870, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15458, 0, 3, 6330, 14648, 2190, 2208,
                                                 6978, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15518, 0, 3, 6348, 14678, 2208, 2226,
                                                 7014, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15578, 0, 3, 6366, 14708, 2226, 2244,
                                                 7050, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15638, 0, 3, 6384, 14738, 2244, 2262,
                                                 7086, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15698, 0, 3, 6402, 14768, 2262, 2280,
                                                 7122, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15758, 0, 3, 6420, 14798, 2280, 2298,
                                                 7158, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15818, 0, 3, 6438, 14828, 2298, 2316,
                                                 7194, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15878, 0, 3, 6456, 14858, 2316, 2334,
                                                 7230, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 15938, 0, 3, 6474, 14888, 2334, 2352,
                                                 7266, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 15998, 0, 3, 6582, 14978, 2388, 2418,
                                                 7362, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16098, 0, 3, 6618, 15038, 2418, 2448,
                                                 7422, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16198, 0, 3, 6654, 15098, 2448, 2478,
                                                 7482, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16298, 0, 3, 6690, 15158, 2478, 2508,
                                                 7542, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16398, 0, 3, 6726, 15218, 2508, 2538,
                                                 7602, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16498, 0, 3, 6762, 15278, 2538, 2568,
                                                 7662, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16598, 0, 3, 6798, 15338, 2568, 2598,
                                                 7722, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16698, 0, 3, 6834, 15398, 2598, 2628,
                                                 7782, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16798, 0, 3, 6978, 15518, 2688, 2718,
                                                 7902, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16898, 0, 3, 7014, 15578, 2718, 2748,
                                                 7962, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 16998, 0, 3, 7050, 15638, 2748, 2778,
                                                 8022, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17098, 0, 3, 7086, 15698, 2778, 2808,
                                                 8082, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17198, 0, 3, 7122, 15758, 2808, 2838,
                                                 8142, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17298, 0, 3, 7158, 15818, 2838, 2868,
                                                 8202, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17398, 0, 3, 7194, 15878, 2868, 2898,
                                                 8262, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 17498, 0, 3, 7230, 15938, 2898, 2928,
                                                 8322, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17598, 0, 3, 14918, 14978, 7362, 16098,
                                                 2988, 3033, 8562, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17748, 0, 3, 14978, 15038, 7422, 16198,
                                                 3033, 3078, 8652, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 17898, 0, 3, 15038, 15098, 7482, 16298,
                                                 3078, 3123, 8742, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18048, 0, 3, 15098, 15158, 7542, 16398,
                                                 3123, 3168, 8832, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18198, 0, 3, 15158, 15218, 7602, 16498,
                                                 3168, 3213, 8922, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18348, 0, 3, 15218, 15278, 7662, 16598,
                                                 3213, 3258, 9012, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18498, 0, 3, 15278, 15338, 7722, 16698,
                                                 3258, 3303, 9102, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18648, 0, 3, 15458, 15518, 7902, 16898,
                                                 3393, 3438, 9372, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18798, 0, 3, 15518, 15578, 7962, 16998,
                                                 3438, 3483, 9462, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 18948, 0, 3, 15578, 15638, 8022, 17098,
                                                 3483, 3528, 9552, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19098, 0, 3, 15638, 15698, 8082, 17198,
                                                 3528, 3573, 9642, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19248, 0, 3, 15698, 15758, 8142, 17298,
                                                 3573, 3618, 9732, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19398, 0, 3, 15758, 15818, 8202, 17398,
                                                 3618, 3663, 9822, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 19548, 0, 3, 15818, 15878, 8262, 17498,
                                                 3663, 3708, 9912, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19698, 0, 3, 15998, 16098, 8562, 17748,
                                                 3798, 3861, 10128, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 19908, 0, 3, 16098, 16198, 8652, 17898,
                                                 3861, 3924, 10254, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20118, 0, 3, 16198, 16298, 8742, 18048,
                                                 3924, 3987, 10380, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20328, 0, 3, 16298, 16398, 8832, 18198,
                                                 3987, 4050, 10506, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20538, 0, 3, 16398, 16498, 8922, 18348,
                                                 4050, 4113, 10632, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20748, 0, 3, 16498, 16598, 9012, 18498,
                                                 4113, 4176, 10758, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 20958, 0, 3, 16798, 16898, 9372, 18798,
                                                 4302, 4365, 11010, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21168, 0, 3, 16898, 16998, 9462, 18948,
                                                 4365, 4428, 11136, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21378, 0, 3, 16998, 17098, 9552, 19098,
                                                 4428, 4491, 11262, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21588, 0, 3, 17098, 17198, 9642, 19248,
                                                 4491, 4554, 11388, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 21798, 0, 3, 17198, 17298, 9732, 19398,
                                                 4554, 4617, 11514, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 22008, 0, 3, 17298, 17398, 9822, 19548,
                                                 4617, 4680, 11640, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22218, 0, 3, 17598, 17748, 10128, 19908,
                                                 4806, 4890, 12102, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22498, 0, 3, 17748, 17898, 10254, 20118,
                                                 4890, 4974, 12270, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 22778, 0, 3, 17898, 18048, 10380, 20328,
                                                 4974, 5058, 12438, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23058, 0, 3, 18048, 18198, 10506, 20538,
                                                 5058, 5142, 12606, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23338, 0, 3, 18198, 18348, 10632, 20748,
                                                 5142, 5226, 12774, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23618, 0, 3, 18648, 18798, 11010, 21168,
                                                 5394, 5478, 13278, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 23898, 0, 3, 18798, 18948, 11136, 21378,
                                                 5478, 5562, 13446, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24178, 0, 3, 18948, 19098, 11262, 21588,
                                                 5562, 5646, 13614, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24458, 0, 3, 19098, 19248, 11388, 21798,
                                                 5646, 5730, 13782, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 24738, 0, 3, 19248, 19398, 11514, 22008,
                                                 5730, 5814, 13950, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25018, 3, 5982, 5988, 14128, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25033, 3, 5988, 5994, 14138, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25048, 3, 5994, 6000, 14148, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25063, 3, 6000, 6006, 14158, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25078, 3, 6006, 6012, 14168, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25093, 3, 6012, 6018, 14178, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25108, 3, 6018, 6024, 14188, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25123, 3, 6024, 6030, 14198, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25138, 3, 6030, 6036, 14208, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25153, 3, 6048, 6054, 14228, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25168, 3, 6054, 6060, 14238, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25183, 3, 6060, 6066, 14248, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25198, 3, 6066, 6072, 14258, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25213, 3, 6072, 6078, 14268, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25228, 3, 6078, 6084, 14278, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25243, 3, 6084, 6090, 14288, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25258, 3, 6090, 6096, 14298, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 25273, 3, 6096, 6102, 14308, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 25288, 0, 3, 14118, 25018, 14348, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25333, 0, 3, 14128, 25033, 14378, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25378, 0, 3, 14138, 25048, 14408, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25423, 0, 3, 14148, 25063, 14438, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25468, 0, 3, 14158, 25078, 14468, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25513, 0, 3, 14168, 25093, 14498, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25558, 0, 3, 14178, 25108, 14528, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25603, 0, 3, 14188, 25123, 14558, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25648, 0, 3, 14198, 25138, 14588, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25693, 0, 3, 14218, 25153, 14648, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25738, 0, 3, 14228, 25168, 14678, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25783, 0, 3, 14238, 25183, 14708, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25828, 0, 3, 14248, 25198, 14738, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25873, 0, 3, 14258, 25213, 14768, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25918, 0, 3, 14268, 25228, 14798, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 25963, 0, 3, 14278, 25243, 14828, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26008, 0, 3, 14288, 25258, 14858, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 26053, 0, 3, 14298, 25273, 14888, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 26098, 0, 3, 14318, 25288, 6510, 6546,
                                                 14918, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26188, 0, 3, 14348, 25333, 6546, 6582,
                                                 14978, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26278, 0, 3, 14378, 25378, 6582, 6618,
                                                 15038, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26368, 0, 3, 14408, 25423, 6618, 6654,
                                                 15098, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26458, 0, 3, 14438, 25468, 6654, 6690,
                                                 15158, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26548, 0, 3, 14468, 25513, 6690, 6726,
                                                 15218, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26638, 0, 3, 14498, 25558, 6726, 6762,
                                                 15278, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26728, 0, 3, 14528, 25603, 6762, 6798,
                                                 15338, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26818, 0, 3, 14558, 25648, 6798, 6834,
                                                 15398, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26908, 0, 3, 14618, 25693, 6906, 6942,
                                                 15458, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 26998, 0, 3, 14648, 25738, 6942, 6978,
                                                 15518, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27088, 0, 3, 14678, 25783, 6978, 7014,
                                                 15578, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27178, 0, 3, 14708, 25828, 7014, 7050,
                                                 15638, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27268, 0, 3, 14738, 25873, 7050, 7086,
                                                 15698, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27358, 0, 3, 14768, 25918, 7086, 7122,
                                                 15758, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27448, 0, 3, 14798, 25963, 7122, 7158,
                                                 15818, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27538, 0, 3, 14828, 26008, 7158, 7194,
                                                 15878, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 27628, 0, 3, 14858, 26053, 7194, 7230,
                                                 15938, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 27718, 0, 3, 14978, 26278, 7302, 7362,
                                                 16098, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 27868, 0, 3, 15038, 26368, 7362, 7422,
                                                 16198, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28018, 0, 3, 15098, 26458, 7422, 7482,
                                                 16298, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28168, 0, 3, 15158, 26548, 7482, 7542,
                                                 16398, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28318, 0, 3, 15218, 26638, 7542, 7602,
                                                 16498, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28468, 0, 3, 15278, 26728, 7602, 7662,
                                                 16598, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28618, 0, 3, 15338, 26818, 7662, 7722,
                                                 16698, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28768, 0, 3, 15518, 27088, 7842, 7902,
                                                 16898, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 28918, 0, 3, 15578, 27178, 7902, 7962,
                                                 16998, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29068, 0, 3, 15638, 27268, 7962, 8022,
                                                 17098, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29218, 0, 3, 15698, 27358, 8022, 8082,
                                                 17198, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29368, 0, 3, 15758, 27448, 8082, 8142,
                                                 17298, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29518, 0, 3, 15818, 27538, 8142, 8202,
                                                 17398, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 29668, 0, 3, 15878, 27628, 8202, 8262,
                                                 17498, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 29818, 0, 3, 26098, 26188, 15998, 27718,
                                                 8382, 8472, 17598, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30043, 0, 3, 26188, 26278, 16098, 27868,
                                                 8472, 8562, 17748, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30268, 0, 3, 26278, 26368, 16198, 28018,
                                                 8562, 8652, 17898, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30493, 0, 3, 26368, 26458, 16298, 28168,
                                                 8652, 8742, 18048, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30718, 0, 3, 26458, 26548, 16398, 28318,
                                                 8742, 8832, 18198, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 30943, 0, 3, 26548, 26638, 16498, 28468,
                                                 8832, 8922, 18348, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31168, 0, 3, 26638, 26728, 16598, 28618,
                                                 8922, 9012, 18498, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31393, 0, 3, 26908, 26998, 16798, 28768,
                                                 9192, 9282, 18648, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31618, 0, 3, 26998, 27088, 16898, 28918,
                                                 9282, 9372, 18798, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 31843, 0, 3, 27088, 27178, 16998, 29068,
                                                 9372, 9462, 18948, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32068, 0, 3, 27178, 27268, 17098, 29218,
                                                 9462, 9552, 19098, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32293, 0, 3, 27268, 27358, 17198, 29368,
                                                 9552, 9642, 19248, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32518, 0, 3, 27358, 27448, 17298, 29518,
                                                 9642, 9732, 19398, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 32743, 0, 3, 27448, 27538, 17398, 29668,
                                                 9732, 9822, 19548, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 32968, 0, 3, 27718, 27868, 17748, 30268,
                                                 10002, 10128, 19908, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33283, 0, 3, 27868, 28018, 17898, 30493,
                                                 10128, 10254, 20118, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33598, 0, 3, 28018, 28168, 18048, 30718,
                                                 10254, 10380, 20328, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 33913, 0, 3, 28168, 28318, 18198, 30943,
                                                 10380, 10506, 20538, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 34228, 0, 3, 28318, 28468, 18348, 31168,
                                                 10506, 10632, 20748, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 34543, 0, 3, 28768, 28918, 18798, 31843,
                                                 10884, 11010, 21168, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 34858, 0, 3, 28918, 29068, 18948, 32068,
                                                 11010, 11136, 21378, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 35173, 0, 3, 29068, 29218, 19098, 32293,
                                                 11136, 11262, 21588, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 35488, 0, 3, 29218, 29368, 19248, 32518,
                                                 11262, 11388, 21798, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 35803, 0, 3, 29368, 29518, 19398, 32743,
                                                 11388, 11514, 22008, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 36118, 0, 3, 29818, 30043, 19698, 32968,
                                                 11766, 11934, 22218, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 36538, 0, 3, 30043, 30268, 19908, 33283,
                                                 11934, 12102, 22498, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 36958, 0, 3, 30268, 30493, 20118, 33598,
                                                 12102, 12270, 22778, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 37378, 0, 3, 30493, 30718, 20328, 33913,
                                                 12270, 12438, 23058, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 37798, 0, 3, 30718, 30943, 20538, 34228,
                                                 12438, 12606, 23338, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 38218, 0, 3, 31393, 31618, 20958, 34543,
                                                 12942, 13110, 23618, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 38638, 0, 3, 31618, 31843, 21168, 34858,
                                                 13110, 13278, 23898, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 39058, 0, 3, 31843, 32068, 21378, 35173,
                                                 13278, 13446, 24178, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 39478, 0, 3, 32068, 32293, 21588, 35488,
                                                 13446, 13614, 24458, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 39898, 0, 3, 32293, 32518, 21798, 35803,
                                                 13614, 13782, 24738, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40318, 3, 14118, 14128, 25033, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40339, 3, 14128, 14138, 25048, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40360, 3, 14138, 14148, 25063, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40381, 3, 14148, 14158, 25078, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40402, 3, 14158, 14168, 25093, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40423, 3, 14168, 14178, 25108, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40444, 3, 14178, 14188, 25123, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40465, 3, 14188, 14198, 25138, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40486, 3, 14218, 14228, 25168, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40507, 3, 14228, 14238, 25183, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40528, 3, 14238, 14248, 25198, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40549, 3, 14248, 14258, 25213, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40570, 3, 14258, 14268, 25228, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40591, 3, 14268, 14278, 25243, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40612, 3, 14278, 14288, 25258, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 40633, 3, 14288, 14298, 25273, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 40654, 0, 3, 25018, 40318, 25333, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 40717, 0, 3, 25033, 40339, 25378, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 40780, 0, 3, 25048, 40360, 25423, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 40843, 0, 3, 25063, 40381, 25468, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 40906, 0, 3, 25078, 40402, 25513, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 40969, 0, 3, 25093, 40423, 25558, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41032, 0, 3, 25108, 40444, 25603, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41095, 0, 3, 25123, 40465, 25648, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41158, 0, 3, 25153, 40486, 25738, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41221, 0, 3, 25168, 40507, 25783, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41284, 0, 3, 25183, 40528, 25828, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41347, 0, 3, 25198, 40549, 25873, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41410, 0, 3, 25213, 40570, 25918, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41473, 0, 3, 25228, 40591, 25963, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41536, 0, 3, 25243, 40612, 26008, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 41599, 0, 3, 25258, 40633, 26053, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 41662, 0, 3, 25333, 40717, 14918, 14978,
                                                 26278, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 41788, 0, 3, 25378, 40780, 14978, 15038,
                                                 26368, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 41914, 0, 3, 25423, 40843, 15038, 15098,
                                                 26458, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 42040, 0, 3, 25468, 40906, 15098, 15158,
                                                 26548, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 42166, 0, 3, 25513, 40969, 15158, 15218,
                                                 26638, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 42292, 0, 3, 25558, 41032, 15218, 15278,
                                                 26728, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 42418, 0, 3, 25603, 41095, 15278, 15338,
                                                 26818, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 42544, 0, 3, 25738, 41221, 15458, 15518,
                                                 27088, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 42670, 0, 3, 25783, 41284, 15518, 15578,
                                                 27178, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 42796, 0, 3, 25828, 41347, 15578, 15638,
                                                 27268, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 42922, 0, 3, 25873, 41410, 15638, 15698,
                                                 27358, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43048, 0, 3, 25918, 41473, 15698, 15758,
                                                 27448, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43174, 0, 3, 25963, 41536, 15758, 15818,
                                                 27538, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 43300, 0, 3, 26008, 41599, 15818, 15878,
                                                 27628, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 43426, 0, 3, 26278, 41788, 15998, 16098,
                                                 27868, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 43636, 0, 3, 26368, 41914, 16098, 16198,
                                                 28018, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 43846, 0, 3, 26458, 42040, 16198, 16298,
                                                 28168, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44056, 0, 3, 26548, 42166, 16298, 16398,
                                                 28318, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44266, 0, 3, 26638, 42292, 16398, 16498,
                                                 28468, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44476, 0, 3, 26728, 42418, 16498, 16598,
                                                 28618, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44686, 0, 3, 27088, 42670, 16798, 16898,
                                                 28918, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 44896, 0, 3, 27178, 42796, 16898, 16998,
                                                 29068, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45106, 0, 3, 27268, 42922, 16998, 17098,
                                                 29218, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45316, 0, 3, 27358, 43048, 17098, 17198,
                                                 29368, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45526, 0, 3, 27448, 43174, 17198, 17298,
                                                 29518, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 45736, 0, 3, 27538, 43300, 17298, 17398,
                                                 29668, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 45946, 0, 3, 41662, 41788, 27868, 43636,
                                                 17598, 17748, 30268, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46261, 0, 3, 41788, 41914, 28018, 43846,
                                                 17748, 17898, 30493, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46576, 0, 3, 41914, 42040, 28168, 44056,
                                                 17898, 18048, 30718, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 46891, 0, 3, 42040, 42166, 28318, 44266,
                                                 18048, 18198, 30943, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47206, 0, 3, 42166, 42292, 28468, 44476,
                                                 18198, 18348, 31168, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47521, 0, 3, 42544, 42670, 28918, 44896,
                                                 18648, 18798, 31843, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 47836, 0, 3, 42670, 42796, 29068, 45106,
                                                 18798, 18948, 32068, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 48151, 0, 3, 42796, 42922, 29218, 45316,
                                                 18948, 19098, 32293, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 48466, 0, 3, 42922, 43048, 29368, 45526,
                                                 19098, 19248, 32518, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 48781, 0, 3, 43048, 43174, 29518, 45736,
                                                 19248, 19398, 32743, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 49096, 0, 3, 43426, 43636, 30268, 46261,
                                                 19698, 19908, 33283, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 49537, 0, 3, 43636, 43846, 30493, 46576,
                                                 19908, 20118, 33598, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 49978, 0, 3, 43846, 44056, 30718, 46891,
                                                 20118, 20328, 33913, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 50419, 0, 3, 44056, 44266, 30943, 47206,
                                                 20328, 20538, 34228, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 50860, 0, 3, 44686, 44896, 31843, 47836,
                                                 20958, 21168, 34858, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 51301, 0, 3, 44896, 45106, 32068, 48151,
                                                 21168, 21378, 35173, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 51742, 0, 3, 45106, 45316, 32293, 48466,
                                                 21378, 21588, 35488, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 52183, 0, 3, 45316, 45526, 32518, 48781,
                                                 21588, 21798, 35803, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 52624, 0, 3, 45946, 46261, 33283, 49537,
                                                 22218, 22498, 36958, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 53212, 0, 3, 46261, 46576, 33598, 49978,
                                                 22498, 22778, 37378, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 53800, 0, 3, 46576, 46891, 33913, 50419,
                                                 22778, 23058, 37798, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 54388, 0, 3, 47521, 47836, 34858, 51301,
                                                 23618, 23898, 39058, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 54976, 0, 3, 47836, 48151, 35173, 51742,
                                                 23898, 24178, 39478, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 55564, 0, 3, 48151, 48466, 35488, 52183,
                                                 24178, 24458, 39898, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56152, 3, 25018, 25033, 40339, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56180, 3, 25033, 25048, 40360, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56208, 3, 25048, 25063, 40381, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56236, 3, 25063, 25078, 40402, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56264, 3, 25078, 25093, 40423, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56292, 3, 25093, 25108, 40444, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56320, 3, 25108, 25123, 40465, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56348, 3, 25153, 25168, 40507, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56376, 3, 25168, 25183, 40528, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56404, 3, 25183, 25198, 40549, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56432, 3, 25198, 25213, 40570, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56460, 3, 25213, 25228, 40591, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56488, 3, 25228, 25243, 40612, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 56516, 3, 25243, 25258, 40633, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 56544, 0, 3, 40318, 56152, 40717, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56628, 0, 3, 40339, 56180, 40780, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56712, 0, 3, 40360, 56208, 40843, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56796, 0, 3, 40381, 56236, 40906, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56880, 0, 3, 40402, 56264, 40969, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 56964, 0, 3, 40423, 56292, 41032, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57048, 0, 3, 40444, 56320, 41095, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57132, 0, 3, 40486, 56348, 41221, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57216, 0, 3, 40507, 56376, 41284, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57300, 0, 3, 40528, 56404, 41347, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57384, 0, 3, 40549, 56432, 41410, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57468, 0, 3, 40570, 56460, 41473, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57552, 0, 3, 40591, 56488, 41536, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 57636, 0, 3, 40612, 56516, 41599, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 57720, 0, 3, 40654, 56544, 26098, 26188,
                                                 41662, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 57888, 0, 3, 40717, 56628, 26188, 26278,
                                                 41788, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58056, 0, 3, 40780, 56712, 26278, 26368,
                                                 41914, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58224, 0, 3, 40843, 56796, 26368, 26458,
                                                 42040, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58392, 0, 3, 40906, 56880, 26458, 26548,
                                                 42166, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58560, 0, 3, 40969, 56964, 26548, 26638,
                                                 42292, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58728, 0, 3, 41032, 57048, 26638, 26728,
                                                 42418, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 58896, 0, 3, 41158, 57132, 26908, 26998,
                                                 42544, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59064, 0, 3, 41221, 57216, 26998, 27088,
                                                 42670, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59232, 0, 3, 41284, 57300, 27088, 27178,
                                                 42796, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59400, 0, 3, 41347, 57384, 27178, 27268,
                                                 42922, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59568, 0, 3, 41410, 57468, 27268, 27358,
                                                 43048, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59736, 0, 3, 41473, 57552, 27358, 27448,
                                                 43174, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 59904, 0, 3, 41536, 57636, 27448, 27538,
                                                 43300, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60072, 0, 3, 41788, 58056, 27718, 27868,
                                                 43636, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60352, 0, 3, 41914, 58224, 27868, 28018,
                                                 43846, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60632, 0, 3, 42040, 58392, 28018, 28168,
                                                 44056, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 60912, 0, 3, 42166, 58560, 28168, 28318,
                                                 44266, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 61192, 0, 3, 42292, 58728, 28318, 28468,
                                                 44476, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 61472, 0, 3, 42670, 59232, 28768, 28918,
                                                 44896, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 61752, 0, 3, 42796, 59400, 28918, 29068,
                                                 45106, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 62032, 0, 3, 42922, 59568, 29068, 29218,
                                                 45316, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 62312, 0, 3, 43048, 59736, 29218, 29368,
                                                 45526, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 62592, 0, 3, 43174, 59904, 29368, 29518,
                                                 45736, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 62872, 0, 3, 57720, 57888, 43426, 60072,
                                                 29818, 30043, 45946, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 63292, 0, 3, 57888, 58056, 43636, 60352,
                                                 30043, 30268, 46261, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 63712, 0, 3, 58056, 58224, 43846, 60632,
                                                 30268, 30493, 46576, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 64132, 0, 3, 58224, 58392, 44056, 60912,
                                                 30493, 30718, 46891, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 64552, 0, 3, 58392, 58560, 44266, 61192,
                                                 30718, 30943, 47206, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 64972, 0, 3, 58896, 59064, 44686, 61472,
                                                 31393, 31618, 47521, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 65392, 0, 3, 59064, 59232, 44896, 61752,
                                                 31618, 31843, 47836, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 65812, 0, 3, 59232, 59400, 45106, 62032,
                                                 31843, 32068, 48151, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 66232, 0, 3, 59400, 59568, 45316, 62312,
                                                 32068, 32293, 48466, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 66652, 0, 3, 59568, 59736, 45526, 62592,
                                                 32293, 32518, 48781, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 67072, 0, 3, 60072, 60352, 46261, 63712,
                                                 32968, 33283, 49537, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 67660, 0, 3, 60352, 60632, 46576, 64132,
                                                 33283, 33598, 49978, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 68248, 0, 3, 60632, 60912, 46891, 64552,
                                                 33598, 33913, 50419, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 68836, 0, 3, 61472, 61752, 47836, 65812,
                                                 34543, 34858, 51301, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 69424, 0, 3, 61752, 62032, 48151, 66232,
                                                 34858, 35173, 51742, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 70012, 0, 3, 62032, 62312, 48466, 66652,
                                                 35173, 35488, 52183, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 70600, 0, 3, 62872, 63292, 49096, 67072,
                                                 36118, 36538, 52624, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 71384, 0, 3, 63292, 63712, 49537, 67660,
                                                 36538, 36958, 53212, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 72168, 0, 3, 63712, 64132, 49978, 68248,
                                                 36958, 37378, 53800, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 72952, 0, 3, 64972, 65392, 50860, 68836,
                                                 38218, 38638, 54388, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 73736, 0, 3, 65392, 65812, 51301, 69424,
                                                 38638, 39058, 54976, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 74520, 0, 3, 65812, 66232, 51742, 70012,
                                                 39058, 39478, 55564, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75304, 3, 40318, 40339, 56180, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75340, 3, 40339, 40360, 56208, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75376, 3, 40360, 40381, 56236, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75412, 3, 40381, 40402, 56264, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75448, 3, 40402, 40423, 56292, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75484, 3, 40423, 40444, 56320, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75520, 3, 40486, 40507, 56376, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75556, 3, 40507, 40528, 56404, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75592, 3, 40528, 40549, 56432, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75628, 3, 40549, 40570, 56460, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75664, 3, 40570, 40591, 56488, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 75700, 3, 40591, 40612, 56516, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 75736, 0, 3, 56152, 75304, 56628, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 75844, 0, 3, 56180, 75340, 56712, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 75952, 0, 3, 56208, 75376, 56796, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76060, 0, 3, 56236, 75412, 56880, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76168, 0, 3, 56264, 75448, 56964, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76276, 0, 3, 56292, 75484, 57048, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76384, 0, 3, 56348, 75520, 57216, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76492, 0, 3, 56376, 75556, 57300, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76600, 0, 3, 56404, 75592, 57384, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76708, 0, 3, 56432, 75628, 57468, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76816, 0, 3, 56460, 75664, 57552, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 76924, 0, 3, 56488, 75700, 57636, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 77032, 0, 3, 56628, 75844, 41662, 41788,
                                                 58056, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 77248, 0, 3, 56712, 75952, 41788, 41914,
                                                 58224, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 77464, 0, 3, 56796, 76060, 41914, 42040,
                                                 58392, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 77680, 0, 3, 56880, 76168, 42040, 42166,
                                                 58560, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 77896, 0, 3, 56964, 76276, 42166, 42292,
                                                 58728, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 78112, 0, 3, 57216, 76492, 42544, 42670,
                                                 59232, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 78328, 0, 3, 57300, 76600, 42670, 42796,
                                                 59400, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 78544, 0, 3, 57384, 76708, 42796, 42922,
                                                 59568, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 78760, 0, 3, 57468, 76816, 42922, 43048,
                                                 59736, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 78976, 0, 3, 57552, 76924, 43048, 43174,
                                                 59904, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 79192, 0, 3, 58056, 77248, 43426, 43636,
                                                 60352, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 79552, 0, 3, 58224, 77464, 43636, 43846,
                                                 60632, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 79912, 0, 3, 58392, 77680, 43846, 44056,
                                                 60912, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 80272, 0, 3, 58560, 77896, 44056, 44266,
                                                 61192, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 80632, 0, 3, 59232, 78328, 44686, 44896,
                                                 61752, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 80992, 0, 3, 59400, 78544, 44896, 45106,
                                                 62032, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 81352, 0, 3, 59568, 78760, 45106, 45316,
                                                 62312, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 81712, 0, 3, 59736, 78976, 45316, 45526,
                                                 62592, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 82072, 0, 3, 77032, 77248, 60352, 79552,
                                                 45946, 46261, 63712, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 82612, 0, 3, 77248, 77464, 60632, 79912,
                                                 46261, 46576, 64132, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 83152, 0, 3, 77464, 77680, 60912, 80272,
                                                 46576, 46891, 64552, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 83692, 0, 3, 78112, 78328, 61752, 80992,
                                                 47521, 47836, 65812, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 84232, 0, 3, 78328, 78544, 62032, 81352,
                                                 47836, 48151, 66232, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 84772, 0, 3, 78544, 78760, 62312, 81712,
                                                 48151, 48466, 66652, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 85312, 0, 3, 79192, 79552, 63712, 82612,
                                                 49096, 49537, 67660, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 86068, 0, 3, 79552, 79912, 64132, 83152,
                                                 49537, 49978, 68248, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 86824, 0, 3, 80632, 80992, 65812, 84232,
                                                 50860, 51301, 69424, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 87580, 0, 3, 80992, 81352, 66232, 84772,
                                                 51301, 51742, 70012, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 88336, 0, 3, 82072, 82612, 67660, 86068,
                                                 52624, 53212, 72168, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 89344, 0, 3, 83692, 84232, 69424, 87580,
                                                 54388, 54976, 74520, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90352, 3, 56152, 56180, 75340, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90397, 3, 56180, 56208, 75376, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90442, 3, 56208, 56236, 75412, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90487, 3, 56236, 56264, 75448, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90532, 3, 56264, 56292, 75484, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90577, 3, 56348, 56376, 75556, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90622, 3, 56376, 56404, 75592, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90667, 3, 56404, 56432, 75628, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90712, 3, 56432, 56460, 75664, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 90757, 3, 56460, 56488, 75700, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 90802, 0, 3, 75304, 90352, 75844, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 90937, 0, 3, 75340, 90397, 75952, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 91072, 0, 3, 75376, 90442, 76060, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 91207, 0, 3, 75412, 90487, 76168, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 91342, 0, 3, 75448, 90532, 76276, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 91477, 0, 3, 75520, 90577, 76492, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 91612, 0, 3, 75556, 90622, 76600, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 91747, 0, 3, 75592, 90667, 76708, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 91882, 0, 3, 75628, 90712, 76816, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 92017, 0, 3, 75664, 90757, 76924, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 92152, 0, 3, 75736, 90802, 57720, 57888,
                                                 77032, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 92422, 0, 3, 75844, 90937, 57888, 58056,
                                                 77248, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 92692, 0, 3, 75952, 91072, 58056, 58224,
                                                 77464, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 92962, 0, 3, 76060, 91207, 58224, 58392,
                                                 77680, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 93232, 0, 3, 76168, 91342, 58392, 58560,
                                                 77896, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 93502, 0, 3, 76384, 91477, 58896, 59064,
                                                 78112, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 93772, 0, 3, 76492, 91612, 59064, 59232,
                                                 78328, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 94042, 0, 3, 76600, 91747, 59232, 59400,
                                                 78544, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 94312, 0, 3, 76708, 91882, 59400, 59568,
                                                 78760, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 94582, 0, 3, 76816, 92017, 59568, 59736,
                                                 78976, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 94852, 0, 3, 77248, 92692, 60072, 60352,
                                                 79552, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 95302, 0, 3, 77464, 92962, 60352, 60632,
                                                 79912, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 95752, 0, 3, 77680, 93232, 60632, 60912,
                                                 80272, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 96202, 0, 3, 78328, 94042, 61472, 61752,
                                                 80992, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 96652, 0, 3, 78544, 94312, 61752, 62032,
                                                 81352, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 97102, 0, 3, 78760, 94582, 62032, 62312,
                                                 81712, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 97552, 0, 3, 92152, 92422, 79192, 94852,
                                                 62872, 63292, 82072, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 98227, 0, 3, 92422, 92692, 79552, 95302,
                                                 63292, 63712, 82612, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 98902, 0, 3, 92692, 92962, 79912, 95752,
                                                 63712, 64132, 83152, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 99577, 0, 3, 93502, 93772, 80632, 96202,
                                                 64972, 65392, 83692, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 100252, 0, 3, 93772, 94042, 80992,
                                                 96652, 65392, 65812, 84232, ncols, alpha, beta,
                                                 p);

            compute_prim_gl_electron_repulsion_0(buffer, 100927, 0, 3, 94042, 94312, 81352,
                                                 97102, 65812, 66232, 84772, ncols, alpha, beta,
                                                 p);

            compute_prim_hl_electron_repulsion_0(buffer, 101602, 0, 3, 94852, 95302, 82612,
                                                 98902, 67072, 67660, 86068, ncols, alpha, beta,
                                                 p);

            compute_prim_hl_electron_repulsion_0(buffer, 102547, 0, 3, 96202, 96652, 84232,
                                                 100927, 68836, 69424, 87580, ncols, alpha, beta,
                                                 p);

            compute_prim_il_electron_repulsion_0(buffer, 103492, 0, 3, 97552, 98227, 85312,
                                                 101602, 70600, 71384, 88336, ncols, alpha, beta,
                                                 p);

            compute_prim_il_electron_repulsion_0(buffer, 104752, 0, 3, 99577, 100252, 86824,
                                                 102547, 72952, 73736, 89344, ncols, alpha, beta,
                                                 p);

            simdgeo::geom_h_x(buffer, 106012, 99577, 104752, 1, 45, ncols, alpha);

            simdgeo::geom_h_y(buffer, 106957, 99577, 104752, 1, 45, ncols, alpha);

            simdgeo::geom_h_z(buffer, 107902, 99577, 104752, 1, 45, ncols, alpha);

            simdgeo::geom_h_x(buffer, 108847, 97552, 103492, 1, 45, ncols, alpha);

            simdgeo::geom_h_y(buffer, 109792, 97552, 103492, 1, 45, ncols, alpha);

            simdgeo::geom_h_z(buffer, 110737, 97552, 103492, 1, 45, ncols, alpha);

            simdfunc::contract_primitives(buffer, 111682, 108847, 2835, ncols);

            simdfunc::contract_primitives(buffer, 114517, 106012, 2835, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 117352, 114517, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 117352, 17, nmax);

    simdtrf::transform_l_inner(buffer, 117352, 115462, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 187 * nvalues, nvalues, buffer, 117352, 17, nmax);

    simdtrf::transform_l_inner(buffer, 117352, 116407, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 374 * nvalues, nvalues, buffer, 117352, 17, nmax);

    simdtrf::transform_l_inner(buffer, 117352, 111682, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 561 * nvalues, nvalues, buffer, 117352, 17, nmax);

    simdtrf::transform_l_inner(buffer, 117352, 112627, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 748 * nvalues, nvalues, buffer, 117352, 17, nmax);

    simdtrf::transform_l_inner(buffer, 117352, 113572, 21, 1, nmax);

    simdtrf::transform_h_outer(values + 935 * nvalues, nvalues, buffer, 117352, 17, nmax);
}

}  // namespace simdt2ceri
