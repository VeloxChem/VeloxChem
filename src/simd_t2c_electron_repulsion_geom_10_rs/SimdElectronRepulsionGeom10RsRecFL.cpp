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


#include "SimdElectronRepulsionGeom10RsRecFL.hpp"

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
#include "SimdGeometryF1.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_fl_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_10_fl_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 43892, 41022, 2700, nvalues);

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

            simdfunc::compute_full_erf_boys_function(buffer, coordinates, 6, 12, ncols, fj, mu,
                                                     omega);

            simdfunc::compute_full_boys_function(buffer, coordinates, 20, 12, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 64, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 67, 0, 23, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 70, 0, 24, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 73, 0, 25, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 76, 0, 26, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 79, 0, 27, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 82, 0, 28, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 85, 0, 29, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 88, 0, 30, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 91, 0, 31, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 94, 0, 32, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 97, 0, 33, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 7, 8, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 8, 9, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 9, 10, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 10, 11, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 11, 12, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 12, 13, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 13, 14, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 142, 0, 14, 15, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 148, 0, 15, 16, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 154, 0, 16, 17, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 160, 0, 17, 18, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 166, 0, 21, 22, 67, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 172, 0, 22, 23, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 178, 0, 23, 24, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 184, 0, 24, 25, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 190, 0, 25, 26, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 196, 0, 26, 27, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 202, 0, 27, 28, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 208, 0, 28, 29, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 214, 0, 29, 30, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 220, 0, 30, 31, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 226, 0, 31, 32, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 34, 37, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 242, 0, 37, 40, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 252, 0, 40, 43, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 262, 0, 43, 46, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 272, 0, 46, 49, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 282, 0, 49, 52, 142, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 292, 0, 52, 55, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 302, 0, 55, 58, 154, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 312, 0, 58, 61, 160, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 322, 0, 67, 70, 178, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 332, 0, 70, 73, 184, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 342, 0, 73, 76, 190, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 352, 0, 76, 79, 196, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 362, 0, 79, 82, 202, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 372, 0, 82, 85, 208, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 382, 0, 85, 88, 214, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 392, 0, 88, 91, 220, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 402, 0, 91, 94, 226, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 412, 0, 100, 106, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 427, 0, 106, 112, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 442, 0, 112, 118, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 457, 0, 118, 124, 262, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 472, 0, 124, 130, 272, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 487, 0, 130, 136, 282, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 502, 0, 136, 142, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 517, 0, 142, 148, 302, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 532, 0, 148, 154, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 547, 0, 166, 172, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 562, 0, 172, 178, 332, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 577, 0, 178, 184, 342, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 592, 0, 184, 190, 352, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 607, 0, 190, 196, 362, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 622, 0, 196, 202, 372, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 637, 0, 202, 208, 382, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 652, 0, 208, 214, 392, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 667, 0, 214, 220, 402, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 682, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 685, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 688, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 691, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 694, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 697, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 700, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 703, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 706, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 709, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 712, 3, 24, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 715, 3, 25, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 718, 3, 26, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 721, 3, 27, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 724, 3, 28, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 727, 3, 29, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 730, 3, 30, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 733, 3, 31, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 736, 3, 32, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 739, 3, 33, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 742, 3, 9, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 751, 3, 10, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 760, 3, 11, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 769, 3, 12, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 778, 3, 13, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 787, 3, 14, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 796, 3, 15, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 805, 3, 16, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 814, 3, 17, 61, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 823, 3, 18, 64, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 832, 3, 23, 70, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 841, 3, 24, 73, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 850, 3, 25, 76, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 859, 3, 26, 79, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 868, 3, 27, 82, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 877, 3, 28, 85, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 886, 3, 29, 88, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 895, 3, 30, 91, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 904, 3, 31, 94, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 913, 3, 32, 97, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 922, 0, 3, 37, 751, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 940, 0, 3, 40, 760, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 958, 0, 3, 43, 769, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 976, 0, 3, 46, 778, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 994, 0, 3, 49, 787, 136, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1012, 0, 3, 52, 796, 142, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1030, 0, 3, 55, 805, 148, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1048, 0, 3, 58, 814, 154, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1066, 0, 3, 61, 823, 160, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1084, 0, 3, 70, 841, 178, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1102, 0, 3, 73, 850, 184, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1120, 0, 3, 76, 859, 190, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1138, 0, 3, 79, 868, 196, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1156, 0, 3, 82, 877, 202, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1174, 0, 3, 85, 886, 208, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1192, 0, 3, 88, 895, 214, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1210, 0, 3, 91, 904, 220, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1228, 0, 3, 94, 913, 226, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1246, 0, 3, 112, 940, 242, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1276, 0, 3, 118, 958, 252, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1306, 0, 3, 124, 976, 262, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1336, 0, 3, 130, 994, 272, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1366, 0, 3, 136, 1012, 282, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1396, 0, 3, 142, 1030, 292, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1426, 0, 3, 148, 1048, 302, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1456, 0, 3, 154, 1066, 312, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1486, 0, 3, 178, 1102, 332, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1516, 0, 3, 184, 1120, 342, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1546, 0, 3, 190, 1138, 352, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1576, 0, 3, 196, 1156, 362, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1606, 0, 3, 202, 1174, 372, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1636, 0, 3, 208, 1192, 382, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1666, 0, 3, 214, 1210, 392, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1696, 0, 3, 220, 1228, 402, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1726, 0, 3, 242, 1276, 442, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1771, 0, 3, 252, 1306, 457, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1816, 0, 3, 262, 1336, 472, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1861, 0, 3, 272, 1366, 487, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1906, 0, 3, 282, 1396, 502, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1951, 0, 3, 292, 1426, 517, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 1996, 0, 3, 302, 1456, 532, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2041, 0, 3, 332, 1516, 577, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2086, 0, 3, 342, 1546, 592, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2131, 0, 3, 352, 1576, 607, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2176, 0, 3, 362, 1606, 622, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2221, 0, 3, 372, 1636, 637, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2266, 0, 3, 382, 1666, 652, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2311, 0, 3, 392, 1696, 667, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2356, 3, 9, 10, 685, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 2362, 3, 10, 11, 688, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2368, 3, 11, 12, 691, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2374, 3, 12, 13, 694, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2380, 3, 13, 14, 697, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2386, 3, 14, 15, 700, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2392, 3, 15, 16, 703, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2398, 3, 16, 17, 706, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2404, 3, 17, 18, 709, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2410, 3, 23, 24, 715, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2416, 3, 24, 25, 718, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2422, 3, 25, 26, 721, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2428, 3, 26, 27, 724, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2434, 3, 27, 28, 727, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2440, 3, 28, 29, 730, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2446, 3, 29, 30, 733, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2452, 3, 30, 31, 736, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 2458, 3, 31, 32, 739, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2464, 0, 3, 682, 2356, 751, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2482, 0, 3, 685, 2362, 760, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2500, 0, 3, 688, 2368, 769, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2518, 0, 3, 691, 2374, 778, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2536, 0, 3, 694, 2380, 787, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2554, 0, 3, 697, 2386, 796, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2572, 0, 3, 700, 2392, 805, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2590, 0, 3, 703, 2398, 814, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2608, 0, 3, 706, 2404, 823, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2626, 0, 3, 712, 2410, 841, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2644, 0, 3, 715, 2416, 850, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2662, 0, 3, 718, 2422, 859, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2680, 0, 3, 721, 2428, 868, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2698, 0, 3, 724, 2434, 877, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2716, 0, 3, 727, 2440, 886, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2734, 0, 3, 730, 2446, 895, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2752, 0, 3, 733, 2452, 904, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 2770, 0, 3, 736, 2458, 913, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2788, 0, 3, 742, 2464, 100, 106, 922,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2824, 0, 3, 751, 2482, 106, 112, 940,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2860, 0, 3, 760, 2500, 112, 118, 958,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2896, 0, 3, 769, 2518, 118, 124, 976,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2932, 0, 3, 778, 2536, 124, 130, 994,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 2968, 0, 3, 787, 2554, 130, 136, 1012,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3004, 0, 3, 796, 2572, 136, 142, 1030,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3040, 0, 3, 805, 2590, 142, 148, 1048,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3076, 0, 3, 814, 2608, 148, 154, 1066,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3112, 0, 3, 832, 2626, 166, 172, 1084,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3148, 0, 3, 841, 2644, 172, 178, 1102,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3184, 0, 3, 850, 2662, 178, 184, 1120,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3220, 0, 3, 859, 2680, 184, 190, 1138,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3256, 0, 3, 868, 2698, 190, 196, 1156,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3292, 0, 3, 877, 2716, 196, 202, 1174,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3328, 0, 3, 886, 2734, 202, 208, 1192,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3364, 0, 3, 895, 2752, 208, 214, 1210,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 3400, 0, 3, 904, 2770, 214, 220, 1228,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3436, 0, 3, 940, 2860, 232, 242, 1276,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3496, 0, 3, 958, 2896, 242, 252, 1306,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3556, 0, 3, 976, 2932, 252, 262, 1336,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3616, 0, 3, 994, 2968, 262, 272, 1366,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3676, 0, 3, 1012, 3004, 272, 282, 1396,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3736, 0, 3, 1030, 3040, 282, 292, 1426,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3796, 0, 3, 1048, 3076, 292, 302, 1456,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3856, 0, 3, 1102, 3184, 322, 332, 1516,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3916, 0, 3, 1120, 3220, 332, 342, 1546,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 3976, 0, 3, 1138, 3256, 342, 352, 1576,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4036, 0, 3, 1156, 3292, 352, 362, 1606,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4096, 0, 3, 1174, 3328, 362, 372, 1636,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4156, 0, 3, 1192, 3364, 372, 382, 1666,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 4216, 0, 3, 1210, 3400, 382, 392, 1696,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4276, 0, 3, 2788, 2824, 1246, 3436, 412,
                                                 427, 1726, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4366, 0, 3, 2824, 2860, 1276, 3496, 427,
                                                 442, 1771, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4456, 0, 3, 2860, 2896, 1306, 3556, 442,
                                                 457, 1816, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4546, 0, 3, 2896, 2932, 1336, 3616, 457,
                                                 472, 1861, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4636, 0, 3, 2932, 2968, 1366, 3676, 472,
                                                 487, 1906, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4726, 0, 3, 2968, 3004, 1396, 3736, 487,
                                                 502, 1951, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4816, 0, 3, 3004, 3040, 1426, 3796, 502,
                                                 517, 1996, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4906, 0, 3, 3112, 3148, 1486, 3856, 547,
                                                 562, 2041, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 4996, 0, 3, 3148, 3184, 1516, 3916, 562,
                                                 577, 2086, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5086, 0, 3, 3184, 3220, 1546, 3976, 577,
                                                 592, 2131, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5176, 0, 3, 3220, 3256, 1576, 4036, 592,
                                                 607, 2176, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5266, 0, 3, 3256, 3292, 1606, 4096, 607,
                                                 622, 2221, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5356, 0, 3, 3292, 3328, 1636, 4156, 622,
                                                 637, 2266, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 5446, 0, 3, 3328, 3364, 1666, 4216, 637,
                                                 652, 2311, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5536, 3, 682, 685, 2362, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5546, 3, 685, 688, 2368, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5556, 3, 688, 691, 2374, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5566, 3, 691, 694, 2380, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5576, 3, 694, 697, 2386, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5586, 3, 697, 700, 2392, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5596, 3, 700, 703, 2398, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5606, 3, 703, 706, 2404, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5616, 3, 712, 715, 2416, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5626, 3, 715, 718, 2422, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5636, 3, 718, 721, 2428, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5646, 3, 721, 724, 2434, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5656, 3, 724, 727, 2440, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5666, 3, 727, 730, 2446, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5676, 3, 730, 733, 2452, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 5686, 3, 733, 736, 2458, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 5696, 0, 3, 2356, 5536, 2482, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5726, 0, 3, 2362, 5546, 2500, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5756, 0, 3, 2368, 5556, 2518, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5786, 0, 3, 2374, 5566, 2536, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5816, 0, 3, 2380, 5576, 2554, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5846, 0, 3, 2386, 5586, 2572, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5876, 0, 3, 2392, 5596, 2590, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5906, 0, 3, 2398, 5606, 2608, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5936, 0, 3, 2410, 5616, 2644, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5966, 0, 3, 2416, 5626, 2662, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 5996, 0, 3, 2422, 5636, 2680, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6026, 0, 3, 2428, 5646, 2698, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6056, 0, 3, 2434, 5656, 2716, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6086, 0, 3, 2440, 5666, 2734, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6116, 0, 3, 2446, 5676, 2752, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 6146, 0, 3, 2452, 5686, 2770, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 6176, 0, 3, 2482, 5726, 922, 940, 2860,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6236, 0, 3, 2500, 5756, 940, 958, 2896,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6296, 0, 3, 2518, 5786, 958, 976, 2932,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6356, 0, 3, 2536, 5816, 976, 994, 2968,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6416, 0, 3, 2554, 5846, 994, 1012, 3004,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6476, 0, 3, 2572, 5876, 1012, 1030,
                                                 3040, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6536, 0, 3, 2590, 5906, 1030, 1048,
                                                 3076, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6596, 0, 3, 2644, 5966, 1084, 1102,
                                                 3184, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6656, 0, 3, 2662, 5996, 1102, 1120,
                                                 3220, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6716, 0, 3, 2680, 6026, 1120, 1138,
                                                 3256, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6776, 0, 3, 2698, 6056, 1138, 1156,
                                                 3292, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6836, 0, 3, 2716, 6086, 1156, 1174,
                                                 3328, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6896, 0, 3, 2734, 6116, 1174, 1192,
                                                 3364, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 6956, 0, 3, 2752, 6146, 1192, 1210,
                                                 3400, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7016, 0, 3, 2860, 6236, 1246, 1276,
                                                 3496, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7116, 0, 3, 2896, 6296, 1276, 1306,
                                                 3556, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7216, 0, 3, 2932, 6356, 1306, 1336,
                                                 3616, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7316, 0, 3, 2968, 6416, 1336, 1366,
                                                 3676, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7416, 0, 3, 3004, 6476, 1366, 1396,
                                                 3736, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7516, 0, 3, 3040, 6536, 1396, 1426,
                                                 3796, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7616, 0, 3, 3184, 6656, 1486, 1516,
                                                 3916, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7716, 0, 3, 3220, 6716, 1516, 1546,
                                                 3976, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7816, 0, 3, 3256, 6776, 1546, 1576,
                                                 4036, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 7916, 0, 3, 3292, 6836, 1576, 1606,
                                                 4096, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8016, 0, 3, 3328, 6896, 1606, 1636,
                                                 4156, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 8116, 0, 3, 3364, 6956, 1636, 1666,
                                                 4216, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8216, 0, 3, 6176, 6236, 3496, 7116,
                                                 1726, 1771, 4456, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8366, 0, 3, 6236, 6296, 3556, 7216,
                                                 1771, 1816, 4546, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8516, 0, 3, 6296, 6356, 3616, 7316,
                                                 1816, 1861, 4636, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8666, 0, 3, 6356, 6416, 3676, 7416,
                                                 1861, 1906, 4726, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8816, 0, 3, 6416, 6476, 3736, 7516,
                                                 1906, 1951, 4816, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 8966, 0, 3, 6596, 6656, 3916, 7716,
                                                 2041, 2086, 5086, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9116, 0, 3, 6656, 6716, 3976, 7816,
                                                 2086, 2131, 5176, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9266, 0, 3, 6716, 6776, 4036, 7916,
                                                 2131, 2176, 5266, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9416, 0, 3, 6776, 6836, 4096, 8016,
                                                 2176, 2221, 5356, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 9566, 0, 3, 6836, 6896, 4156, 8116,
                                                 2221, 2266, 5446, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9716, 3, 2356, 2362, 5546, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9731, 3, 2362, 2368, 5556, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9746, 3, 2368, 2374, 5566, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9761, 3, 2374, 2380, 5576, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9776, 3, 2380, 2386, 5586, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9791, 3, 2386, 2392, 5596, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9806, 3, 2392, 2398, 5606, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9821, 3, 2410, 2416, 5626, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9836, 3, 2416, 2422, 5636, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9851, 3, 2422, 2428, 5646, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9866, 3, 2428, 2434, 5656, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9881, 3, 2434, 2440, 5666, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9896, 3, 2440, 2446, 5676, ncols, alpha,
                                                 beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 9911, 3, 2446, 2452, 5686, ncols, alpha,
                                                 beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 9926, 0, 3, 5536, 9716, 5726, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 9971, 0, 3, 5546, 9731, 5756, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10016, 0, 3, 5556, 9746, 5786, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10061, 0, 3, 5566, 9761, 5816, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10106, 0, 3, 5576, 9776, 5846, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10151, 0, 3, 5586, 9791, 5876, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10196, 0, 3, 5596, 9806, 5906, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10241, 0, 3, 5616, 9821, 5966, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10286, 0, 3, 5626, 9836, 5996, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10331, 0, 3, 5636, 9851, 6026, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10376, 0, 3, 5646, 9866, 6056, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10421, 0, 3, 5656, 9881, 6086, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10466, 0, 3, 5666, 9896, 6116, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 10511, 0, 3, 5676, 9911, 6146, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 10556, 0, 3, 5696, 9926, 2788, 2824,
                                                 6176, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10646, 0, 3, 5726, 9971, 2824, 2860,
                                                 6236, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10736, 0, 3, 5756, 10016, 2860, 2896,
                                                 6296, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10826, 0, 3, 5786, 10061, 2896, 2932,
                                                 6356, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 10916, 0, 3, 5816, 10106, 2932, 2968,
                                                 6416, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11006, 0, 3, 5846, 10151, 2968, 3004,
                                                 6476, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11096, 0, 3, 5876, 10196, 3004, 3040,
                                                 6536, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11186, 0, 3, 5936, 10241, 3112, 3148,
                                                 6596, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11276, 0, 3, 5966, 10286, 3148, 3184,
                                                 6656, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11366, 0, 3, 5996, 10331, 3184, 3220,
                                                 6716, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11456, 0, 3, 6026, 10376, 3220, 3256,
                                                 6776, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11546, 0, 3, 6056, 10421, 3256, 3292,
                                                 6836, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11636, 0, 3, 6086, 10466, 3292, 3328,
                                                 6896, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 11726, 0, 3, 6116, 10511, 3328, 3364,
                                                 6956, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 11816, 0, 3, 6236, 10736, 3436, 3496,
                                                 7116, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 11966, 0, 3, 6296, 10826, 3496, 3556,
                                                 7216, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12116, 0, 3, 6356, 10916, 3556, 3616,
                                                 7316, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12266, 0, 3, 6416, 11006, 3616, 3676,
                                                 7416, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12416, 0, 3, 6476, 11096, 3676, 3736,
                                                 7516, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12566, 0, 3, 6656, 11366, 3856, 3916,
                                                 7716, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12716, 0, 3, 6716, 11456, 3916, 3976,
                                                 7816, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 12866, 0, 3, 6776, 11546, 3976, 4036,
                                                 7916, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13016, 0, 3, 6836, 11636, 4036, 4096,
                                                 8016, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 13166, 0, 3, 6896, 11726, 4096, 4156,
                                                 8116, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13316, 0, 3, 10556, 10646, 7016, 11816,
                                                 4276, 4366, 8216, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13541, 0, 3, 10646, 10736, 7116, 11966,
                                                 4366, 4456, 8366, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13766, 0, 3, 10736, 10826, 7216, 12116,
                                                 4456, 4546, 8516, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 13991, 0, 3, 10826, 10916, 7316, 12266,
                                                 4546, 4636, 8666, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14216, 0, 3, 10916, 11006, 7416, 12416,
                                                 4636, 4726, 8816, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14441, 0, 3, 11186, 11276, 7616, 12566,
                                                 4906, 4996, 8966, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14666, 0, 3, 11276, 11366, 7716, 12716,
                                                 4996, 5086, 9116, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 14891, 0, 3, 11366, 11456, 7816, 12866,
                                                 5086, 5176, 9266, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15116, 0, 3, 11456, 11546, 7916, 13016,
                                                 5176, 5266, 9416, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 15341, 0, 3, 11546, 11636, 8016, 13166,
                                                 5266, 5356, 9566, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15566, 3, 5536, 5546, 9731, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15587, 3, 5546, 5556, 9746, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15608, 3, 5556, 5566, 9761, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15629, 3, 5566, 5576, 9776, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15650, 3, 5576, 5586, 9791, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15671, 3, 5586, 5596, 9806, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15692, 3, 5616, 5626, 9836, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15713, 3, 5626, 5636, 9851, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15734, 3, 5636, 5646, 9866, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15755, 3, 5646, 5656, 9881, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15776, 3, 5656, 5666, 9896, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 15797, 3, 5666, 5676, 9911, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 15818, 0, 3, 9716, 15566, 9971, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 15881, 0, 3, 9731, 15587, 10016, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 15944, 0, 3, 9746, 15608, 10061, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16007, 0, 3, 9761, 15629, 10106, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16070, 0, 3, 9776, 15650, 10151, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16133, 0, 3, 9791, 15671, 10196, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16196, 0, 3, 9821, 15692, 10286, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16259, 0, 3, 9836, 15713, 10331, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16322, 0, 3, 9851, 15734, 10376, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16385, 0, 3, 9866, 15755, 10421, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16448, 0, 3, 9881, 15776, 10466, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 16511, 0, 3, 9896, 15797, 10511, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 16574, 0, 3, 9971, 15881, 6176, 6236,
                                                 10736, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 16700, 0, 3, 10016, 15944, 6236, 6296,
                                                 10826, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 16826, 0, 3, 10061, 16007, 6296, 6356,
                                                 10916, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 16952, 0, 3, 10106, 16070, 6356, 6416,
                                                 11006, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17078, 0, 3, 10151, 16133, 6416, 6476,
                                                 11096, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17204, 0, 3, 10286, 16259, 6596, 6656,
                                                 11366, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17330, 0, 3, 10331, 16322, 6656, 6716,
                                                 11456, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17456, 0, 3, 10376, 16385, 6716, 6776,
                                                 11546, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17582, 0, 3, 10421, 16448, 6776, 6836,
                                                 11636, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 17708, 0, 3, 10466, 16511, 6836, 6896,
                                                 11726, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 17834, 0, 3, 10736, 16700, 7016, 7116,
                                                 11966, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18044, 0, 3, 10826, 16826, 7116, 7216,
                                                 12116, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18254, 0, 3, 10916, 16952, 7216, 7316,
                                                 12266, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18464, 0, 3, 11006, 17078, 7316, 7416,
                                                 12416, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18674, 0, 3, 11366, 17330, 7616, 7716,
                                                 12716, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 18884, 0, 3, 11456, 17456, 7716, 7816,
                                                 12866, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 19094, 0, 3, 11546, 17582, 7816, 7916,
                                                 13016, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 19304, 0, 3, 11636, 17708, 7916, 8016,
                                                 13166, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 19514, 0, 3, 16574, 16700, 11966, 18044,
                                                 8216, 8366, 13766, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 19829, 0, 3, 16700, 16826, 12116, 18254,
                                                 8366, 8516, 13991, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 20144, 0, 3, 16826, 16952, 12266, 18464,
                                                 8516, 8666, 14216, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 20459, 0, 3, 17204, 17330, 12716, 18884,
                                                 8966, 9116, 14891, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 20774, 0, 3, 17330, 17456, 12866, 19094,
                                                 9116, 9266, 15116, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 21089, 0, 3, 17456, 17582, 13016, 19304,
                                                 9266, 9416, 15341, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21404, 3, 9716, 9731, 15587, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21432, 3, 9731, 9746, 15608, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21460, 3, 9746, 9761, 15629, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21488, 3, 9761, 9776, 15650, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21516, 3, 9776, 9791, 15671, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21544, 3, 9821, 9836, 15713, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21572, 3, 9836, 9851, 15734, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21600, 3, 9851, 9866, 15755, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21628, 3, 9866, 9881, 15776, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 21656, 3, 9881, 9896, 15797, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 21684, 0, 3, 15566, 21404, 15881, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 21768, 0, 3, 15587, 21432, 15944, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 21852, 0, 3, 15608, 21460, 16007, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 21936, 0, 3, 15629, 21488, 16070, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 22020, 0, 3, 15650, 21516, 16133, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 22104, 0, 3, 15692, 21544, 16259, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 22188, 0, 3, 15713, 21572, 16322, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 22272, 0, 3, 15734, 21600, 16385, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 22356, 0, 3, 15755, 21628, 16448, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 22440, 0, 3, 15776, 21656, 16511, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 22524, 0, 3, 15818, 21684, 10556, 10646,
                                                 16574, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 22692, 0, 3, 15881, 21768, 10646, 10736,
                                                 16700, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 22860, 0, 3, 15944, 21852, 10736, 10826,
                                                 16826, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23028, 0, 3, 16007, 21936, 10826, 10916,
                                                 16952, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23196, 0, 3, 16070, 22020, 10916, 11006,
                                                 17078, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23364, 0, 3, 16196, 22104, 11186, 11276,
                                                 17204, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23532, 0, 3, 16259, 22188, 11276, 11366,
                                                 17330, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23700, 0, 3, 16322, 22272, 11366, 11456,
                                                 17456, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 23868, 0, 3, 16385, 22356, 11456, 11546,
                                                 17582, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 24036, 0, 3, 16448, 22440, 11546, 11636,
                                                 17708, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 24204, 0, 3, 16700, 22860, 11816, 11966,
                                                 18044, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 24484, 0, 3, 16826, 23028, 11966, 12116,
                                                 18254, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 24764, 0, 3, 16952, 23196, 12116, 12266,
                                                 18464, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25044, 0, 3, 17330, 23700, 12566, 12716,
                                                 18884, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25324, 0, 3, 17456, 23868, 12716, 12866,
                                                 19094, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 25604, 0, 3, 17582, 24036, 12866, 13016,
                                                 19304, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 25884, 0, 3, 22524, 22692, 17834, 24204,
                                                 13316, 13541, 19514, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 26304, 0, 3, 22692, 22860, 18044, 24484,
                                                 13541, 13766, 19829, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 26724, 0, 3, 22860, 23028, 18254, 24764,
                                                 13766, 13991, 20144, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 27144, 0, 3, 23364, 23532, 18674, 25044,
                                                 14441, 14666, 20459, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 27564, 0, 3, 23532, 23700, 18884, 25324,
                                                 14666, 14891, 20774, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 27984, 0, 3, 23700, 23868, 19094, 25604,
                                                 14891, 15116, 21089, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 28404, 3, 15566, 15587, 21432, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 28440, 3, 15587, 15608, 21460, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 28476, 3, 15608, 15629, 21488, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 28512, 3, 15629, 15650, 21516, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 28548, 3, 15692, 15713, 21572, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 28584, 3, 15713, 15734, 21600, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 28620, 3, 15734, 15755, 21628, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 28656, 3, 15755, 15776, 21656, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 28692, 0, 3, 21404, 28404, 21768, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 28800, 0, 3, 21432, 28440, 21852, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 28908, 0, 3, 21460, 28476, 21936, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 29016, 0, 3, 21488, 28512, 22020, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 29124, 0, 3, 21544, 28548, 22188, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 29232, 0, 3, 21572, 28584, 22272, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 29340, 0, 3, 21600, 28620, 22356, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 29448, 0, 3, 21628, 28656, 22440, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 29556, 0, 3, 21768, 28800, 16574, 16700,
                                                 22860, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 29772, 0, 3, 21852, 28908, 16700, 16826,
                                                 23028, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 29988, 0, 3, 21936, 29016, 16826, 16952,
                                                 23196, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 30204, 0, 3, 22188, 29232, 17204, 17330,
                                                 23700, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 30420, 0, 3, 22272, 29340, 17330, 17456,
                                                 23868, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 30636, 0, 3, 22356, 29448, 17456, 17582,
                                                 24036, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 30852, 0, 3, 22860, 29772, 17834, 18044,
                                                 24484, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 31212, 0, 3, 23028, 29988, 18044, 18254,
                                                 24764, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 31572, 0, 3, 23700, 30420, 18674, 18884,
                                                 25324, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 31932, 0, 3, 23868, 30636, 18884, 19094,
                                                 25604, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 32292, 0, 3, 29556, 29772, 24484, 31212,
                                                 19514, 19829, 26724, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 32832, 0, 3, 30204, 30420, 25324, 31932,
                                                 20459, 20774, 27984, ncols, alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 33372, 3, 21404, 21432, 28440, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 33417, 3, 21432, 21460, 28476, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 33462, 3, 21460, 21488, 28512, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 33507, 3, 21544, 21572, 28584, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 33552, 3, 21572, 21600, 28620, ncols,
                                                 alpha, beta, p);

            compute_prim_sl_electron_repulsion_0(buffer, 33597, 3, 21600, 21628, 28656, ncols,
                                                 alpha, beta, p);

            compute_prim_pl_electron_repulsion_0(buffer, 33642, 0, 3, 28404, 33372, 28800, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 33777, 0, 3, 28440, 33417, 28908, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 33912, 0, 3, 28476, 33462, 29016, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 34047, 0, 3, 28548, 33507, 29232, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 34182, 0, 3, 28584, 33552, 29340, ncols,
                                                 p);

            compute_prim_pl_electron_repulsion_0(buffer, 34317, 0, 3, 28620, 33597, 29448, ncols,
                                                 p);

            compute_prim_dl_electron_repulsion_0(buffer, 34452, 0, 3, 28692, 33642, 22524, 22692,
                                                 29556, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 34722, 0, 3, 28800, 33777, 22692, 22860,
                                                 29772, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 34992, 0, 3, 28908, 33912, 22860, 23028,
                                                 29988, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 35262, 0, 3, 29124, 34047, 23364, 23532,
                                                 30204, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 35532, 0, 3, 29232, 34182, 23532, 23700,
                                                 30420, ncols, alpha, beta, p);

            compute_prim_dl_electron_repulsion_0(buffer, 35802, 0, 3, 29340, 34317, 23700, 23868,
                                                 30636, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 36072, 0, 3, 29772, 34992, 24204, 24484,
                                                 31212, ncols, alpha, beta, p);

            compute_prim_fl_electron_repulsion_0(buffer, 36522, 0, 3, 30420, 35802, 25044, 25324,
                                                 31932, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 36972, 0, 3, 34452, 34722, 30852, 36072,
                                                 25884, 26304, 32292, ncols, alpha, beta, p);

            compute_prim_gl_electron_repulsion_0(buffer, 37647, 0, 3, 35262, 35532, 31572, 36522,
                                                 27144, 27564, 32832, ncols, alpha, beta, p);

            simdgeo::geom_f_x(buffer, 38322, 35262, 37647, 1, 45, ncols, alpha);

            simdgeo::geom_f_y(buffer, 38772, 35262, 37647, 1, 45, ncols, alpha);

            simdgeo::geom_f_z(buffer, 39222, 35262, 37647, 1, 45, ncols, alpha);

            simdgeo::geom_f_x(buffer, 39672, 34452, 36972, 1, 45, ncols, alpha);

            simdgeo::geom_f_y(buffer, 40122, 34452, 36972, 1, 45, ncols, alpha);

            simdgeo::geom_f_z(buffer, 40572, 34452, 36972, 1, 45, ncols, alpha);

            simdfunc::contract_primitives(buffer, 41022, 39672, 1350, ncols);

            simdfunc::contract_primitives(buffer, 42372, 38322, 1350, ncols);
        }
    }

    simdtrf::transform_l_inner(buffer, 43722, 42372, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 43722, 17, nmax);

    simdtrf::transform_l_inner(buffer, 43722, 42822, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 119 * nvalues, nvalues, buffer, 43722, 17, nmax);

    simdtrf::transform_l_inner(buffer, 43722, 43272, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 238 * nvalues, nvalues, buffer, 43722, 17, nmax);

    simdtrf::transform_l_inner(buffer, 43722, 41022, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 357 * nvalues, nvalues, buffer, 43722, 17, nmax);

    simdtrf::transform_l_inner(buffer, 43722, 41472, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 476 * nvalues, nvalues, buffer, 43722, 17, nmax);

    simdtrf::transform_l_inner(buffer, 43722, 41922, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 595 * nvalues, nvalues, buffer, 43722, 17, nmax);
}

}  // namespace simdt2ceri
