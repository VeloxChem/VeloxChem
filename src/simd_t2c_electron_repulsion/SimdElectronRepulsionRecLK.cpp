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


#include "SimdElectronRepulsionRecLK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

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
#include "SimdTransformK.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_lk_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_lk_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(83013, nvalues);

    buffer.zero();

    const auto nmax = nvalues;

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

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            10, 11, 12, 13, 14, 15}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 61, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 7, 8, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 8, 9, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 9, 10, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 10, 11, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 11, 12, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 12, 13, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 13, 14, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 14, 15, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 15, 16, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 118, 0, 16, 17, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 124, 0, 17, 18, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 130, 0, 18, 19, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 136, 0, 19, 20, 61, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 142, 0, 22, 25, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 152, 0, 25, 28, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 162, 0, 28, 31, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 172, 0, 31, 34, 88, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 182, 0, 34, 37, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 192, 0, 37, 40, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 202, 0, 40, 43, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 212, 0, 43, 46, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 222, 0, 46, 49, 118, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 232, 0, 49, 52, 124, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 242, 0, 52, 55, 130, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 252, 0, 55, 58, 136, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 262, 0, 64, 70, 152, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 277, 0, 70, 76, 162, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 292, 0, 76, 82, 172, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 307, 0, 82, 88, 182, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 322, 0, 88, 94, 192, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 337, 0, 94, 100, 202, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 352, 0, 100, 106, 212, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 367, 0, 106, 112, 222, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 382, 0, 112, 118, 232, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 397, 0, 118, 124, 242, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 412, 0, 124, 130, 252, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 427, 0, 142, 152, 277, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 448, 0, 152, 162, 292, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 469, 0, 162, 172, 307, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 490, 0, 172, 182, 322, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 511, 0, 182, 192, 337, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 532, 0, 192, 202, 352, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 553, 0, 202, 212, 367, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 574, 0, 212, 222, 382, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 595, 0, 222, 232, 397, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 616, 0, 232, 242, 412, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 637, 0, 262, 277, 448, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 665, 0, 277, 292, 469, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 693, 0, 292, 307, 490, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 721, 0, 307, 322, 511, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 749, 0, 322, 337, 532, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 777, 0, 337, 352, 553, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 805, 0, 352, 367, 574, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 833, 0, 367, 382, 595, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 861, 0, 382, 397, 616, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 889, 0, 427, 448, 665, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 925, 0, 448, 469, 693, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 961, 0, 469, 490, 721, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 997, 0, 490, 511, 749, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1033, 0, 511, 532, 777, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1069, 0, 532, 553, 805, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1105, 0, 553, 574, 833, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1141, 0, 574, 595, 861, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1177, 0, 637, 665, 925, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1222, 0, 665, 693, 961, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1267, 0, 693, 721, 997, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1312, 0, 721, 749, 1033, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1357, 0, 749, 777, 1069, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1402, 0, 777, 805, 1105, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1447, 0, 805, 833, 1141, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1492, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1495, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1498, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1501, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1504, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1507, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1510, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1513, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1516, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1519, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1522, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1525, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1528, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1531, 3, 8, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1540, 3, 9, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1549, 3, 10, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1558, 3, 11, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1567, 3, 12, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1576, 3, 13, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1585, 3, 14, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1594, 3, 15, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1603, 3, 16, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1612, 3, 17, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1621, 3, 18, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1630, 3, 19, 58, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1639, 3, 20, 61, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1648, 0, 3, 22, 1531, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1666, 0, 3, 25, 1540, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1684, 0, 3, 28, 1549, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1702, 0, 3, 31, 1558, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1720, 0, 3, 34, 1567, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1738, 0, 3, 37, 1576, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1756, 0, 3, 40, 1585, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1774, 0, 3, 43, 1594, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1792, 0, 3, 46, 1603, 112, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1810, 0, 3, 49, 1612, 118, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1828, 0, 3, 52, 1621, 124, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1846, 0, 3, 55, 1630, 130, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1864, 0, 3, 58, 1639, 136, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1882, 0, 3, 70, 1684, 152, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1912, 0, 3, 76, 1702, 162, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1942, 0, 3, 82, 1720, 172, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1972, 0, 3, 88, 1738, 182, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2002, 0, 3, 94, 1756, 192, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2032, 0, 3, 100, 1774, 202, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2062, 0, 3, 106, 1792, 212, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2092, 0, 3, 112, 1810, 222, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2122, 0, 3, 118, 1828, 232, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2152, 0, 3, 124, 1846, 242, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2182, 0, 3, 130, 1864, 252, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2212, 0, 3, 142, 1882, 262, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2257, 0, 3, 152, 1912, 277, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2302, 0, 3, 162, 1942, 292, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2347, 0, 3, 172, 1972, 307, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2392, 0, 3, 182, 2002, 322, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2437, 0, 3, 192, 2032, 337, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2482, 0, 3, 202, 2062, 352, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2527, 0, 3, 212, 2092, 367, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2572, 0, 3, 222, 2122, 382, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2617, 0, 3, 232, 2152, 397, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2662, 0, 3, 242, 2182, 412, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2707, 0, 3, 277, 2302, 448, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2770, 0, 3, 292, 2347, 469, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2833, 0, 3, 307, 2392, 490, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2896, 0, 3, 322, 2437, 511, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2959, 0, 3, 337, 2482, 532, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3022, 0, 3, 352, 2527, 553, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3085, 0, 3, 367, 2572, 574, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3148, 0, 3, 382, 2617, 595, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 3211, 0, 3, 397, 2662, 616, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3274, 0, 3, 427, 2707, 637, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3358, 0, 3, 448, 2770, 665, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3442, 0, 3, 469, 2833, 693, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3526, 0, 3, 490, 2896, 721, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3610, 0, 3, 511, 2959, 749, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3694, 0, 3, 532, 3022, 777, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3778, 0, 3, 553, 3085, 805, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3862, 0, 3, 574, 3148, 833, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3946, 0, 3, 595, 3211, 861, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4030, 0, 3, 665, 3442, 925, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4138, 0, 3, 693, 3526, 961, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4246, 0, 3, 721, 3610, 997, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4354, 0, 3, 749, 3694, 1033, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4462, 0, 3, 777, 3778, 1069, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4570, 0, 3, 805, 3862, 1105, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4678, 0, 3, 833, 3946, 1141, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4786, 0, 3, 889, 4030, 1177, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4921, 0, 3, 925, 4138, 1222, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 5056, 0, 3, 961, 4246, 1267, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 5191, 0, 3, 997, 4354, 1312, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 5326, 0, 3, 1033, 4462, 1357, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5461, 0, 3, 1069, 4570, 1402, ncols,
                                                 p);

            compute_prim_lp_electron_repulsion_0(buffer, 5596, 0, 3, 1105, 4678, 1447, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 5731, 3, 8, 9, 1495, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 5737, 3, 9, 10, 1498, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5743, 3, 10, 11, 1501, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5749, 3, 11, 12, 1504, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5755, 3, 12, 13, 1507, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5761, 3, 13, 14, 1510, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5767, 3, 14, 15, 1513, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5773, 3, 15, 16, 1516, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5779, 3, 16, 17, 1519, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5785, 3, 17, 18, 1522, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5791, 3, 18, 19, 1525, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 5797, 3, 19, 20, 1528, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 5803, 0, 3, 1492, 5731, 1540, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5821, 0, 3, 1495, 5737, 1549, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5839, 0, 3, 1498, 5743, 1558, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5857, 0, 3, 1501, 5749, 1567, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5875, 0, 3, 1504, 5755, 1576, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5893, 0, 3, 1507, 5761, 1585, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5911, 0, 3, 1510, 5767, 1594, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5929, 0, 3, 1513, 5773, 1603, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5947, 0, 3, 1516, 5779, 1612, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5965, 0, 3, 1519, 5785, 1621, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5983, 0, 3, 1522, 5791, 1630, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 6001, 0, 3, 1525, 5797, 1639, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 6019, 0, 3, 1540, 5821, 64, 70, 1684,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6055, 0, 3, 1549, 5839, 70, 76, 1702,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6091, 0, 3, 1558, 5857, 76, 82, 1720,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6127, 0, 3, 1567, 5875, 82, 88, 1738,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6163, 0, 3, 1576, 5893, 88, 94, 1756,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6199, 0, 3, 1585, 5911, 94, 100, 1774,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6235, 0, 3, 1594, 5929, 100, 106, 1792,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6271, 0, 3, 1603, 5947, 106, 112, 1810,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6307, 0, 3, 1612, 5965, 112, 118, 1828,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6343, 0, 3, 1621, 5983, 118, 124, 1846,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 6379, 0, 3, 1630, 6001, 124, 130, 1864,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6415, 0, 3, 1684, 6055, 142, 152, 1912,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6475, 0, 3, 1702, 6091, 152, 162, 1942,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6535, 0, 3, 1720, 6127, 162, 172, 1972,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6595, 0, 3, 1738, 6163, 172, 182, 2002,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6655, 0, 3, 1756, 6199, 182, 192, 2032,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6715, 0, 3, 1774, 6235, 192, 202, 2062,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6775, 0, 3, 1792, 6271, 202, 212, 2092,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6835, 0, 3, 1810, 6307, 212, 222, 2122,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6895, 0, 3, 1828, 6343, 222, 232, 2152,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6955, 0, 3, 1846, 6379, 232, 242, 2182,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7015, 0, 3, 6019, 6055, 1912, 6475, 262,
                                                 277, 2302, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7105, 0, 3, 6055, 6091, 1942, 6535, 277,
                                                 292, 2347, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7195, 0, 3, 6091, 6127, 1972, 6595, 292,
                                                 307, 2392, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7285, 0, 3, 6127, 6163, 2002, 6655, 307,
                                                 322, 2437, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7375, 0, 3, 6163, 6199, 2032, 6715, 322,
                                                 337, 2482, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7465, 0, 3, 6199, 6235, 2062, 6775, 337,
                                                 352, 2527, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7555, 0, 3, 6235, 6271, 2092, 6835, 352,
                                                 367, 2572, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7645, 0, 3, 6271, 6307, 2122, 6895, 367,
                                                 382, 2617, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 7735, 0, 3, 6307, 6343, 2152, 6955, 382,
                                                 397, 2662, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7825, 0, 3, 6415, 6475, 2302, 7105, 427,
                                                 448, 2770, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7951, 0, 3, 6475, 6535, 2347, 7195, 448,
                                                 469, 2833, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8077, 0, 3, 6535, 6595, 2392, 7285, 469,
                                                 490, 2896, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8203, 0, 3, 6595, 6655, 2437, 7375, 490,
                                                 511, 2959, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8329, 0, 3, 6655, 6715, 2482, 7465, 511,
                                                 532, 3022, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8455, 0, 3, 6715, 6775, 2527, 7555, 532,
                                                 553, 3085, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8581, 0, 3, 6775, 6835, 2572, 7645, 553,
                                                 574, 3148, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 8707, 0, 3, 6835, 6895, 2617, 7735, 574,
                                                 595, 3211, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8833, 0, 3, 7015, 7105, 2770, 7951, 637,
                                                 665, 3442, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9001, 0, 3, 7105, 7195, 2833, 8077, 665,
                                                 693, 3526, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9169, 0, 3, 7195, 7285, 2896, 8203, 693,
                                                 721, 3610, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9337, 0, 3, 7285, 7375, 2959, 8329, 721,
                                                 749, 3694, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9505, 0, 3, 7375, 7465, 3022, 8455, 749,
                                                 777, 3778, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9673, 0, 3, 7465, 7555, 3085, 8581, 777,
                                                 805, 3862, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 9841, 0, 3, 7555, 7645, 3148, 8707, 805,
                                                 833, 3946, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10009, 0, 3, 7825, 7951, 3442, 9001,
                                                 889, 925, 4138, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10225, 0, 3, 7951, 8077, 3526, 9169,
                                                 925, 961, 4246, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10441, 0, 3, 8077, 8203, 3610, 9337,
                                                 961, 997, 4354, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10657, 0, 3, 8203, 8329, 3694, 9505,
                                                 997, 1033, 4462, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 10873, 0, 3, 8329, 8455, 3778, 9673,
                                                 1033, 1069, 4570, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 11089, 0, 3, 8455, 8581, 3862, 9841,
                                                 1069, 1105, 4678, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11305, 0, 3, 8833, 9001, 4138, 10225,
                                                 1177, 1222, 5056, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11575, 0, 3, 9001, 9169, 4246, 10441,
                                                 1222, 1267, 5191, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11845, 0, 3, 9169, 9337, 4354, 10657,
                                                 1267, 1312, 5326, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 12115, 0, 3, 9337, 9505, 4462, 10873,
                                                 1312, 1357, 5461, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 12385, 0, 3, 9505, 9673, 4570, 11089,
                                                 1357, 1402, 5596, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12655, 3, 1492, 1495, 5737, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12665, 3, 1495, 1498, 5743, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12675, 3, 1498, 1501, 5749, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12685, 3, 1501, 1504, 5755, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12695, 3, 1504, 1507, 5761, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12705, 3, 1507, 1510, 5767, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12715, 3, 1510, 1513, 5773, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12725, 3, 1513, 1516, 5779, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12735, 3, 1516, 1519, 5785, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12745, 3, 1519, 1522, 5791, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 12755, 3, 1522, 1525, 5797, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 12765, 0, 3, 5731, 12655, 5821, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12795, 0, 3, 5737, 12665, 5839, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12825, 0, 3, 5743, 12675, 5857, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12855, 0, 3, 5749, 12685, 5875, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12885, 0, 3, 5755, 12695, 5893, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12915, 0, 3, 5761, 12705, 5911, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12945, 0, 3, 5767, 12715, 5929, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 12975, 0, 3, 5773, 12725, 5947, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13005, 0, 3, 5779, 12735, 5965, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13035, 0, 3, 5785, 12745, 5983, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 13065, 0, 3, 5791, 12755, 6001, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 13095, 0, 3, 5803, 12765, 1648, 1666,
                                                 6019, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13155, 0, 3, 5821, 12795, 1666, 1684,
                                                 6055, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13215, 0, 3, 5839, 12825, 1684, 1702,
                                                 6091, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13275, 0, 3, 5857, 12855, 1702, 1720,
                                                 6127, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13335, 0, 3, 5875, 12885, 1720, 1738,
                                                 6163, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13395, 0, 3, 5893, 12915, 1738, 1756,
                                                 6199, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13455, 0, 3, 5911, 12945, 1756, 1774,
                                                 6235, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13515, 0, 3, 5929, 12975, 1774, 1792,
                                                 6271, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13575, 0, 3, 5947, 13005, 1792, 1810,
                                                 6307, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13635, 0, 3, 5965, 13035, 1810, 1828,
                                                 6343, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 13695, 0, 3, 5983, 13065, 1828, 1846,
                                                 6379, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13755, 0, 3, 6055, 13215, 1882, 1912,
                                                 6475, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13855, 0, 3, 6091, 13275, 1912, 1942,
                                                 6535, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13955, 0, 3, 6127, 13335, 1942, 1972,
                                                 6595, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14055, 0, 3, 6163, 13395, 1972, 2002,
                                                 6655, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14155, 0, 3, 6199, 13455, 2002, 2032,
                                                 6715, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14255, 0, 3, 6235, 13515, 2032, 2062,
                                                 6775, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14355, 0, 3, 6271, 13575, 2062, 2092,
                                                 6835, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14455, 0, 3, 6307, 13635, 2092, 2122,
                                                 6895, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 14555, 0, 3, 6343, 13695, 2122, 2152,
                                                 6955, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14655, 0, 3, 13095, 13155, 6415, 13755,
                                                 2212, 2257, 7015, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14805, 0, 3, 13155, 13215, 6475, 13855,
                                                 2257, 2302, 7105, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14955, 0, 3, 13215, 13275, 6535, 13955,
                                                 2302, 2347, 7195, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15105, 0, 3, 13275, 13335, 6595, 14055,
                                                 2347, 2392, 7285, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15255, 0, 3, 13335, 13395, 6655, 14155,
                                                 2392, 2437, 7375, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15405, 0, 3, 13395, 13455, 6715, 14255,
                                                 2437, 2482, 7465, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15555, 0, 3, 13455, 13515, 6775, 14355,
                                                 2482, 2527, 7555, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15705, 0, 3, 13515, 13575, 6835, 14455,
                                                 2527, 2572, 7645, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 15855, 0, 3, 13575, 13635, 6895, 14555,
                                                 2572, 2617, 7735, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16005, 0, 3, 13755, 13855, 7105, 14955,
                                                 2707, 2770, 7951, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16215, 0, 3, 13855, 13955, 7195, 15105,
                                                 2770, 2833, 8077, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16425, 0, 3, 13955, 14055, 7285, 15255,
                                                 2833, 2896, 8203, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16635, 0, 3, 14055, 14155, 7375, 15405,
                                                 2896, 2959, 8329, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 16845, 0, 3, 14155, 14255, 7465, 15555,
                                                 2959, 3022, 8455, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17055, 0, 3, 14255, 14355, 7555, 15705,
                                                 3022, 3085, 8581, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 17265, 0, 3, 14355, 14455, 7645, 15855,
                                                 3085, 3148, 8707, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17475, 0, 3, 14655, 14805, 7825, 16005,
                                                 3274, 3358, 8833, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 17755, 0, 3, 14805, 14955, 7951, 16215,
                                                 3358, 3442, 9001, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18035, 0, 3, 14955, 15105, 8077, 16425,
                                                 3442, 3526, 9169, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18315, 0, 3, 15105, 15255, 8203, 16635,
                                                 3526, 3610, 9337, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18595, 0, 3, 15255, 15405, 8329, 16845,
                                                 3610, 3694, 9505, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 18875, 0, 3, 15405, 15555, 8455, 17055,
                                                 3694, 3778, 9673, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 19155, 0, 3, 15555, 15705, 8581, 17265,
                                                 3778, 3862, 9841, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 19435, 0, 3, 16005, 16215, 9001, 18035,
                                                 4030, 4138, 10225, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 19795, 0, 3, 16215, 16425, 9169, 18315,
                                                 4138, 4246, 10441, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 20155, 0, 3, 16425, 16635, 9337, 18595,
                                                 4246, 4354, 10657, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 20515, 0, 3, 16635, 16845, 9505, 18875,
                                                 4354, 4462, 10873, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 20875, 0, 3, 16845, 17055, 9673, 19155,
                                                 4462, 4570, 11089, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 21235, 0, 3, 17475, 17755, 10009, 19435,
                                                 4786, 4921, 11305, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 21685, 0, 3, 17755, 18035, 10225, 19795,
                                                 4921, 5056, 11575, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 22135, 0, 3, 18035, 18315, 10441, 20155,
                                                 5056, 5191, 11845, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 22585, 0, 3, 18315, 18595, 10657, 20515,
                                                 5191, 5326, 12115, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 23035, 0, 3, 18595, 18875, 10873, 20875,
                                                 5326, 5461, 12385, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23485, 3, 5731, 5737, 12665, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23500, 3, 5737, 5743, 12675, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23515, 3, 5743, 5749, 12685, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23530, 3, 5749, 5755, 12695, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23545, 3, 5755, 5761, 12705, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23560, 3, 5761, 5767, 12715, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23575, 3, 5767, 5773, 12725, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23590, 3, 5773, 5779, 12735, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23605, 3, 5779, 5785, 12745, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 23620, 3, 5785, 5791, 12755, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 23635, 0, 3, 12655, 23485, 12795, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23680, 0, 3, 12665, 23500, 12825, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23725, 0, 3, 12675, 23515, 12855, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23770, 0, 3, 12685, 23530, 12885, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23815, 0, 3, 12695, 23545, 12915, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23860, 0, 3, 12705, 23560, 12945, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23905, 0, 3, 12715, 23575, 12975, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23950, 0, 3, 12725, 23590, 13005, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 23995, 0, 3, 12735, 23605, 13035, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 24040, 0, 3, 12745, 23620, 13065, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 24085, 0, 3, 12795, 23680, 6019, 6055,
                                                 13215, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24175, 0, 3, 12825, 23725, 6055, 6091,
                                                 13275, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24265, 0, 3, 12855, 23770, 6091, 6127,
                                                 13335, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24355, 0, 3, 12885, 23815, 6127, 6163,
                                                 13395, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24445, 0, 3, 12915, 23860, 6163, 6199,
                                                 13455, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24535, 0, 3, 12945, 23905, 6199, 6235,
                                                 13515, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24625, 0, 3, 12975, 23950, 6235, 6271,
                                                 13575, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24715, 0, 3, 13005, 23995, 6271, 6307,
                                                 13635, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 24805, 0, 3, 13035, 24040, 6307, 6343,
                                                 13695, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 24895, 0, 3, 13215, 24175, 6415, 6475,
                                                 13855, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25045, 0, 3, 13275, 24265, 6475, 6535,
                                                 13955, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25195, 0, 3, 13335, 24355, 6535, 6595,
                                                 14055, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25345, 0, 3, 13395, 24445, 6595, 6655,
                                                 14155, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25495, 0, 3, 13455, 24535, 6655, 6715,
                                                 14255, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25645, 0, 3, 13515, 24625, 6715, 6775,
                                                 14355, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25795, 0, 3, 13575, 24715, 6775, 6835,
                                                 14455, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 25945, 0, 3, 13635, 24805, 6835, 6895,
                                                 14555, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26095, 0, 3, 24085, 24175, 13855, 25045,
                                                 7015, 7105, 14955, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26320, 0, 3, 24175, 24265, 13955, 25195,
                                                 7105, 7195, 15105, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26545, 0, 3, 24265, 24355, 14055, 25345,
                                                 7195, 7285, 15255, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26770, 0, 3, 24355, 24445, 14155, 25495,
                                                 7285, 7375, 15405, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 26995, 0, 3, 24445, 24535, 14255, 25645,
                                                 7375, 7465, 15555, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27220, 0, 3, 24535, 24625, 14355, 25795,
                                                 7465, 7555, 15705, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 27445, 0, 3, 24625, 24715, 14455, 25945,
                                                 7555, 7645, 15855, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 27670, 0, 3, 24895, 25045, 14955, 26320,
                                                 7825, 7951, 16215, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 27985, 0, 3, 25045, 25195, 15105, 26545,
                                                 7951, 8077, 16425, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28300, 0, 3, 25195, 25345, 15255, 26770,
                                                 8077, 8203, 16635, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28615, 0, 3, 25345, 25495, 15405, 26995,
                                                 8203, 8329, 16845, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 28930, 0, 3, 25495, 25645, 15555, 27220,
                                                 8329, 8455, 17055, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 29245, 0, 3, 25645, 25795, 15705, 27445,
                                                 8455, 8581, 17265, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 29560, 0, 3, 26095, 26320, 16215, 27985,
                                                 8833, 9001, 18035, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 29980, 0, 3, 26320, 26545, 16425, 28300,
                                                 9001, 9169, 18315, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 30400, 0, 3, 26545, 26770, 16635, 28615,
                                                 9169, 9337, 18595, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 30820, 0, 3, 26770, 26995, 16845, 28930,
                                                 9337, 9505, 18875, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 31240, 0, 3, 26995, 27220, 17055, 29245,
                                                 9505, 9673, 19155, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 31660, 0, 3, 27670, 27985, 18035, 29980,
                                                 10009, 10225, 19795, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 32200, 0, 3, 27985, 28300, 18315, 30400,
                                                 10225, 10441, 20155, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 32740, 0, 3, 28300, 28615, 18595, 30820,
                                                 10441, 10657, 20515, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 33280, 0, 3, 28615, 28930, 18875, 31240,
                                                 10657, 10873, 20875, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 33820, 0, 3, 29560, 29980, 19795, 32200,
                                                 11305, 11575, 22135, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 34495, 0, 3, 29980, 30400, 20155, 32740,
                                                 11575, 11845, 22585, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 35170, 0, 3, 30400, 30820, 20515, 33280,
                                                 11845, 12115, 23035, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 35845, 3, 12655, 12665, 23500, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 35866, 3, 12665, 12675, 23515, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 35887, 3, 12675, 12685, 23530, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 35908, 3, 12685, 12695, 23545, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 35929, 3, 12695, 12705, 23560, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 35950, 3, 12705, 12715, 23575, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 35971, 3, 12715, 12725, 23590, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 35992, 3, 12725, 12735, 23605, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 36013, 3, 12735, 12745, 23620, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 36034, 0, 3, 23485, 35845, 23680, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 36097, 0, 3, 23500, 35866, 23725, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 36160, 0, 3, 23515, 35887, 23770, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 36223, 0, 3, 23530, 35908, 23815, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 36286, 0, 3, 23545, 35929, 23860, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 36349, 0, 3, 23560, 35950, 23905, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 36412, 0, 3, 23575, 35971, 23950, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 36475, 0, 3, 23590, 35992, 23995, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 36538, 0, 3, 23605, 36013, 24040, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 36601, 0, 3, 23635, 36034, 13095, 13155,
                                                 24085, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 36727, 0, 3, 23680, 36097, 13155, 13215,
                                                 24175, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 36853, 0, 3, 23725, 36160, 13215, 13275,
                                                 24265, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 36979, 0, 3, 23770, 36223, 13275, 13335,
                                                 24355, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 37105, 0, 3, 23815, 36286, 13335, 13395,
                                                 24445, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 37231, 0, 3, 23860, 36349, 13395, 13455,
                                                 24535, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 37357, 0, 3, 23905, 36412, 13455, 13515,
                                                 24625, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 37483, 0, 3, 23950, 36475, 13515, 13575,
                                                 24715, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 37609, 0, 3, 23995, 36538, 13575, 13635,
                                                 24805, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 37735, 0, 3, 24175, 36853, 13755, 13855,
                                                 25045, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 37945, 0, 3, 24265, 36979, 13855, 13955,
                                                 25195, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 38155, 0, 3, 24355, 37105, 13955, 14055,
                                                 25345, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 38365, 0, 3, 24445, 37231, 14055, 14155,
                                                 25495, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 38575, 0, 3, 24535, 37357, 14155, 14255,
                                                 25645, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 38785, 0, 3, 24625, 37483, 14255, 14355,
                                                 25795, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 38995, 0, 3, 24715, 37609, 14355, 14455,
                                                 25945, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 39205, 0, 3, 36601, 36727, 24895, 37735,
                                                 14655, 14805, 26095, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 39520, 0, 3, 36727, 36853, 25045, 37945,
                                                 14805, 14955, 26320, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 39835, 0, 3, 36853, 36979, 25195, 38155,
                                                 14955, 15105, 26545, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 40150, 0, 3, 36979, 37105, 25345, 38365,
                                                 15105, 15255, 26770, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 40465, 0, 3, 37105, 37231, 25495, 38575,
                                                 15255, 15405, 26995, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 40780, 0, 3, 37231, 37357, 25645, 38785,
                                                 15405, 15555, 27220, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 41095, 0, 3, 37357, 37483, 25795, 38995,
                                                 15555, 15705, 27445, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 41410, 0, 3, 37735, 37945, 26320, 39835,
                                                 16005, 16215, 27985, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 41851, 0, 3, 37945, 38155, 26545, 40150,
                                                 16215, 16425, 28300, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 42292, 0, 3, 38155, 38365, 26770, 40465,
                                                 16425, 16635, 28615, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 42733, 0, 3, 38365, 38575, 26995, 40780,
                                                 16635, 16845, 28930, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 43174, 0, 3, 38575, 38785, 27220, 41095,
                                                 16845, 17055, 29245, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 43615, 0, 3, 39205, 39520, 27670, 41410,
                                                 17475, 17755, 29560, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 44203, 0, 3, 39520, 39835, 27985, 41851,
                                                 17755, 18035, 29980, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 44791, 0, 3, 39835, 40150, 28300, 42292,
                                                 18035, 18315, 30400, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 45379, 0, 3, 40150, 40465, 28615, 42733,
                                                 18315, 18595, 30820, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 45967, 0, 3, 40465, 40780, 28930, 43174,
                                                 18595, 18875, 31240, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 46555, 0, 3, 41410, 41851, 29980, 44791,
                                                 19435, 19795, 32200, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 47311, 0, 3, 41851, 42292, 30400, 45379,
                                                 19795, 20155, 32740, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 48067, 0, 3, 42292, 42733, 30820, 45967,
                                                 20155, 20515, 33280, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 48823, 0, 3, 43615, 44203, 31660, 46555,
                                                 21235, 21685, 33820, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 49768, 0, 3, 44203, 44791, 32200, 47311,
                                                 21685, 22135, 34495, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 50713, 0, 3, 44791, 45379, 32740, 48067,
                                                 22135, 22585, 35170, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 51658, 3, 23485, 23500, 35866, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 51686, 3, 23500, 23515, 35887, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 51714, 3, 23515, 23530, 35908, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 51742, 3, 23530, 23545, 35929, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 51770, 3, 23545, 23560, 35950, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 51798, 3, 23560, 23575, 35971, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 51826, 3, 23575, 23590, 35992, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 51854, 3, 23590, 23605, 36013, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 51882, 0, 3, 35845, 51658, 36097, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 51966, 0, 3, 35866, 51686, 36160, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 52050, 0, 3, 35887, 51714, 36223, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 52134, 0, 3, 35908, 51742, 36286, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 52218, 0, 3, 35929, 51770, 36349, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 52302, 0, 3, 35950, 51798, 36412, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 52386, 0, 3, 35971, 51826, 36475, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 52470, 0, 3, 35992, 51854, 36538, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 52554, 0, 3, 36097, 51966, 24085, 24175,
                                                 36853, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 52722, 0, 3, 36160, 52050, 24175, 24265,
                                                 36979, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 52890, 0, 3, 36223, 52134, 24265, 24355,
                                                 37105, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 53058, 0, 3, 36286, 52218, 24355, 24445,
                                                 37231, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 53226, 0, 3, 36349, 52302, 24445, 24535,
                                                 37357, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 53394, 0, 3, 36412, 52386, 24535, 24625,
                                                 37483, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 53562, 0, 3, 36475, 52470, 24625, 24715,
                                                 37609, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 53730, 0, 3, 36853, 52722, 24895, 25045,
                                                 37945, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 54010, 0, 3, 36979, 52890, 25045, 25195,
                                                 38155, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 54290, 0, 3, 37105, 53058, 25195, 25345,
                                                 38365, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 54570, 0, 3, 37231, 53226, 25345, 25495,
                                                 38575, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 54850, 0, 3, 37357, 53394, 25495, 25645,
                                                 38785, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 55130, 0, 3, 37483, 53562, 25645, 25795,
                                                 38995, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 55410, 0, 3, 52554, 52722, 37945, 54010,
                                                 26095, 26320, 39835, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 55830, 0, 3, 52722, 52890, 38155, 54290,
                                                 26320, 26545, 40150, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 56250, 0, 3, 52890, 53058, 38365, 54570,
                                                 26545, 26770, 40465, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 56670, 0, 3, 53058, 53226, 38575, 54850,
                                                 26770, 26995, 40780, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 57090, 0, 3, 53226, 53394, 38785, 55130,
                                                 26995, 27220, 41095, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 57510, 0, 3, 53730, 54010, 39835, 55830,
                                                 27670, 27985, 41851, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 58098, 0, 3, 54010, 54290, 40150, 56250,
                                                 27985, 28300, 42292, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 58686, 0, 3, 54290, 54570, 40465, 56670,
                                                 28300, 28615, 42733, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 59274, 0, 3, 54570, 54850, 40780, 57090,
                                                 28615, 28930, 43174, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 59862, 0, 3, 55410, 55830, 41851, 58098,
                                                 29560, 29980, 44791, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 60646, 0, 3, 55830, 56250, 42292, 58686,
                                                 29980, 30400, 45379, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 61430, 0, 3, 56250, 56670, 42733, 59274,
                                                 30400, 30820, 45967, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 62214, 0, 3, 57510, 58098, 44791, 60646,
                                                 31660, 32200, 47311, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 63222, 0, 3, 58098, 58686, 45379, 61430,
                                                 32200, 32740, 48067, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 64230, 0, 3, 59862, 60646, 47311, 63222,
                                                 33820, 34495, 50713, ncols, alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 65490, 3, 35845, 35866, 51686, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 65526, 3, 35866, 35887, 51714, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 65562, 3, 35887, 35908, 51742, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 65598, 3, 35908, 35929, 51770, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 65634, 3, 35929, 35950, 51798, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 65670, 3, 35950, 35971, 51826, ncols,
                                                 alpha, beta, p);

            compute_prim_sk_electron_repulsion_0(buffer, 65706, 3, 35971, 35992, 51854, ncols,
                                                 alpha, beta, p);

            compute_prim_pk_electron_repulsion_0(buffer, 65742, 0, 3, 51658, 65490, 51966, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 65850, 0, 3, 51686, 65526, 52050, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 65958, 0, 3, 51714, 65562, 52134, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 66066, 0, 3, 51742, 65598, 52218, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 66174, 0, 3, 51770, 65634, 52302, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 66282, 0, 3, 51798, 65670, 52386, ncols,
                                                 p);

            compute_prim_pk_electron_repulsion_0(buffer, 66390, 0, 3, 51826, 65706, 52470, ncols,
                                                 p);

            compute_prim_dk_electron_repulsion_0(buffer, 66498, 0, 3, 51882, 65742, 36601, 36727,
                                                 52554, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 66714, 0, 3, 51966, 65850, 36727, 36853,
                                                 52722, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 66930, 0, 3, 52050, 65958, 36853, 36979,
                                                 52890, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 67146, 0, 3, 52134, 66066, 36979, 37105,
                                                 53058, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 67362, 0, 3, 52218, 66174, 37105, 37231,
                                                 53226, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 67578, 0, 3, 52302, 66282, 37231, 37357,
                                                 53394, ncols, alpha, beta, p);

            compute_prim_dk_electron_repulsion_0(buffer, 67794, 0, 3, 52386, 66390, 37357, 37483,
                                                 53562, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 68010, 0, 3, 52722, 66930, 37735, 37945,
                                                 54010, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 68370, 0, 3, 52890, 67146, 37945, 38155,
                                                 54290, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 68730, 0, 3, 53058, 67362, 38155, 38365,
                                                 54570, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 69090, 0, 3, 53226, 67578, 38365, 38575,
                                                 54850, ncols, alpha, beta, p);

            compute_prim_fk_electron_repulsion_0(buffer, 69450, 0, 3, 53394, 67794, 38575, 38785,
                                                 55130, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 69810, 0, 3, 66498, 66714, 53730, 68010,
                                                 39205, 39520, 55410, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 70350, 0, 3, 66714, 66930, 54010, 68370,
                                                 39520, 39835, 55830, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 70890, 0, 3, 66930, 67146, 54290, 68730,
                                                 39835, 40150, 56250, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 71430, 0, 3, 67146, 67362, 54570, 69090,
                                                 40150, 40465, 56670, ncols, alpha, beta, p);

            compute_prim_gk_electron_repulsion_0(buffer, 71970, 0, 3, 67362, 67578, 54850, 69450,
                                                 40465, 40780, 57090, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 72510, 0, 3, 68010, 68370, 55830, 70890,
                                                 41410, 41851, 58098, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 73266, 0, 3, 68370, 68730, 56250, 71430,
                                                 41851, 42292, 58686, ncols, alpha, beta, p);

            compute_prim_hk_electron_repulsion_0(buffer, 74022, 0, 3, 68730, 69090, 56670, 71970,
                                                 42292, 42733, 59274, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 74778, 0, 3, 69810, 70350, 57510, 72510,
                                                 43615, 44203, 59862, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 75786, 0, 3, 70350, 70890, 58098, 73266,
                                                 44203, 44791, 60646, ncols, alpha, beta, p);

            compute_prim_ik_electron_repulsion_0(buffer, 76794, 0, 3, 70890, 71430, 58686, 74022,
                                                 44791, 45379, 61430, ncols, alpha, beta, p);

            compute_prim_kk_electron_repulsion_0(buffer, 77802, 0, 3, 72510, 73266, 60646, 76794,
                                                 46555, 47311, 63222, ncols, alpha, beta, p);

            compute_prim_lk_electron_repulsion_0(buffer, 79098, 0, 3, 74778, 75786, 62214, 77802,
                                                 48823, 49768, 64230, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 80718, 79098, 1620, ncols);
        }
    }

    simdtrf::transform_k_inner(buffer, 82338, 80718, 45, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 82338, 15, nmax);
}

}  // namespace simdt2ceri
