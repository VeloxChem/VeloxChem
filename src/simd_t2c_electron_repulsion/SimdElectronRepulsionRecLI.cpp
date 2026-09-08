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


#include "SimdElectronRepulsionRecLI.hpp"

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
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_li_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_li_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(54980, nvalues);

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

            simdfunc::compute_full_boys_function(buffer, coordinates, 6, 14, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 58, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 61, 0, 7, 8, 22, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 67, 0, 8, 9, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 73, 0, 9, 10, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 79, 0, 10, 11, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 85, 0, 11, 12, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 91, 0, 12, 13, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 97, 0, 13, 14, 40, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 103, 0, 14, 15, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 109, 0, 15, 16, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 115, 0, 16, 17, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 121, 0, 17, 18, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 127, 0, 18, 19, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 133, 0, 19, 20, 58, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 139, 0, 22, 25, 73, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 149, 0, 25, 28, 79, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 159, 0, 28, 31, 85, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 169, 0, 31, 34, 91, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 179, 0, 34, 37, 97, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 189, 0, 37, 40, 103, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 199, 0, 40, 43, 109, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 209, 0, 43, 46, 115, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 219, 0, 46, 49, 121, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 229, 0, 49, 52, 127, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 239, 0, 52, 55, 133, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 249, 0, 61, 67, 139, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 264, 0, 67, 73, 149, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 279, 0, 73, 79, 159, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 294, 0, 79, 85, 169, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 309, 0, 85, 91, 179, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 324, 0, 91, 97, 189, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 339, 0, 97, 103, 199, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 354, 0, 103, 109, 209, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 369, 0, 109, 115, 219, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 384, 0, 115, 121, 229, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 399, 0, 121, 127, 239, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 414, 0, 139, 149, 279, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 435, 0, 149, 159, 294, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 456, 0, 159, 169, 309, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 477, 0, 169, 179, 324, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 498, 0, 179, 189, 339, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 519, 0, 189, 199, 354, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 540, 0, 199, 209, 369, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 561, 0, 209, 219, 384, ncols, alpha,
                                                 beta, p);

            compute_prim_hs_electron_repulsion_0(buffer, 582, 0, 219, 229, 399, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 603, 0, 249, 264, 414, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 631, 0, 264, 279, 435, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 659, 0, 279, 294, 456, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 687, 0, 294, 309, 477, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 715, 0, 309, 324, 498, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 743, 0, 324, 339, 519, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 771, 0, 339, 354, 540, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 799, 0, 354, 369, 561, ncols, alpha,
                                                 beta, p);

            compute_prim_is_electron_repulsion_0(buffer, 827, 0, 369, 384, 582, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 855, 0, 414, 435, 659, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 891, 0, 435, 456, 687, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 927, 0, 456, 477, 715, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 963, 0, 477, 498, 743, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 999, 0, 498, 519, 771, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1035, 0, 519, 540, 799, ncols, alpha,
                                                 beta, p);

            compute_prim_ks_electron_repulsion_0(buffer, 1071, 0, 540, 561, 827, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1107, 0, 603, 631, 855, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1152, 0, 631, 659, 891, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1197, 0, 659, 687, 927, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1242, 0, 687, 715, 963, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1287, 0, 715, 743, 999, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1332, 0, 743, 771, 1035, ncols, alpha,
                                                 beta, p);

            compute_prim_ls_electron_repulsion_0(buffer, 1377, 0, 771, 799, 1071, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 1422, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1425, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1428, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1431, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1434, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1437, 3, 15, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1440, 3, 16, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1443, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1446, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1449, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1452, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 1455, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 1458, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1467, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1476, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1485, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1494, 3, 13, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1503, 3, 14, 40, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1512, 3, 15, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1521, 3, 16, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1530, 3, 17, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1539, 3, 18, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1548, 3, 19, 55, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 1557, 3, 20, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1566, 0, 3, 25, 1467, 73, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1584, 0, 3, 28, 1476, 79, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1602, 0, 3, 31, 1485, 85, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1620, 0, 3, 34, 1494, 91, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1638, 0, 3, 37, 1503, 97, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1656, 0, 3, 40, 1512, 103, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1674, 0, 3, 43, 1521, 109, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1692, 0, 3, 46, 1530, 115, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1710, 0, 3, 49, 1539, 121, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1728, 0, 3, 52, 1548, 127, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 1746, 0, 3, 55, 1557, 133, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1764, 0, 3, 73, 1584, 149, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1794, 0, 3, 79, 1602, 159, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1824, 0, 3, 85, 1620, 169, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1854, 0, 3, 91, 1638, 179, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1884, 0, 3, 97, 1656, 189, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1914, 0, 3, 103, 1674, 199, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1944, 0, 3, 109, 1692, 209, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 1974, 0, 3, 115, 1710, 219, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2004, 0, 3, 121, 1728, 229, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 2034, 0, 3, 127, 1746, 239, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2064, 0, 3, 149, 1794, 279, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2109, 0, 3, 159, 1824, 294, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2154, 0, 3, 169, 1854, 309, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2199, 0, 3, 179, 1884, 324, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2244, 0, 3, 189, 1914, 339, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2289, 0, 3, 199, 1944, 354, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2334, 0, 3, 209, 1974, 369, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2379, 0, 3, 219, 2004, 384, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 2424, 0, 3, 229, 2034, 399, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2469, 0, 3, 279, 2109, 435, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2532, 0, 3, 294, 2154, 456, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2595, 0, 3, 309, 2199, 477, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2658, 0, 3, 324, 2244, 498, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2721, 0, 3, 339, 2289, 519, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2784, 0, 3, 354, 2334, 540, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2847, 0, 3, 369, 2379, 561, ncols, p);

            compute_prim_hp_electron_repulsion_0(buffer, 2910, 0, 3, 384, 2424, 582, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 2973, 0, 3, 435, 2532, 659, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3057, 0, 3, 456, 2595, 687, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3141, 0, 3, 477, 2658, 715, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3225, 0, 3, 498, 2721, 743, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3309, 0, 3, 519, 2784, 771, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3393, 0, 3, 540, 2847, 799, ncols, p);

            compute_prim_ip_electron_repulsion_0(buffer, 3477, 0, 3, 561, 2910, 827, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3561, 0, 3, 659, 3057, 891, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3669, 0, 3, 687, 3141, 927, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3777, 0, 3, 715, 3225, 963, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3885, 0, 3, 743, 3309, 999, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 3993, 0, 3, 771, 3393, 1035, ncols, p);

            compute_prim_kp_electron_repulsion_0(buffer, 4101, 0, 3, 799, 3477, 1071, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4209, 0, 3, 891, 3669, 1197, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4344, 0, 3, 927, 3777, 1242, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4479, 0, 3, 963, 3885, 1287, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4614, 0, 3, 999, 3993, 1332, ncols, p);

            compute_prim_lp_electron_repulsion_0(buffer, 4749, 0, 3, 1035, 4101, 1377, ncols,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 4884, 3, 9, 10, 1425, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4890, 3, 10, 11, 1428, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4896, 3, 11, 12, 1431, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4902, 3, 12, 13, 1434, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4908, 3, 13, 14, 1437, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4914, 3, 14, 15, 1440, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4920, 3, 15, 16, 1443, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4926, 3, 16, 17, 1446, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4932, 3, 17, 18, 1449, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4938, 3, 18, 19, 1452, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 4944, 3, 19, 20, 1455, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 4950, 0, 3, 1422, 4884, 1467, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4968, 0, 3, 1425, 4890, 1476, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 4986, 0, 3, 1428, 4896, 1485, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5004, 0, 3, 1431, 4902, 1494, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5022, 0, 3, 1434, 4908, 1503, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5040, 0, 3, 1437, 4914, 1512, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5058, 0, 3, 1440, 4920, 1521, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5076, 0, 3, 1443, 4926, 1530, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5094, 0, 3, 1446, 4932, 1539, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5112, 0, 3, 1449, 4938, 1548, ncols,
                                                 p);

            compute_prim_pd_electron_repulsion_0(buffer, 5130, 0, 3, 1452, 4944, 1557, ncols,
                                                 p);

            compute_prim_dd_electron_repulsion_0(buffer, 5148, 0, 3, 1458, 4950, 61, 67, 1566,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5184, 0, 3, 1467, 4968, 67, 73, 1584,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5220, 0, 3, 1476, 4986, 73, 79, 1602,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5256, 0, 3, 1485, 5004, 79, 85, 1620,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5292, 0, 3, 1494, 5022, 85, 91, 1638,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5328, 0, 3, 1503, 5040, 91, 97, 1656,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5364, 0, 3, 1512, 5058, 97, 103, 1674,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5400, 0, 3, 1521, 5076, 103, 109, 1692,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5436, 0, 3, 1530, 5094, 109, 115, 1710,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5472, 0, 3, 1539, 5112, 115, 121, 1728,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 5508, 0, 3, 1548, 5130, 121, 127, 1746,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5544, 0, 3, 1584, 5220, 139, 149, 1794,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5604, 0, 3, 1602, 5256, 149, 159, 1824,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5664, 0, 3, 1620, 5292, 159, 169, 1854,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5724, 0, 3, 1638, 5328, 169, 179, 1884,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5784, 0, 3, 1656, 5364, 179, 189, 1914,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5844, 0, 3, 1674, 5400, 189, 199, 1944,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5904, 0, 3, 1692, 5436, 199, 209, 1974,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 5964, 0, 3, 1710, 5472, 209, 219, 2004,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 6024, 0, 3, 1728, 5508, 219, 229, 2034,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6084, 0, 3, 5148, 5184, 1764, 5544, 249,
                                                 264, 2064, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6174, 0, 3, 5184, 5220, 1794, 5604, 264,
                                                 279, 2109, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6264, 0, 3, 5220, 5256, 1824, 5664, 279,
                                                 294, 2154, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6354, 0, 3, 5256, 5292, 1854, 5724, 294,
                                                 309, 2199, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6444, 0, 3, 5292, 5328, 1884, 5784, 309,
                                                 324, 2244, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6534, 0, 3, 5328, 5364, 1914, 5844, 324,
                                                 339, 2289, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6624, 0, 3, 5364, 5400, 1944, 5904, 339,
                                                 354, 2334, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6714, 0, 3, 5400, 5436, 1974, 5964, 354,
                                                 369, 2379, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 6804, 0, 3, 5436, 5472, 2004, 6024, 369,
                                                 384, 2424, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 6894, 0, 3, 5544, 5604, 2109, 6264, 414,
                                                 435, 2532, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7020, 0, 3, 5604, 5664, 2154, 6354, 435,
                                                 456, 2595, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7146, 0, 3, 5664, 5724, 2199, 6444, 456,
                                                 477, 2658, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7272, 0, 3, 5724, 5784, 2244, 6534, 477,
                                                 498, 2721, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7398, 0, 3, 5784, 5844, 2289, 6624, 498,
                                                 519, 2784, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7524, 0, 3, 5844, 5904, 2334, 6714, 519,
                                                 540, 2847, ncols, alpha, beta, p);

            compute_prim_hd_electron_repulsion_0(buffer, 7650, 0, 3, 5904, 5964, 2379, 6804, 540,
                                                 561, 2910, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7776, 0, 3, 6084, 6174, 2469, 6894, 603,
                                                 631, 2973, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 7944, 0, 3, 6174, 6264, 2532, 7020, 631,
                                                 659, 3057, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8112, 0, 3, 6264, 6354, 2595, 7146, 659,
                                                 687, 3141, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8280, 0, 3, 6354, 6444, 2658, 7272, 687,
                                                 715, 3225, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8448, 0, 3, 6444, 6534, 2721, 7398, 715,
                                                 743, 3309, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8616, 0, 3, 6534, 6624, 2784, 7524, 743,
                                                 771, 3393, ncols, alpha, beta, p);

            compute_prim_id_electron_repulsion_0(buffer, 8784, 0, 3, 6624, 6714, 2847, 7650, 771,
                                                 799, 3477, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 8952, 0, 3, 6894, 7020, 3057, 8112, 855,
                                                 891, 3669, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9168, 0, 3, 7020, 7146, 3141, 8280, 891,
                                                 927, 3777, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9384, 0, 3, 7146, 7272, 3225, 8448, 927,
                                                 963, 3885, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9600, 0, 3, 7272, 7398, 3309, 8616, 963,
                                                 999, 3993, ncols, alpha, beta, p);

            compute_prim_kd_electron_repulsion_0(buffer, 9816, 0, 3, 7398, 7524, 3393, 8784, 999,
                                                 1035, 4101, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 10032, 0, 3, 7776, 7944, 3561, 8952,
                                                 1107, 1152, 4209, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 10302, 0, 3, 7944, 8112, 3669, 9168,
                                                 1152, 1197, 4344, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 10572, 0, 3, 8112, 8280, 3777, 9384,
                                                 1197, 1242, 4479, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 10842, 0, 3, 8280, 8448, 3885, 9600,
                                                 1242, 1287, 4614, ncols, alpha, beta, p);

            compute_prim_ld_electron_repulsion_0(buffer, 11112, 0, 3, 8448, 8616, 3993, 9816,
                                                 1287, 1332, 4749, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11382, 3, 1422, 1425, 4890, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11392, 3, 1425, 1428, 4896, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11402, 3, 1428, 1431, 4902, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11412, 3, 1431, 1434, 4908, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11422, 3, 1434, 1437, 4914, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11432, 3, 1437, 1440, 4920, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11442, 3, 1440, 1443, 4926, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11452, 3, 1443, 1446, 4932, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11462, 3, 1446, 1449, 4938, ncols,
                                                 alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 11472, 3, 1449, 1452, 4944, ncols,
                                                 alpha, beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 11482, 0, 3, 4884, 11382, 4968, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11512, 0, 3, 4890, 11392, 4986, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11542, 0, 3, 4896, 11402, 5004, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11572, 0, 3, 4902, 11412, 5022, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11602, 0, 3, 4908, 11422, 5040, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11632, 0, 3, 4914, 11432, 5058, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11662, 0, 3, 4920, 11442, 5076, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11692, 0, 3, 4926, 11452, 5094, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11722, 0, 3, 4932, 11462, 5112, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 11752, 0, 3, 4938, 11472, 5130, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 11782, 0, 3, 4968, 11512, 1566, 1584,
                                                 5220, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11842, 0, 3, 4986, 11542, 1584, 1602,
                                                 5256, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11902, 0, 3, 5004, 11572, 1602, 1620,
                                                 5292, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 11962, 0, 3, 5022, 11602, 1620, 1638,
                                                 5328, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12022, 0, 3, 5040, 11632, 1638, 1656,
                                                 5364, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12082, 0, 3, 5058, 11662, 1656, 1674,
                                                 5400, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12142, 0, 3, 5076, 11692, 1674, 1692,
                                                 5436, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12202, 0, 3, 5094, 11722, 1692, 1710,
                                                 5472, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 12262, 0, 3, 5112, 11752, 1710, 1728,
                                                 5508, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12322, 0, 3, 5220, 11842, 1764, 1794,
                                                 5604, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12422, 0, 3, 5256, 11902, 1794, 1824,
                                                 5664, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12522, 0, 3, 5292, 11962, 1824, 1854,
                                                 5724, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12622, 0, 3, 5328, 12022, 1854, 1884,
                                                 5784, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12722, 0, 3, 5364, 12082, 1884, 1914,
                                                 5844, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12822, 0, 3, 5400, 12142, 1914, 1944,
                                                 5904, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 12922, 0, 3, 5436, 12202, 1944, 1974,
                                                 5964, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 13022, 0, 3, 5472, 12262, 1974, 2004,
                                                 6024, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13122, 0, 3, 11782, 11842, 5604, 12422,
                                                 2064, 2109, 6264, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13272, 0, 3, 11842, 11902, 5664, 12522,
                                                 2109, 2154, 6354, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13422, 0, 3, 11902, 11962, 5724, 12622,
                                                 2154, 2199, 6444, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13572, 0, 3, 11962, 12022, 5784, 12722,
                                                 2199, 2244, 6534, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13722, 0, 3, 12022, 12082, 5844, 12822,
                                                 2244, 2289, 6624, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 13872, 0, 3, 12082, 12142, 5904, 12922,
                                                 2289, 2334, 6714, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 14022, 0, 3, 12142, 12202, 5964, 13022,
                                                 2334, 2379, 6804, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14172, 0, 3, 12322, 12422, 6264, 13272,
                                                 2469, 2532, 7020, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14382, 0, 3, 12422, 12522, 6354, 13422,
                                                 2532, 2595, 7146, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14592, 0, 3, 12522, 12622, 6444, 13572,
                                                 2595, 2658, 7272, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 14802, 0, 3, 12622, 12722, 6534, 13722,
                                                 2658, 2721, 7398, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15012, 0, 3, 12722, 12822, 6624, 13872,
                                                 2721, 2784, 7524, ncols, alpha, beta, p);

            compute_prim_hf_electron_repulsion_0(buffer, 15222, 0, 3, 12822, 12922, 6714, 14022,
                                                 2784, 2847, 7650, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15432, 0, 3, 13122, 13272, 7020, 14382,
                                                 2973, 3057, 8112, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15712, 0, 3, 13272, 13422, 7146, 14592,
                                                 3057, 3141, 8280, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 15992, 0, 3, 13422, 13572, 7272, 14802,
                                                 3141, 3225, 8448, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16272, 0, 3, 13572, 13722, 7398, 15012,
                                                 3225, 3309, 8616, ncols, alpha, beta, p);

            compute_prim_if_electron_repulsion_0(buffer, 16552, 0, 3, 13722, 13872, 7524, 15222,
                                                 3309, 3393, 8784, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 16832, 0, 3, 14172, 14382, 8112, 15712,
                                                 3561, 3669, 9168, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 17192, 0, 3, 14382, 14592, 8280, 15992,
                                                 3669, 3777, 9384, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 17552, 0, 3, 14592, 14802, 8448, 16272,
                                                 3777, 3885, 9600, ncols, alpha, beta, p);

            compute_prim_kf_electron_repulsion_0(buffer, 17912, 0, 3, 14802, 15012, 8616, 16552,
                                                 3885, 3993, 9816, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 18272, 0, 3, 15432, 15712, 9168, 17192,
                                                 4209, 4344, 10572, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 18722, 0, 3, 15712, 15992, 9384, 17552,
                                                 4344, 4479, 10842, ncols, alpha, beta, p);

            compute_prim_lf_electron_repulsion_0(buffer, 19172, 0, 3, 15992, 16272, 9600, 17912,
                                                 4479, 4614, 11112, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19622, 3, 4884, 4890, 11392, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19637, 3, 4890, 4896, 11402, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19652, 3, 4896, 4902, 11412, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19667, 3, 4902, 4908, 11422, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19682, 3, 4908, 4914, 11432, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19697, 3, 4914, 4920, 11442, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19712, 3, 4920, 4926, 11452, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19727, 3, 4926, 4932, 11462, ncols,
                                                 alpha, beta, p);

            compute_prim_sg_electron_repulsion_0(buffer, 19742, 3, 4932, 4938, 11472, ncols,
                                                 alpha, beta, p);

            compute_prim_pg_electron_repulsion_0(buffer, 19757, 0, 3, 11382, 19622, 11512, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 19802, 0, 3, 11392, 19637, 11542, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 19847, 0, 3, 11402, 19652, 11572, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 19892, 0, 3, 11412, 19667, 11602, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 19937, 0, 3, 11422, 19682, 11632, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 19982, 0, 3, 11432, 19697, 11662, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 20027, 0, 3, 11442, 19712, 11692, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 20072, 0, 3, 11452, 19727, 11722, ncols,
                                                 p);

            compute_prim_pg_electron_repulsion_0(buffer, 20117, 0, 3, 11462, 19742, 11752, ncols,
                                                 p);

            compute_prim_dg_electron_repulsion_0(buffer, 20162, 0, 3, 11482, 19757, 5148, 5184,
                                                 11782, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20252, 0, 3, 11512, 19802, 5184, 5220,
                                                 11842, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20342, 0, 3, 11542, 19847, 5220, 5256,
                                                 11902, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20432, 0, 3, 11572, 19892, 5256, 5292,
                                                 11962, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20522, 0, 3, 11602, 19937, 5292, 5328,
                                                 12022, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20612, 0, 3, 11632, 19982, 5328, 5364,
                                                 12082, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20702, 0, 3, 11662, 20027, 5364, 5400,
                                                 12142, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20792, 0, 3, 11692, 20072, 5400, 5436,
                                                 12202, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_0(buffer, 20882, 0, 3, 11722, 20117, 5436, 5472,
                                                 12262, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 20972, 0, 3, 11842, 20342, 5544, 5604,
                                                 12422, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21122, 0, 3, 11902, 20432, 5604, 5664,
                                                 12522, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21272, 0, 3, 11962, 20522, 5664, 5724,
                                                 12622, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21422, 0, 3, 12022, 20612, 5724, 5784,
                                                 12722, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21572, 0, 3, 12082, 20702, 5784, 5844,
                                                 12822, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21722, 0, 3, 12142, 20792, 5844, 5904,
                                                 12922, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_0(buffer, 21872, 0, 3, 12202, 20882, 5904, 5964,
                                                 13022, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22022, 0, 3, 20162, 20252, 12322, 20972,
                                                 6084, 6174, 13122, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22247, 0, 3, 20252, 20342, 12422, 21122,
                                                 6174, 6264, 13272, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22472, 0, 3, 20342, 20432, 12522, 21272,
                                                 6264, 6354, 13422, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22697, 0, 3, 20432, 20522, 12622, 21422,
                                                 6354, 6444, 13572, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 22922, 0, 3, 20522, 20612, 12722, 21572,
                                                 6444, 6534, 13722, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 23147, 0, 3, 20612, 20702, 12822, 21722,
                                                 6534, 6624, 13872, ncols, alpha, beta, p);

            compute_prim_gg_electron_repulsion_0(buffer, 23372, 0, 3, 20702, 20792, 12922, 21872,
                                                 6624, 6714, 14022, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23597, 0, 3, 20972, 21122, 13272, 22472,
                                                 6894, 7020, 14382, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 23912, 0, 3, 21122, 21272, 13422, 22697,
                                                 7020, 7146, 14592, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24227, 0, 3, 21272, 21422, 13572, 22922,
                                                 7146, 7272, 14802, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24542, 0, 3, 21422, 21572, 13722, 23147,
                                                 7272, 7398, 15012, ncols, alpha, beta, p);

            compute_prim_hg_electron_repulsion_0(buffer, 24857, 0, 3, 21572, 21722, 13872, 23372,
                                                 7398, 7524, 15222, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 25172, 0, 3, 22022, 22247, 14172, 23597,
                                                 7776, 7944, 15432, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 25592, 0, 3, 22247, 22472, 14382, 23912,
                                                 7944, 8112, 15712, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 26012, 0, 3, 22472, 22697, 14592, 24227,
                                                 8112, 8280, 15992, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 26432, 0, 3, 22697, 22922, 14802, 24542,
                                                 8280, 8448, 16272, ncols, alpha, beta, p);

            compute_prim_ig_electron_repulsion_0(buffer, 26852, 0, 3, 22922, 23147, 15012, 24857,
                                                 8448, 8616, 16552, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 27272, 0, 3, 23597, 23912, 15712, 26012,
                                                 8952, 9168, 17192, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 27812, 0, 3, 23912, 24227, 15992, 26432,
                                                 9168, 9384, 17552, ncols, alpha, beta, p);

            compute_prim_kg_electron_repulsion_0(buffer, 28352, 0, 3, 24227, 24542, 16272, 26852,
                                                 9384, 9600, 17912, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 28892, 0, 3, 25172, 25592, 16832, 27272,
                                                 10032, 10302, 18272, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 29567, 0, 3, 25592, 26012, 17192, 27812,
                                                 10302, 10572, 18722, ncols, alpha, beta, p);

            compute_prim_lg_electron_repulsion_0(buffer, 30242, 0, 3, 26012, 26432, 17552, 28352,
                                                 10572, 10842, 19172, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 30917, 3, 11382, 11392, 19637, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 30938, 3, 11392, 11402, 19652, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 30959, 3, 11402, 11412, 19667, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 30980, 3, 11412, 11422, 19682, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 31001, 3, 11422, 11432, 19697, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 31022, 3, 11432, 11442, 19712, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 31043, 3, 11442, 11452, 19727, ncols,
                                                 alpha, beta, p);

            compute_prim_sh_electron_repulsion_0(buffer, 31064, 3, 11452, 11462, 19742, ncols,
                                                 alpha, beta, p);

            compute_prim_ph_electron_repulsion_0(buffer, 31085, 0, 3, 19622, 30917, 19802, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 31148, 0, 3, 19637, 30938, 19847, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 31211, 0, 3, 19652, 30959, 19892, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 31274, 0, 3, 19667, 30980, 19937, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 31337, 0, 3, 19682, 31001, 19982, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 31400, 0, 3, 19697, 31022, 20027, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 31463, 0, 3, 19712, 31043, 20072, ncols,
                                                 p);

            compute_prim_ph_electron_repulsion_0(buffer, 31526, 0, 3, 19727, 31064, 20117, ncols,
                                                 p);

            compute_prim_dh_electron_repulsion_0(buffer, 31589, 0, 3, 19802, 31148, 11782, 11842,
                                                 20342, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 31715, 0, 3, 19847, 31211, 11842, 11902,
                                                 20432, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 31841, 0, 3, 19892, 31274, 11902, 11962,
                                                 20522, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 31967, 0, 3, 19937, 31337, 11962, 12022,
                                                 20612, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 32093, 0, 3, 19982, 31400, 12022, 12082,
                                                 20702, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 32219, 0, 3, 20027, 31463, 12082, 12142,
                                                 20792, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_0(buffer, 32345, 0, 3, 20072, 31526, 12142, 12202,
                                                 20882, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 32471, 0, 3, 20342, 31715, 12322, 12422,
                                                 21122, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 32681, 0, 3, 20432, 31841, 12422, 12522,
                                                 21272, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 32891, 0, 3, 20522, 31967, 12522, 12622,
                                                 21422, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 33101, 0, 3, 20612, 32093, 12622, 12722,
                                                 21572, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 33311, 0, 3, 20702, 32219, 12722, 12822,
                                                 21722, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_0(buffer, 33521, 0, 3, 20792, 32345, 12822, 12922,
                                                 21872, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 33731, 0, 3, 31589, 31715, 21122, 32681,
                                                 13122, 13272, 22472, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 34046, 0, 3, 31715, 31841, 21272, 32891,
                                                 13272, 13422, 22697, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 34361, 0, 3, 31841, 31967, 21422, 33101,
                                                 13422, 13572, 22922, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 34676, 0, 3, 31967, 32093, 21572, 33311,
                                                 13572, 13722, 23147, ncols, alpha, beta, p);

            compute_prim_gh_electron_repulsion_0(buffer, 34991, 0, 3, 32093, 32219, 21722, 33521,
                                                 13722, 13872, 23372, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 35306, 0, 3, 32471, 32681, 22472, 34046,
                                                 14172, 14382, 23912, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 35747, 0, 3, 32681, 32891, 22697, 34361,
                                                 14382, 14592, 24227, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 36188, 0, 3, 32891, 33101, 22922, 34676,
                                                 14592, 14802, 24542, ncols, alpha, beta, p);

            compute_prim_hh_electron_repulsion_0(buffer, 36629, 0, 3, 33101, 33311, 23147, 34991,
                                                 14802, 15012, 24857, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 37070, 0, 3, 33731, 34046, 23912, 35747,
                                                 15432, 15712, 26012, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 37658, 0, 3, 34046, 34361, 24227, 36188,
                                                 15712, 15992, 26432, ncols, alpha, beta, p);

            compute_prim_ih_electron_repulsion_0(buffer, 38246, 0, 3, 34361, 34676, 24542, 36629,
                                                 15992, 16272, 26852, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 38834, 0, 3, 35306, 35747, 26012, 37658,
                                                 16832, 17192, 27812, ncols, alpha, beta, p);

            compute_prim_kh_electron_repulsion_0(buffer, 39590, 0, 3, 35747, 36188, 26432, 38246,
                                                 17192, 17552, 28352, ncols, alpha, beta, p);

            compute_prim_lh_electron_repulsion_0(buffer, 40346, 0, 3, 37070, 37658, 27812, 39590,
                                                 18272, 18722, 30242, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41291, 3, 19622, 19637, 30938, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41319, 3, 19637, 19652, 30959, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41347, 3, 19652, 19667, 30980, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41375, 3, 19667, 19682, 31001, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41403, 3, 19682, 19697, 31022, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41431, 3, 19697, 19712, 31043, ncols,
                                                 alpha, beta, p);

            compute_prim_si_electron_repulsion_0(buffer, 41459, 3, 19712, 19727, 31064, ncols,
                                                 alpha, beta, p);

            compute_prim_pi_electron_repulsion_0(buffer, 41487, 0, 3, 30917, 41291, 31148, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41571, 0, 3, 30938, 41319, 31211, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41655, 0, 3, 30959, 41347, 31274, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41739, 0, 3, 30980, 41375, 31337, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41823, 0, 3, 31001, 41403, 31400, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41907, 0, 3, 31022, 41431, 31463, ncols,
                                                 p);

            compute_prim_pi_electron_repulsion_0(buffer, 41991, 0, 3, 31043, 41459, 31526, ncols,
                                                 p);

            compute_prim_di_electron_repulsion_0(buffer, 42075, 0, 3, 31085, 41487, 20162, 20252,
                                                 31589, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42243, 0, 3, 31148, 41571, 20252, 20342,
                                                 31715, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42411, 0, 3, 31211, 41655, 20342, 20432,
                                                 31841, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42579, 0, 3, 31274, 41739, 20432, 20522,
                                                 31967, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42747, 0, 3, 31337, 41823, 20522, 20612,
                                                 32093, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 42915, 0, 3, 31400, 41907, 20612, 20702,
                                                 32219, ncols, alpha, beta, p);

            compute_prim_di_electron_repulsion_0(buffer, 43083, 0, 3, 31463, 41991, 20702, 20792,
                                                 32345, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 43251, 0, 3, 31715, 42411, 20972, 21122,
                                                 32681, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 43531, 0, 3, 31841, 42579, 21122, 21272,
                                                 32891, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 43811, 0, 3, 31967, 42747, 21272, 21422,
                                                 33101, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 44091, 0, 3, 32093, 42915, 21422, 21572,
                                                 33311, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 44371, 0, 3, 32219, 43083, 21572, 21722,
                                                 33521, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 44651, 0, 3, 42075, 42243, 32471, 43251,
                                                 22022, 22247, 33731, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 45071, 0, 3, 42243, 42411, 32681, 43531,
                                                 22247, 22472, 34046, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 45491, 0, 3, 42411, 42579, 32891, 43811,
                                                 22472, 22697, 34361, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 45911, 0, 3, 42579, 42747, 33101, 44091,
                                                 22697, 22922, 34676, ncols, alpha, beta, p);

            compute_prim_gi_electron_repulsion_0(buffer, 46331, 0, 3, 42747, 42915, 33311, 44371,
                                                 22922, 23147, 34991, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 46751, 0, 3, 43251, 43531, 34046, 45491,
                                                 23597, 23912, 35747, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 47339, 0, 3, 43531, 43811, 34361, 45911,
                                                 23912, 24227, 36188, ncols, alpha, beta, p);

            compute_prim_hi_electron_repulsion_0(buffer, 47927, 0, 3, 43811, 44091, 34676, 46331,
                                                 24227, 24542, 36629, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 48515, 0, 3, 44651, 45071, 35306, 46751,
                                                 25172, 25592, 37070, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 49299, 0, 3, 45071, 45491, 35747, 47339,
                                                 25592, 26012, 37658, ncols, alpha, beta, p);

            compute_prim_ii_electron_repulsion_0(buffer, 50083, 0, 3, 45491, 45911, 36188, 47927,
                                                 26012, 26432, 38246, ncols, alpha, beta, p);

            compute_prim_ki_electron_repulsion_0(buffer, 50867, 0, 3, 46751, 47339, 37658, 50083,
                                                 27272, 27812, 39590, ncols, alpha, beta, p);

            compute_prim_li_electron_repulsion_0(buffer, 51875, 0, 3, 48515, 49299, 38834, 50867,
                                                 28892, 29567, 40346, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 53135, 51875, 1260, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 54395, 53135, 45, nmax);

    simdtrf::transform_l_outer(values, nvalues, buffer, 54395, 13, nmax);
}

}  // namespace simdt2ceri
