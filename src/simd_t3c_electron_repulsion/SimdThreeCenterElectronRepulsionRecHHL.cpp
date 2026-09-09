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


#include "SimdThreeCenterElectronRepulsionRecHHL.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHG.hpp"
#include "SimdTransferHH.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hhl_three_center_electron_repulsion(double               *values,
                                            const size_t          npairs,
                                            const size_t          natoms,
                                            const CBasisFunction &a_function,
                                            const CBasisFunction &b_function,
                                            const CBasisFunction &c_function,
                                            const CSimdMatrix    &coordinates,
                                            const CSimdMatrix    &c_coordinates,
                                            CSimdMatrix          &buffer,
                                            const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hhl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (npairs == 0 || natoms == 0) return;

    const auto &a_exps = a_function.exponents();

    const auto &b_exps = b_function.exponents();

    const auto &c_exps = c_function.exponents();

    const auto &a_norms = a_function.normalization_factors();

    const auto &b_norms = b_function.normalization_factors();

    const auto &c_norms = c_function.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprim_c = c_exps.size();

    const auto nprims = nprim_a * nprim_b * nprim_c;

    // NOTE: the bound neglects the position of the atom on the ket side, so the
    // columns that survive are the same for every one of them and are counted
    // once here rather than inside the loop over them.

    const auto dimensions = simdfunc::make_column_dimensions(
        a_function, b_function, c_function, npairs, coordinates,
        screenfunc::three_center_electron_repulsion_primitive_bound,
        threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 316158, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2057 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 316158, 239532, 14440, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 18,
                                                             ncols, fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 7, 8,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 8, 9,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 9, 10,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 10, 11,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 11, 12,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 12, 13,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 13, 14,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 14, 15,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 15, 16,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 16, 17,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 17, 18,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 18, 19,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 19, 20,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 20, 21,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 164, 0, 3, 21, 22,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 170, 0, 3, 22, 23,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 176, 0, 3, 23, 24,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 26, 29,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 29, 32,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 32, 35,
                                                                       92, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 35, 38,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 38, 41,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 41, 44,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 44, 47,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 47, 50,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 50, 53,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 53, 56,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 282, 0, 3, 56, 59,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 292, 0, 3, 59, 62,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 62, 65,
                                                                       152, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 312, 0, 3, 65, 68,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 322, 0, 3, 68, 71,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 332, 0, 3, 71, 74,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 342, 0, 3, 80, 86,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 357, 0, 3, 86, 92,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 372, 0, 3, 92, 98,
                                                                       202, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 387, 0, 3, 98,
                                                                       104, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 402, 0, 3, 104,
                                                                       110, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 417, 0, 3, 110,
                                                                       116, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 432, 0, 3, 116,
                                                                       122, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 447, 0, 3, 122,
                                                                       128, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 462, 0, 3, 128,
                                                                       134, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 477, 0, 3, 134,
                                                                       140, 272, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 492, 0, 3, 140,
                                                                       146, 282, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 507, 0, 3, 146,
                                                                       152, 292, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 522, 0, 3, 152,
                                                                       158, 302, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 537, 0, 3, 158,
                                                                       164, 312, 322, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 552, 0, 3, 164,
                                                                       170, 322, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 567, 0, 3, 182,
                                                                       192, 342, 357, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 192,
                                                                       202, 357, 372, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 609, 0, 3, 202,
                                                                       212, 372, 387, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 630, 0, 3, 212,
                                                                       222, 387, 402, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 651, 0, 3, 222,
                                                                       232, 402, 417, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 672, 0, 3, 232,
                                                                       242, 417, 432, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 242,
                                                                       252, 432, 447, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 714, 0, 3, 252,
                                                                       262, 447, 462, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 735, 0, 3, 262,
                                                                       272, 462, 477, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 756, 0, 3, 272,
                                                                       282, 477, 492, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 777, 0, 3, 282,
                                                                       292, 492, 507, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 798, 0, 3, 292,
                                                                       302, 507, 522, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 819, 0, 3, 302,
                                                                       312, 522, 537, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 840, 0, 3, 312,
                                                                       322, 537, 552, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 861, 0, 3, 342,
                                                                       357, 567, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 889, 0, 3, 357,
                                                                       372, 588, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 917, 0, 3, 372,
                                                                       387, 609, 630, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 945, 0, 3, 387,
                                                                       402, 630, 651, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 973, 0, 3, 402,
                                                                       417, 651, 672, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 417,
                                                                       432, 672, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 432,
                                                                       447, 693, 714, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 447,
                                                                       462, 714, 735, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1085, 0, 3, 462,
                                                                       477, 735, 756, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1113, 0, 3, 477,
                                                                       492, 756, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 492,
                                                                       507, 777, 798, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1169, 0, 3, 507,
                                                                       522, 798, 819, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1197, 0, 3, 522,
                                                                       537, 819, 840, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1225, 0, 3, 567,
                                                                       588, 861, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1261, 0, 3, 588,
                                                                       609, 889, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1297, 0, 3, 609,
                                                                       630, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1333, 0, 3, 630,
                                                                       651, 945, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1369, 0, 3, 651,
                                                                       672, 973, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1405, 0, 3, 672,
                                                                       693, 1001, 1029, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1441, 0, 3, 693,
                                                                       714, 1029, 1057, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1477, 0, 3, 714,
                                                                       735, 1057, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1513, 0, 3, 735,
                                                                       756, 1085, 1113, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1549, 0, 3, 756,
                                                                       777, 1113, 1141, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1585, 0, 3, 777,
                                                                       798, 1141, 1169, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1621, 0, 3, 798,
                                                                       819, 1169, 1197, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1657, 0, 3, 861,
                                                                       889, 1225, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1702, 0, 3, 889,
                                                                       917, 1261, 1297, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1747, 0, 3, 917,
                                                                       945, 1297, 1333, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1792, 0, 3, 945,
                                                                       973, 1333, 1369, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1837, 0, 3, 973,
                                                                       1001, 1369, 1405, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1882, 0, 3, 1001,
                                                                       1029, 1405, 1441, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1927, 0, 3, 1029,
                                                                       1057, 1441, 1477, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1972, 0, 3, 1057,
                                                                       1085, 1477, 1513, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2017, 0, 3, 1085,
                                                                       1113, 1513, 1549, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2062, 0, 3, 1113,
                                                                       1141, 1549, 1585, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2107, 0, 3, 1141,
                                                                       1169, 1585, 1621, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2152, 0, 3, 1225,
                                                                       1261, 1657, 1702, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2207, 0, 3, 1261,
                                                                       1297, 1702, 1747, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2262, 0, 3, 1297,
                                                                       1333, 1747, 1792, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2317, 0, 3, 1333,
                                                                       1369, 1792, 1837, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1369,
                                                                       1405, 1837, 1882, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2427, 0, 3, 1405,
                                                                       1441, 1882, 1927, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2482, 0, 3, 1441,
                                                                       1477, 1927, 1972, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2537, 0, 3, 1477,
                                                                       1513, 1972, 2017, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2592, 0, 3, 1513,
                                                                       1549, 2017, 2062, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2647, 0, 3, 1549,
                                                                       1585, 2062, 2107, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2702, 0, 3, 1657,
                                                                       1702, 2152, 2207, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2768, 0, 3, 1702,
                                                                       1747, 2207, 2262, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2834, 0, 3, 1747,
                                                                       1792, 2262, 2317, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2900, 0, 3, 1792,
                                                                       1837, 2317, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2966, 0, 3, 1837,
                                                                       1882, 2372, 2427, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3032, 0, 3, 1882,
                                                                       1927, 2427, 2482, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3098, 0, 3, 1927,
                                                                       1972, 2482, 2537, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3164, 0, 3, 1972,
                                                                       2017, 2537, 2592, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3230, 0, 3, 2017,
                                                                       2062, 2592, 2647, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3296, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3299, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3302, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3305, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3308, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3311, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3314, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3317, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3320, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3323, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3326, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3329, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3332, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3335, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3338, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3341, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3344, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3347, 3, 9, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3356, 3, 10, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3365, 3, 11, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3374, 3, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3383, 3, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3392, 3, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3401, 3, 15, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3410, 3, 16, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3419, 3, 17, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3428, 3, 18, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3437, 3, 19, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3446, 3, 20, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3455, 3, 21, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3464, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3473, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3482, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3491, 3, 32, 92,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3509, 3, 35, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3527, 3, 38, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3545, 3, 41, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3563, 3, 44, 116,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3581, 3, 47, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3599, 3, 50, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3617, 3, 53, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3635, 3, 56, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3653, 3, 59, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3671, 3, 62, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3689, 3, 65, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3707, 3, 68, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3725, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3743, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3761, 3, 92, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3791, 3, 98, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3821, 3, 104, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3851, 3, 110, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3881, 3, 116, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3911, 3, 122, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3941, 3, 128, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3971, 3, 134, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4001, 3, 140, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4031, 3, 146, 292,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4061, 3, 152, 302,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4091, 3, 158, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4121, 3, 164, 322,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4151, 3, 170, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4181, 3, 202, 372,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4226, 3, 212, 387,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4271, 3, 222, 402,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4316, 3, 232, 417,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4361, 3, 242, 432,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4406, 3, 252, 447,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4451, 3, 262, 462,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4496, 3, 272, 477,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4541, 3, 282, 492,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4586, 3, 292, 507,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4631, 3, 302, 522,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4676, 3, 312, 537,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4721, 3, 322, 552,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4766, 3, 372, 609,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4829, 3, 387, 630,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4892, 3, 402, 651,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4955, 3, 417, 672,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5018, 3, 432, 693,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5081, 3, 447, 714,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5144, 3, 462, 735,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5207, 3, 477, 756,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5270, 3, 492, 777,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5333, 3, 507, 798,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5396, 3, 522, 819,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5459, 3, 537, 840,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5522, 3, 609, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5606, 3, 630, 945,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5690, 3, 651, 973,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5774, 3, 672,
                                                                       1001, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5858, 3, 693,
                                                                       1029, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5942, 3, 714,
                                                                       1057, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6026, 3, 735,
                                                                       1085, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6110, 3, 756,
                                                                       1113, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6194, 3, 777,
                                                                       1141, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6278, 3, 798,
                                                                       1169, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6362, 3, 819,
                                                                       1197, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6446, 3, 917,
                                                                       1297, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6554, 3, 945,
                                                                       1333, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6662, 3, 973,
                                                                       1369, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6770, 3, 1001,
                                                                       1405, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6878, 3, 1029,
                                                                       1441, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6986, 3, 1057,
                                                                       1477, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7094, 3, 1085,
                                                                       1513, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7202, 3, 1113,
                                                                       1549, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7310, 3, 1141,
                                                                       1585, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7418, 3, 1169,
                                                                       1621, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7526, 3, 1297,
                                                                       1747, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7661, 3, 1333,
                                                                       1792, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7796, 3, 1369,
                                                                       1837, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7931, 3, 1405,
                                                                       1882, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8066, 3, 1441,
                                                                       1927, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8201, 3, 1477,
                                                                       1972, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8336, 3, 1513,
                                                                       2017, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8471, 3, 1549,
                                                                       2062, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8606, 3, 1585,
                                                                       2107, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8741, 3, 1747,
                                                                       2262, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8906, 3, 1792,
                                                                       2317, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9071, 3, 1837,
                                                                       2372, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9236, 3, 1882,
                                                                       2427, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9401, 3, 1927,
                                                                       2482, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9566, 3, 1972,
                                                                       2537, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9731, 3, 2017,
                                                                       2592, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 9896, 3, 2062,
                                                                       2647, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10061, 3, 2262,
                                                                       2834, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10259, 3, 2317,
                                                                       2900, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10457, 3, 2372,
                                                                       2966, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10655, 3, 2427,
                                                                       3032, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 10853, 3, 2482,
                                                                       3098, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11051, 3, 2537,
                                                                       3164, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11249, 3, 2592,
                                                                       3230, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11447, 3, 7, 8,
                                                                       3296, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11453, 3, 8, 9,
                                                                       3299, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11459, 3, 9, 10,
                                                                       3302, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11465, 3, 10, 11,
                                                                       3305, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11471, 3, 11, 12,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11477, 3, 12, 13,
                                                                       3311, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11483, 3, 13, 14,
                                                                       3314, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11489, 3, 14, 15,
                                                                       3317, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11495, 3, 15, 16,
                                                                       3320, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11501, 3, 16, 17,
                                                                       3323, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11507, 3, 17, 18,
                                                                       3326, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11513, 3, 18, 19,
                                                                       3329, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11519, 3, 19, 20,
                                                                       3332, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11525, 3, 20, 21,
                                                                       3335, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11531, 3, 21, 22,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11537, 3, 22, 23,
                                                                       3341, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11543, 3, 23, 24,
                                                                       3344, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11549, 0, 3,
                                                                       11447, 3296, 11453, 3347,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11567, 0, 3,
                                                                       11453, 3299, 11459, 3356,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11585, 0, 3,
                                                                       11459, 3302, 11465, 3365,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11603, 0, 3,
                                                                       11465, 3305, 11471, 3374,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11621, 0, 3,
                                                                       11471, 3308, 11477, 3383,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11639, 0, 3,
                                                                       11477, 3311, 11483, 3392,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11657, 0, 3,
                                                                       11483, 3314, 11489, 3401,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11675, 0, 3,
                                                                       11489, 3317, 11495, 3410,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11693, 0, 3,
                                                                       11495, 3320, 11501, 3419,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11711, 0, 3,
                                                                       11501, 3323, 11507, 3428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11729, 0, 3,
                                                                       11507, 3326, 11513, 3437,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11747, 0, 3,
                                                                       11513, 3329, 11519, 3446,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11765, 0, 3,
                                                                       11519, 3332, 11525, 3455,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11783, 0, 3,
                                                                       11525, 3335, 11531, 3464,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11801, 0, 3,
                                                                       11531, 3338, 11537, 3473,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11819, 0, 3,
                                                                       11537, 3341, 11543, 3482,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11837, 0, 3,
                                                                       11549, 3347, 11567, 80,
                                                                       86, 3491, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11873, 0, 3,
                                                                       11567, 3356, 11585, 86,
                                                                       92, 3509, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11909, 0, 3,
                                                                       11585, 3365, 11603, 92,
                                                                       98, 3527, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11945, 0, 3,
                                                                       11603, 3374, 11621, 98,
                                                                       104, 3545, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11981, 0, 3,
                                                                       11621, 3383, 11639, 104,
                                                                       110, 3563, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12017, 0, 3,
                                                                       11639, 3392, 11657, 110,
                                                                       116, 3581, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12053, 0, 3,
                                                                       11657, 3401, 11675, 116,
                                                                       122, 3599, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12089, 0, 3,
                                                                       11675, 3410, 11693, 122,
                                                                       128, 3617, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12125, 0, 3,
                                                                       11693, 3419, 11711, 128,
                                                                       134, 3635, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12161, 0, 3,
                                                                       11711, 3428, 11729, 134,
                                                                       140, 3653, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12197, 0, 3,
                                                                       11729, 3437, 11747, 140,
                                                                       146, 3671, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12233, 0, 3,
                                                                       11747, 3446, 11765, 146,
                                                                       152, 3689, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12269, 0, 3,
                                                                       11765, 3455, 11783, 152,
                                                                       158, 3707, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12305, 0, 3,
                                                                       11783, 3464, 11801, 158,
                                                                       164, 3725, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 12341, 0, 3,
                                                                       11801, 3473, 11819, 164,
                                                                       170, 3743, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12377, 0, 3,
                                                                       11837, 3491, 11873, 182,
                                                                       192, 3761, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12437, 0, 3,
                                                                       11873, 3509, 11909, 192,
                                                                       202, 3791, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12497, 0, 3,
                                                                       11909, 3527, 11945, 202,
                                                                       212, 3821, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12557, 0, 3,
                                                                       11945, 3545, 11981, 212,
                                                                       222, 3851, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12617, 0, 3,
                                                                       11981, 3563, 12017, 222,
                                                                       232, 3881, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12677, 0, 3,
                                                                       12017, 3581, 12053, 232,
                                                                       242, 3911, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12737, 0, 3,
                                                                       12053, 3599, 12089, 242,
                                                                       252, 3941, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12797, 0, 3,
                                                                       12089, 3617, 12125, 252,
                                                                       262, 3971, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12857, 0, 3,
                                                                       12125, 3635, 12161, 262,
                                                                       272, 4001, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12917, 0, 3,
                                                                       12161, 3653, 12197, 272,
                                                                       282, 4031, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12977, 0, 3,
                                                                       12197, 3671, 12233, 282,
                                                                       292, 4061, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13037, 0, 3,
                                                                       12233, 3689, 12269, 292,
                                                                       302, 4091, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13097, 0, 3,
                                                                       12269, 3707, 12305, 302,
                                                                       312, 4121, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 13157, 0, 3,
                                                                       12305, 3725, 12341, 312,
                                                                       322, 4151, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13217, 0, 3,
                                                                       12377, 3761, 12437, 342,
                                                                       357, 4181, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13307, 0, 3,
                                                                       12437, 3791, 12497, 357,
                                                                       372, 4226, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13397, 0, 3,
                                                                       12497, 3821, 12557, 372,
                                                                       387, 4271, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13487, 0, 3,
                                                                       12557, 3851, 12617, 387,
                                                                       402, 4316, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13577, 0, 3,
                                                                       12617, 3881, 12677, 402,
                                                                       417, 4361, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13667, 0, 3,
                                                                       12677, 3911, 12737, 417,
                                                                       432, 4406, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13757, 0, 3,
                                                                       12737, 3941, 12797, 432,
                                                                       447, 4451, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13847, 0, 3,
                                                                       12797, 3971, 12857, 447,
                                                                       462, 4496, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13937, 0, 3,
                                                                       12857, 4001, 12917, 462,
                                                                       477, 4541, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14027, 0, 3,
                                                                       12917, 4031, 12977, 477,
                                                                       492, 4586, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14117, 0, 3,
                                                                       12977, 4061, 13037, 492,
                                                                       507, 4631, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14207, 0, 3,
                                                                       13037, 4091, 13097, 507,
                                                                       522, 4676, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14297, 0, 3,
                                                                       13097, 4121, 13157, 522,
                                                                       537, 4721, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14387, 0, 3,
                                                                       13217, 4181, 13307, 567,
                                                                       588, 4766, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14513, 0, 3,
                                                                       13307, 4226, 13397, 588,
                                                                       609, 4829, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14639, 0, 3,
                                                                       13397, 4271, 13487, 609,
                                                                       630, 4892, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14765, 0, 3,
                                                                       13487, 4316, 13577, 630,
                                                                       651, 4955, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14891, 0, 3,
                                                                       13577, 4361, 13667, 651,
                                                                       672, 5018, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15017, 0, 3,
                                                                       13667, 4406, 13757, 672,
                                                                       693, 5081, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15143, 0, 3,
                                                                       13757, 4451, 13847, 693,
                                                                       714, 5144, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15269, 0, 3,
                                                                       13847, 4496, 13937, 714,
                                                                       735, 5207, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15395, 0, 3,
                                                                       13937, 4541, 14027, 735,
                                                                       756, 5270, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15521, 0, 3,
                                                                       14027, 4586, 14117, 756,
                                                                       777, 5333, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15647, 0, 3,
                                                                       14117, 4631, 14207, 777,
                                                                       798, 5396, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15773, 0, 3,
                                                                       14207, 4676, 14297, 798,
                                                                       819, 5459, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15899, 0, 3,
                                                                       14387, 4766, 14513, 861,
                                                                       889, 5522, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16067, 0, 3,
                                                                       14513, 4829, 14639, 889,
                                                                       917, 5606, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16235, 0, 3,
                                                                       14639, 4892, 14765, 917,
                                                                       945, 5690, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16403, 0, 3,
                                                                       14765, 4955, 14891, 945,
                                                                       973, 5774, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16571, 0, 3,
                                                                       14891, 5018, 15017, 973,
                                                                       1001, 5858, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16739, 0, 3,
                                                                       15017, 5081, 15143, 1001,
                                                                       1029, 5942, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16907, 0, 3,
                                                                       15143, 5144, 15269, 1029,
                                                                       1057, 6026, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17075, 0, 3,
                                                                       15269, 5207, 15395, 1057,
                                                                       1085, 6110, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17243, 0, 3,
                                                                       15395, 5270, 15521, 1085,
                                                                       1113, 6194, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17411, 0, 3,
                                                                       15521, 5333, 15647, 1113,
                                                                       1141, 6278, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17579, 0, 3,
                                                                       15647, 5396, 15773, 1141,
                                                                       1169, 6362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17747, 0, 3,
                                                                       15899, 5522, 16067, 1225,
                                                                       1261, 6446, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17963, 0, 3,
                                                                       16067, 5606, 16235, 1261,
                                                                       1297, 6554, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18179, 0, 3,
                                                                       16235, 5690, 16403, 1297,
                                                                       1333, 6662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18395, 0, 3,
                                                                       16403, 5774, 16571, 1333,
                                                                       1369, 6770, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18611, 0, 3,
                                                                       16571, 5858, 16739, 1369,
                                                                       1405, 6878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18827, 0, 3,
                                                                       16739, 5942, 16907, 1405,
                                                                       1441, 6986, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19043, 0, 3,
                                                                       16907, 6026, 17075, 1441,
                                                                       1477, 7094, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19259, 0, 3,
                                                                       17075, 6110, 17243, 1477,
                                                                       1513, 7202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19475, 0, 3,
                                                                       17243, 6194, 17411, 1513,
                                                                       1549, 7310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19691, 0, 3,
                                                                       17411, 6278, 17579, 1549,
                                                                       1585, 7418, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 19907, 0, 3,
                                                                       17747, 6446, 17963, 1657,
                                                                       1702, 7526, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20177, 0, 3,
                                                                       17963, 6554, 18179, 1702,
                                                                       1747, 7661, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20447, 0, 3,
                                                                       18179, 6662, 18395, 1747,
                                                                       1792, 7796, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20717, 0, 3,
                                                                       18395, 6770, 18611, 1792,
                                                                       1837, 7931, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 20987, 0, 3,
                                                                       18611, 6878, 18827, 1837,
                                                                       1882, 8066, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21257, 0, 3,
                                                                       18827, 6986, 19043, 1882,
                                                                       1927, 8201, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21527, 0, 3,
                                                                       19043, 7094, 19259, 1927,
                                                                       1972, 8336, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21797, 0, 3,
                                                                       19259, 7202, 19475, 1972,
                                                                       2017, 8471, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22067, 0, 3,
                                                                       19475, 7310, 19691, 2017,
                                                                       2062, 8606, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22337, 0, 3,
                                                                       19907, 7526, 20177, 2152,
                                                                       2207, 8741, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22667, 0, 3,
                                                                       20177, 7661, 20447, 2207,
                                                                       2262, 8906, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 22997, 0, 3,
                                                                       20447, 7796, 20717, 2262,
                                                                       2317, 9071, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23327, 0, 3,
                                                                       20717, 7931, 20987, 2317,
                                                                       2372, 9236, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23657, 0, 3,
                                                                       20987, 8066, 21257, 2372,
                                                                       2427, 9401, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 23987, 0, 3,
                                                                       21257, 8201, 21527, 2427,
                                                                       2482, 9566, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24317, 0, 3,
                                                                       21527, 8336, 21797, 2482,
                                                                       2537, 9731, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 24647, 0, 3,
                                                                       21797, 8471, 22067, 2537,
                                                                       2592, 9896, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 24977, 0, 3,
                                                                       22337, 8741, 22667, 2702,
                                                                       2768, 10061, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25373, 0, 3,
                                                                       22667, 8906, 22997, 2768,
                                                                       2834, 10259, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 25769, 0, 3,
                                                                       22997, 9071, 23327, 2834,
                                                                       2900, 10457, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26165, 0, 3,
                                                                       23327, 9236, 23657, 2900,
                                                                       2966, 10655, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26561, 0, 3,
                                                                       23657, 9401, 23987, 2966,
                                                                       3032, 10853, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 26957, 0, 3,
                                                                       23987, 9566, 24317, 3032,
                                                                       3098, 11051, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 27353, 0, 3,
                                                                       24317, 9731, 24647, 3098,
                                                                       3164, 11249, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27749, 3, 3296,
                                                                       3299, 11459, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27759, 3, 3299,
                                                                       3302, 11465, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27769, 3, 3302,
                                                                       3305, 11471, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27779, 3, 3305,
                                                                       3308, 11477, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27789, 3, 3308,
                                                                       3311, 11483, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27799, 3, 3311,
                                                                       3314, 11489, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27809, 3, 3314,
                                                                       3317, 11495, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27819, 3, 3317,
                                                                       3320, 11501, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27829, 3, 3320,
                                                                       3323, 11507, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27839, 3, 3323,
                                                                       3326, 11513, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27849, 3, 3326,
                                                                       3329, 11519, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27859, 3, 3329,
                                                                       3332, 11525, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27869, 3, 3332,
                                                                       3335, 11531, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27879, 3, 3335,
                                                                       3338, 11537, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 27889, 3, 3338,
                                                                       3341, 11543, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27899, 0, 3,
                                                                       27749, 11459, 27759, 3347,
                                                                       3356, 11585, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27929, 0, 3,
                                                                       27759, 11465, 27769, 3356,
                                                                       3365, 11603, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27959, 0, 3,
                                                                       27769, 11471, 27779, 3365,
                                                                       3374, 11621, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 27989, 0, 3,
                                                                       27779, 11477, 27789, 3374,
                                                                       3383, 11639, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28019, 0, 3,
                                                                       27789, 11483, 27799, 3383,
                                                                       3392, 11657, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28049, 0, 3,
                                                                       27799, 11489, 27809, 3392,
                                                                       3401, 11675, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28079, 0, 3,
                                                                       27809, 11495, 27819, 3401,
                                                                       3410, 11693, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28109, 0, 3,
                                                                       27819, 11501, 27829, 3410,
                                                                       3419, 11711, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28139, 0, 3,
                                                                       27829, 11507, 27839, 3419,
                                                                       3428, 11729, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28169, 0, 3,
                                                                       27839, 11513, 27849, 3428,
                                                                       3437, 11747, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28199, 0, 3,
                                                                       27849, 11519, 27859, 3437,
                                                                       3446, 11765, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28229, 0, 3,
                                                                       27859, 11525, 27869, 3446,
                                                                       3455, 11783, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28259, 0, 3,
                                                                       27869, 11531, 27879, 3455,
                                                                       3464, 11801, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 28289, 0, 3,
                                                                       27879, 11537, 27889, 3464,
                                                                       3473, 11819, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28319, 0, 3,
                                                                       27899, 11585, 27929, 3491,
                                                                       3509, 11909, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28379, 0, 3,
                                                                       27929, 11603, 27959, 3509,
                                                                       3527, 11945, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28439, 0, 3,
                                                                       27959, 11621, 27989, 3527,
                                                                       3545, 11981, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28499, 0, 3,
                                                                       27989, 11639, 28019, 3545,
                                                                       3563, 12017, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28559, 0, 3,
                                                                       28019, 11657, 28049, 3563,
                                                                       3581, 12053, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28619, 0, 3,
                                                                       28049, 11675, 28079, 3581,
                                                                       3599, 12089, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28679, 0, 3,
                                                                       28079, 11693, 28109, 3599,
                                                                       3617, 12125, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28739, 0, 3,
                                                                       28109, 11711, 28139, 3617,
                                                                       3635, 12161, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28799, 0, 3,
                                                                       28139, 11729, 28169, 3635,
                                                                       3653, 12197, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28859, 0, 3,
                                                                       28169, 11747, 28199, 3653,
                                                                       3671, 12233, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28919, 0, 3,
                                                                       28199, 11765, 28229, 3671,
                                                                       3689, 12269, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 28979, 0, 3,
                                                                       28229, 11783, 28259, 3689,
                                                                       3707, 12305, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 29039, 0, 3,
                                                                       28259, 11801, 28289, 3707,
                                                                       3725, 12341, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29099, 0, 3,
                                                                       28319, 11909, 28379, 3761,
                                                                       3791, 12497, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29199, 0, 3,
                                                                       28379, 11945, 28439, 3791,
                                                                       3821, 12557, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29299, 0, 3,
                                                                       28439, 11981, 28499, 3821,
                                                                       3851, 12617, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29399, 0, 3,
                                                                       28499, 12017, 28559, 3851,
                                                                       3881, 12677, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29499, 0, 3,
                                                                       28559, 12053, 28619, 3881,
                                                                       3911, 12737, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29599, 0, 3,
                                                                       28619, 12089, 28679, 3911,
                                                                       3941, 12797, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29699, 0, 3,
                                                                       28679, 12125, 28739, 3941,
                                                                       3971, 12857, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29799, 0, 3,
                                                                       28739, 12161, 28799, 3971,
                                                                       4001, 12917, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29899, 0, 3,
                                                                       28799, 12197, 28859, 4001,
                                                                       4031, 12977, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 29999, 0, 3,
                                                                       28859, 12233, 28919, 4031,
                                                                       4061, 13037, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30099, 0, 3,
                                                                       28919, 12269, 28979, 4061,
                                                                       4091, 13097, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 30199, 0, 3,
                                                                       28979, 12305, 29039, 4091,
                                                                       4121, 13157, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30299, 0, 3,
                                                                       29099, 12497, 29199, 4181,
                                                                       4226, 13397, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30449, 0, 3,
                                                                       29199, 12557, 29299, 4226,
                                                                       4271, 13487, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30599, 0, 3,
                                                                       29299, 12617, 29399, 4271,
                                                                       4316, 13577, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30749, 0, 3,
                                                                       29399, 12677, 29499, 4316,
                                                                       4361, 13667, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 30899, 0, 3,
                                                                       29499, 12737, 29599, 4361,
                                                                       4406, 13757, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31049, 0, 3,
                                                                       29599, 12797, 29699, 4406,
                                                                       4451, 13847, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31199, 0, 3,
                                                                       29699, 12857, 29799, 4451,
                                                                       4496, 13937, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31349, 0, 3,
                                                                       29799, 12917, 29899, 4496,
                                                                       4541, 14027, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31499, 0, 3,
                                                                       29899, 12977, 29999, 4541,
                                                                       4586, 14117, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31649, 0, 3,
                                                                       29999, 13037, 30099, 4586,
                                                                       4631, 14207, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 31799, 0, 3,
                                                                       30099, 13097, 30199, 4631,
                                                                       4676, 14297, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31949, 0, 3,
                                                                       30299, 13397, 30449, 4766,
                                                                       4829, 14639, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32159, 0, 3,
                                                                       30449, 13487, 30599, 4829,
                                                                       4892, 14765, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32369, 0, 3,
                                                                       30599, 13577, 30749, 4892,
                                                                       4955, 14891, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32579, 0, 3,
                                                                       30749, 13667, 30899, 4955,
                                                                       5018, 15017, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32789, 0, 3,
                                                                       30899, 13757, 31049, 5018,
                                                                       5081, 15143, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 32999, 0, 3,
                                                                       31049, 13847, 31199, 5081,
                                                                       5144, 15269, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33209, 0, 3,
                                                                       31199, 13937, 31349, 5144,
                                                                       5207, 15395, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33419, 0, 3,
                                                                       31349, 14027, 31499, 5207,
                                                                       5270, 15521, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33629, 0, 3,
                                                                       31499, 14117, 31649, 5270,
                                                                       5333, 15647, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 33839, 0, 3,
                                                                       31649, 14207, 31799, 5333,
                                                                       5396, 15773, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34049, 0, 3,
                                                                       31949, 14639, 32159, 5522,
                                                                       5606, 16235, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34329, 0, 3,
                                                                       32159, 14765, 32369, 5606,
                                                                       5690, 16403, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34609, 0, 3,
                                                                       32369, 14891, 32579, 5690,
                                                                       5774, 16571, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34889, 0, 3,
                                                                       32579, 15017, 32789, 5774,
                                                                       5858, 16739, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35169, 0, 3,
                                                                       32789, 15143, 32999, 5858,
                                                                       5942, 16907, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35449, 0, 3,
                                                                       32999, 15269, 33209, 5942,
                                                                       6026, 17075, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 35729, 0, 3,
                                                                       33209, 15395, 33419, 6026,
                                                                       6110, 17243, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36009, 0, 3,
                                                                       33419, 15521, 33629, 6110,
                                                                       6194, 17411, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 36289, 0, 3,
                                                                       33629, 15647, 33839, 6194,
                                                                       6278, 17579, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36569, 0, 3,
                                                                       34049, 16235, 34329, 6446,
                                                                       6554, 18179, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36929, 0, 3,
                                                                       34329, 16403, 34609, 6554,
                                                                       6662, 18395, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37289, 0, 3,
                                                                       34609, 16571, 34889, 6662,
                                                                       6770, 18611, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37649, 0, 3,
                                                                       34889, 16739, 35169, 6770,
                                                                       6878, 18827, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38009, 0, 3,
                                                                       35169, 16907, 35449, 6878,
                                                                       6986, 19043, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38369, 0, 3,
                                                                       35449, 17075, 35729, 6986,
                                                                       7094, 19259, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 38729, 0, 3,
                                                                       35729, 17243, 36009, 7094,
                                                                       7202, 19475, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 39089, 0, 3,
                                                                       36009, 17411, 36289, 7202,
                                                                       7310, 19691, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39449, 0, 3,
                                                                       36569, 18179, 36929, 7526,
                                                                       7661, 20447, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39899, 0, 3,
                                                                       36929, 18395, 37289, 7661,
                                                                       7796, 20717, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 40349, 0, 3,
                                                                       37289, 18611, 37649, 7796,
                                                                       7931, 20987, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 40799, 0, 3,
                                                                       37649, 18827, 38009, 7931,
                                                                       8066, 21257, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 41249, 0, 3,
                                                                       38009, 19043, 38369, 8066,
                                                                       8201, 21527, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 41699, 0, 3,
                                                                       38369, 19259, 38729, 8201,
                                                                       8336, 21797, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 42149, 0, 3,
                                                                       38729, 19475, 39089, 8336,
                                                                       8471, 22067, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 42599, 0, 3,
                                                                       39449, 20447, 39899, 8741,
                                                                       8906, 22997, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 43149, 0, 3,
                                                                       39899, 20717, 40349, 8906,
                                                                       9071, 23327, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 43699, 0, 3,
                                                                       40349, 20987, 40799, 9071,
                                                                       9236, 23657, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 44249, 0, 3,
                                                                       40799, 21257, 41249, 9236,
                                                                       9401, 23987, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 44799, 0, 3,
                                                                       41249, 21527, 41699, 9401,
                                                                       9566, 24317, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 45349, 0, 3,
                                                                       41699, 21797, 42149, 9566,
                                                                       9731, 24647, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 45899, 0, 3,
                                                                       42599, 22997, 43149,
                                                                       10061, 10259, 25769,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 46559, 0, 3,
                                                                       43149, 23327, 43699,
                                                                       10259, 10457, 26165,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 47219, 0, 3,
                                                                       43699, 23657, 44249,
                                                                       10457, 10655, 26561,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 47879, 0, 3,
                                                                       44249, 23987, 44799,
                                                                       10655, 10853, 26957,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 48539, 0, 3,
                                                                       44799, 24317, 45349,
                                                                       10853, 11051, 27353,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49199, 3, 11447,
                                                                       11453, 27749, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49214, 3, 11453,
                                                                       11459, 27759, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49229, 3, 11459,
                                                                       11465, 27769, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49244, 3, 11465,
                                                                       11471, 27779, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49259, 3, 11471,
                                                                       11477, 27789, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49274, 3, 11477,
                                                                       11483, 27799, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49289, 3, 11483,
                                                                       11489, 27809, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49304, 3, 11489,
                                                                       11495, 27819, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49319, 3, 11495,
                                                                       11501, 27829, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49334, 3, 11501,
                                                                       11507, 27839, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49349, 3, 11507,
                                                                       11513, 27849, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49364, 3, 11513,
                                                                       11519, 27859, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49379, 3, 11519,
                                                                       11525, 27869, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49394, 3, 11525,
                                                                       11531, 27879, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49409, 3, 11531,
                                                                       11537, 27889, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49424, 0, 3,
                                                                       49199, 27749, 49214,
                                                                       11549, 11567, 27899,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49469, 0, 3,
                                                                       49214, 27759, 49229,
                                                                       11567, 11585, 27929,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49514, 0, 3,
                                                                       49229, 27769, 49244,
                                                                       11585, 11603, 27959,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49559, 0, 3,
                                                                       49244, 27779, 49259,
                                                                       11603, 11621, 27989,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49604, 0, 3,
                                                                       49259, 27789, 49274,
                                                                       11621, 11639, 28019,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49649, 0, 3,
                                                                       49274, 27799, 49289,
                                                                       11639, 11657, 28049,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49694, 0, 3,
                                                                       49289, 27809, 49304,
                                                                       11657, 11675, 28079,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49739, 0, 3,
                                                                       49304, 27819, 49319,
                                                                       11675, 11693, 28109,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49784, 0, 3,
                                                                       49319, 27829, 49334,
                                                                       11693, 11711, 28139,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49829, 0, 3,
                                                                       49334, 27839, 49349,
                                                                       11711, 11729, 28169,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49874, 0, 3,
                                                                       49349, 27849, 49364,
                                                                       11729, 11747, 28199,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49919, 0, 3,
                                                                       49364, 27859, 49379,
                                                                       11747, 11765, 28229,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 49964, 0, 3,
                                                                       49379, 27869, 49394,
                                                                       11765, 11783, 28259,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 50009, 0, 3,
                                                                       49394, 27879, 49409,
                                                                       11783, 11801, 28289,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50054, 0, 3,
                                                                       49424, 27899, 49469,
                                                                       11837, 11873, 28319,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50144, 0, 3,
                                                                       49469, 27929, 49514,
                                                                       11873, 11909, 28379,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50234, 0, 3,
                                                                       49514, 27959, 49559,
                                                                       11909, 11945, 28439,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50324, 0, 3,
                                                                       49559, 27989, 49604,
                                                                       11945, 11981, 28499,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50414, 0, 3,
                                                                       49604, 28019, 49649,
                                                                       11981, 12017, 28559,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50504, 0, 3,
                                                                       49649, 28049, 49694,
                                                                       12017, 12053, 28619,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50594, 0, 3,
                                                                       49694, 28079, 49739,
                                                                       12053, 12089, 28679,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50684, 0, 3,
                                                                       49739, 28109, 49784,
                                                                       12089, 12125, 28739,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50774, 0, 3,
                                                                       49784, 28139, 49829,
                                                                       12125, 12161, 28799,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50864, 0, 3,
                                                                       49829, 28169, 49874,
                                                                       12161, 12197, 28859,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 50954, 0, 3,
                                                                       49874, 28199, 49919,
                                                                       12197, 12233, 28919,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51044, 0, 3,
                                                                       49919, 28229, 49964,
                                                                       12233, 12269, 28979,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 51134, 0, 3,
                                                                       49964, 28259, 50009,
                                                                       12269, 12305, 29039,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51224, 0, 3,
                                                                       50054, 28319, 50144,
                                                                       12377, 12437, 29099,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51374, 0, 3,
                                                                       50144, 28379, 50234,
                                                                       12437, 12497, 29199,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51524, 0, 3,
                                                                       50234, 28439, 50324,
                                                                       12497, 12557, 29299,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51674, 0, 3,
                                                                       50324, 28499, 50414,
                                                                       12557, 12617, 29399,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51824, 0, 3,
                                                                       50414, 28559, 50504,
                                                                       12617, 12677, 29499,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 51974, 0, 3,
                                                                       50504, 28619, 50594,
                                                                       12677, 12737, 29599,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52124, 0, 3,
                                                                       50594, 28679, 50684,
                                                                       12737, 12797, 29699,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52274, 0, 3,
                                                                       50684, 28739, 50774,
                                                                       12797, 12857, 29799,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52424, 0, 3,
                                                                       50774, 28799, 50864,
                                                                       12857, 12917, 29899,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52574, 0, 3,
                                                                       50864, 28859, 50954,
                                                                       12917, 12977, 29999,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52724, 0, 3,
                                                                       50954, 28919, 51044,
                                                                       12977, 13037, 30099,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 52874, 0, 3,
                                                                       51044, 28979, 51134,
                                                                       13037, 13097, 30199,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53024, 0, 3,
                                                                       51224, 29099, 51374,
                                                                       13217, 13307, 30299,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53249, 0, 3,
                                                                       51374, 29199, 51524,
                                                                       13307, 13397, 30449,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53474, 0, 3,
                                                                       51524, 29299, 51674,
                                                                       13397, 13487, 30599,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53699, 0, 3,
                                                                       51674, 29399, 51824,
                                                                       13487, 13577, 30749,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 53924, 0, 3,
                                                                       51824, 29499, 51974,
                                                                       13577, 13667, 30899,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54149, 0, 3,
                                                                       51974, 29599, 52124,
                                                                       13667, 13757, 31049,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54374, 0, 3,
                                                                       52124, 29699, 52274,
                                                                       13757, 13847, 31199,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54599, 0, 3,
                                                                       52274, 29799, 52424,
                                                                       13847, 13937, 31349,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 54824, 0, 3,
                                                                       52424, 29899, 52574,
                                                                       13937, 14027, 31499,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55049, 0, 3,
                                                                       52574, 29999, 52724,
                                                                       14027, 14117, 31649,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 55274, 0, 3,
                                                                       52724, 30099, 52874,
                                                                       14117, 14207, 31799,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 55499, 0, 3,
                                                                       53024, 30299, 53249,
                                                                       14387, 14513, 31949,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 55814, 0, 3,
                                                                       53249, 30449, 53474,
                                                                       14513, 14639, 32159,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56129, 0, 3,
                                                                       53474, 30599, 53699,
                                                                       14639, 14765, 32369,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56444, 0, 3,
                                                                       53699, 30749, 53924,
                                                                       14765, 14891, 32579,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 56759, 0, 3,
                                                                       53924, 30899, 54149,
                                                                       14891, 15017, 32789,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57074, 0, 3,
                                                                       54149, 31049, 54374,
                                                                       15017, 15143, 32999,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57389, 0, 3,
                                                                       54374, 31199, 54599,
                                                                       15143, 15269, 33209,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 57704, 0, 3,
                                                                       54599, 31349, 54824,
                                                                       15269, 15395, 33419,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58019, 0, 3,
                                                                       54824, 31499, 55049,
                                                                       15395, 15521, 33629,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 58334, 0, 3,
                                                                       55049, 31649, 55274,
                                                                       15521, 15647, 33839,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 58649, 0, 3,
                                                                       55499, 31949, 55814,
                                                                       15899, 16067, 34049,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 59069, 0, 3,
                                                                       55814, 32159, 56129,
                                                                       16067, 16235, 34329,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 59489, 0, 3,
                                                                       56129, 32369, 56444,
                                                                       16235, 16403, 34609,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 59909, 0, 3,
                                                                       56444, 32579, 56759,
                                                                       16403, 16571, 34889,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 60329, 0, 3,
                                                                       56759, 32789, 57074,
                                                                       16571, 16739, 35169,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 60749, 0, 3,
                                                                       57074, 32999, 57389,
                                                                       16739, 16907, 35449,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61169, 0, 3,
                                                                       57389, 33209, 57704,
                                                                       16907, 17075, 35729,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 61589, 0, 3,
                                                                       57704, 33419, 58019,
                                                                       17075, 17243, 36009,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 62009, 0, 3,
                                                                       58019, 33629, 58334,
                                                                       17243, 17411, 36289,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 62429, 0, 3,
                                                                       58649, 34049, 59069,
                                                                       17747, 17963, 36569,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 62969, 0, 3,
                                                                       59069, 34329, 59489,
                                                                       17963, 18179, 36929,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 63509, 0, 3,
                                                                       59489, 34609, 59909,
                                                                       18179, 18395, 37289,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 64049, 0, 3,
                                                                       59909, 34889, 60329,
                                                                       18395, 18611, 37649,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 64589, 0, 3,
                                                                       60329, 35169, 60749,
                                                                       18611, 18827, 38009,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65129, 0, 3,
                                                                       60749, 35449, 61169,
                                                                       18827, 19043, 38369,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 65669, 0, 3,
                                                                       61169, 35729, 61589,
                                                                       19043, 19259, 38729,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 66209, 0, 3,
                                                                       61589, 36009, 62009,
                                                                       19259, 19475, 39089,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 66749, 0, 3,
                                                                       62429, 36569, 62969,
                                                                       19907, 20177, 39449,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 67424, 0, 3,
                                                                       62969, 36929, 63509,
                                                                       20177, 20447, 39899,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 68099, 0, 3,
                                                                       63509, 37289, 64049,
                                                                       20447, 20717, 40349,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 68774, 0, 3,
                                                                       64049, 37649, 64589,
                                                                       20717, 20987, 40799,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 69449, 0, 3,
                                                                       64589, 38009, 65129,
                                                                       20987, 21257, 41249,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 70124, 0, 3,
                                                                       65129, 38369, 65669,
                                                                       21257, 21527, 41699,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 70799, 0, 3,
                                                                       65669, 38729, 66209,
                                                                       21527, 21797, 42149,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 71474, 0, 3,
                                                                       66749, 39449, 67424,
                                                                       22337, 22667, 42599,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 72299, 0, 3,
                                                                       67424, 39899, 68099,
                                                                       22667, 22997, 43149,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 73124, 0, 3,
                                                                       68099, 40349, 68774,
                                                                       22997, 23327, 43699,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 73949, 0, 3,
                                                                       68774, 40799, 69449,
                                                                       23327, 23657, 44249,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 74774, 0, 3,
                                                                       69449, 41249, 70124,
                                                                       23657, 23987, 44799,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 75599, 0, 3,
                                                                       70124, 41699, 70799,
                                                                       23987, 24317, 45349,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 76424, 0, 3,
                                                                       71474, 42599, 72299,
                                                                       24977, 25373, 45899,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 77414, 0, 3,
                                                                       72299, 43149, 73124,
                                                                       25373, 25769, 46559,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 78404, 0, 3,
                                                                       73124, 43699, 73949,
                                                                       25769, 26165, 47219,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 79394, 0, 3,
                                                                       73949, 44249, 74774,
                                                                       26165, 26561, 47879,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 80384, 0, 3,
                                                                       74774, 44799, 75599,
                                                                       26561, 26957, 48539,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81374, 3, 27749,
                                                                       27759, 49229, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81395, 3, 27759,
                                                                       27769, 49244, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81416, 3, 27769,
                                                                       27779, 49259, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81437, 3, 27779,
                                                                       27789, 49274, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81458, 3, 27789,
                                                                       27799, 49289, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81479, 3, 27799,
                                                                       27809, 49304, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81500, 3, 27809,
                                                                       27819, 49319, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81521, 3, 27819,
                                                                       27829, 49334, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81542, 3, 27829,
                                                                       27839, 49349, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81563, 3, 27839,
                                                                       27849, 49364, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81584, 3, 27849,
                                                                       27859, 49379, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81605, 3, 27859,
                                                                       27869, 49394, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81626, 3, 27869,
                                                                       27879, 49409, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81647, 0, 3,
                                                                       81374, 49229, 81395,
                                                                       27899, 27929, 49514,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81710, 0, 3,
                                                                       81395, 49244, 81416,
                                                                       27929, 27959, 49559,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81773, 0, 3,
                                                                       81416, 49259, 81437,
                                                                       27959, 27989, 49604,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81836, 0, 3,
                                                                       81437, 49274, 81458,
                                                                       27989, 28019, 49649,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81899, 0, 3,
                                                                       81458, 49289, 81479,
                                                                       28019, 28049, 49694,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 81962, 0, 3,
                                                                       81479, 49304, 81500,
                                                                       28049, 28079, 49739,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82025, 0, 3,
                                                                       81500, 49319, 81521,
                                                                       28079, 28109, 49784,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82088, 0, 3,
                                                                       81521, 49334, 81542,
                                                                       28109, 28139, 49829,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82151, 0, 3,
                                                                       81542, 49349, 81563,
                                                                       28139, 28169, 49874,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82214, 0, 3,
                                                                       81563, 49364, 81584,
                                                                       28169, 28199, 49919,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82277, 0, 3,
                                                                       81584, 49379, 81605,
                                                                       28199, 28229, 49964,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 82340, 0, 3,
                                                                       81605, 49394, 81626,
                                                                       28229, 28259, 50009,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82403, 0, 3,
                                                                       81647, 49514, 81710,
                                                                       28319, 28379, 50234,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82529, 0, 3,
                                                                       81710, 49559, 81773,
                                                                       28379, 28439, 50324,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82655, 0, 3,
                                                                       81773, 49604, 81836,
                                                                       28439, 28499, 50414,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82781, 0, 3,
                                                                       81836, 49649, 81899,
                                                                       28499, 28559, 50504,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 82907, 0, 3,
                                                                       81899, 49694, 81962,
                                                                       28559, 28619, 50594,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83033, 0, 3,
                                                                       81962, 49739, 82025,
                                                                       28619, 28679, 50684,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83159, 0, 3,
                                                                       82025, 49784, 82088,
                                                                       28679, 28739, 50774,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83285, 0, 3,
                                                                       82088, 49829, 82151,
                                                                       28739, 28799, 50864,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83411, 0, 3,
                                                                       82151, 49874, 82214,
                                                                       28799, 28859, 50954,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83537, 0, 3,
                                                                       82214, 49919, 82277,
                                                                       28859, 28919, 51044,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 83663, 0, 3,
                                                                       82277, 49964, 82340,
                                                                       28919, 28979, 51134,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 83789, 0, 3,
                                                                       82403, 50234, 82529,
                                                                       29099, 29199, 51524,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 83999, 0, 3,
                                                                       82529, 50324, 82655,
                                                                       29199, 29299, 51674,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84209, 0, 3,
                                                                       82655, 50414, 82781,
                                                                       29299, 29399, 51824,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84419, 0, 3,
                                                                       82781, 50504, 82907,
                                                                       29399, 29499, 51974,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84629, 0, 3,
                                                                       82907, 50594, 83033,
                                                                       29499, 29599, 52124,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 84839, 0, 3,
                                                                       83033, 50684, 83159,
                                                                       29599, 29699, 52274,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85049, 0, 3,
                                                                       83159, 50774, 83285,
                                                                       29699, 29799, 52424,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85259, 0, 3,
                                                                       83285, 50864, 83411,
                                                                       29799, 29899, 52574,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85469, 0, 3,
                                                                       83411, 50954, 83537,
                                                                       29899, 29999, 52724,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 85679, 0, 3,
                                                                       83537, 51044, 83663,
                                                                       29999, 30099, 52874,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 85889, 0, 3,
                                                                       83789, 51524, 83999,
                                                                       30299, 30449, 53474,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86204, 0, 3,
                                                                       83999, 51674, 84209,
                                                                       30449, 30599, 53699,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86519, 0, 3,
                                                                       84209, 51824, 84419,
                                                                       30599, 30749, 53924,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 86834, 0, 3,
                                                                       84419, 51974, 84629,
                                                                       30749, 30899, 54149,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 87149, 0, 3,
                                                                       84629, 52124, 84839,
                                                                       30899, 31049, 54374,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 87464, 0, 3,
                                                                       84839, 52274, 85049,
                                                                       31049, 31199, 54599,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 87779, 0, 3,
                                                                       85049, 52424, 85259,
                                                                       31199, 31349, 54824,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 88094, 0, 3,
                                                                       85259, 52574, 85469,
                                                                       31349, 31499, 55049,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 88409, 0, 3,
                                                                       85469, 52724, 85679,
                                                                       31499, 31649, 55274,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 88724, 0, 3,
                                                                       85889, 53474, 86204,
                                                                       31949, 32159, 56129,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 89165, 0, 3,
                                                                       86204, 53699, 86519,
                                                                       32159, 32369, 56444,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 89606, 0, 3,
                                                                       86519, 53924, 86834,
                                                                       32369, 32579, 56759,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90047, 0, 3,
                                                                       86834, 54149, 87149,
                                                                       32579, 32789, 57074,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90488, 0, 3,
                                                                       87149, 54374, 87464,
                                                                       32789, 32999, 57389,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 90929, 0, 3,
                                                                       87464, 54599, 87779,
                                                                       32999, 33209, 57704,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 91370, 0, 3,
                                                                       87779, 54824, 88094,
                                                                       33209, 33419, 58019,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 91811, 0, 3,
                                                                       88094, 55049, 88409,
                                                                       33419, 33629, 58334,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 92252, 0, 3,
                                                                       88724, 56129, 89165,
                                                                       34049, 34329, 59489,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 92840, 0, 3,
                                                                       89165, 56444, 89606,
                                                                       34329, 34609, 59909,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 93428, 0, 3,
                                                                       89606, 56759, 90047,
                                                                       34609, 34889, 60329,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 94016, 0, 3,
                                                                       90047, 57074, 90488,
                                                                       34889, 35169, 60749,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 94604, 0, 3,
                                                                       90488, 57389, 90929,
                                                                       35169, 35449, 61169,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 95192, 0, 3,
                                                                       90929, 57704, 91370,
                                                                       35449, 35729, 61589,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 95780, 0, 3,
                                                                       91370, 58019, 91811,
                                                                       35729, 36009, 62009,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 96368, 0, 3,
                                                                       92252, 59489, 92840,
                                                                       36569, 36929, 63509,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 97124, 0, 3,
                                                                       92840, 59909, 93428,
                                                                       36929, 37289, 64049,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 97880, 0, 3,
                                                                       93428, 60329, 94016,
                                                                       37289, 37649, 64589,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 98636, 0, 3,
                                                                       94016, 60749, 94604,
                                                                       37649, 38009, 65129,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 99392, 0, 3,
                                                                       94604, 61169, 95192,
                                                                       38009, 38369, 65669,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 100148, 0, 3,
                                                                       95192, 61589, 95780,
                                                                       38369, 38729, 66209,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 100904, 0, 3,
                                                                       96368, 63509, 97124,
                                                                       39449, 39899, 68099,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 101849, 0, 3,
                                                                       97124, 64049, 97880,
                                                                       39899, 40349, 68774,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 102794, 0, 3,
                                                                       97880, 64589, 98636,
                                                                       40349, 40799, 69449,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 103739, 0, 3,
                                                                       98636, 65129, 99392,
                                                                       40799, 41249, 70124,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 104684, 0, 3,
                                                                       99392, 65669, 100148,
                                                                       41249, 41699, 70799,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 105629, 0, 3,
                                                                       100904, 68099, 101849,
                                                                       42599, 43149, 73124,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 106784, 0, 3,
                                                                       101849, 68774, 102794,
                                                                       43149, 43699, 73949,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 107939, 0, 3,
                                                                       102794, 69449, 103739,
                                                                       43699, 44249, 74774,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 109094, 0, 3,
                                                                       103739, 70124, 104684,
                                                                       44249, 44799, 75599,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 110249, 0, 3,
                                                                       105629, 73124, 106784,
                                                                       45899, 46559, 78404,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 111635, 0, 3,
                                                                       106784, 73949, 107939,
                                                                       46559, 47219, 79394,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 113021, 0, 3,
                                                                       107939, 74774, 109094,
                                                                       47219, 47879, 80384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114407, 3, 49199,
                                                                       49214, 81374, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114435, 3, 49214,
                                                                       49229, 81395, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114463, 3, 49229,
                                                                       49244, 81416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114491, 3, 49244,
                                                                       49259, 81437, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114519, 3, 49259,
                                                                       49274, 81458, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114547, 3, 49274,
                                                                       49289, 81479, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114575, 3, 49289,
                                                                       49304, 81500, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114603, 3, 49304,
                                                                       49319, 81521, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114631, 3, 49319,
                                                                       49334, 81542, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114659, 3, 49334,
                                                                       49349, 81563, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114687, 3, 49349,
                                                                       49364, 81584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114715, 3, 49364,
                                                                       49379, 81605, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 114743, 3, 49379,
                                                                       49394, 81626, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 114771, 0, 3,
                                                                       114407, 81374, 114435,
                                                                       49424, 49469, 81647,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 114855, 0, 3,
                                                                       114435, 81395, 114463,
                                                                       49469, 49514, 81710,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 114939, 0, 3,
                                                                       114463, 81416, 114491,
                                                                       49514, 49559, 81773,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115023, 0, 3,
                                                                       114491, 81437, 114519,
                                                                       49559, 49604, 81836,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115107, 0, 3,
                                                                       114519, 81458, 114547,
                                                                       49604, 49649, 81899,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115191, 0, 3,
                                                                       114547, 81479, 114575,
                                                                       49649, 49694, 81962,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115275, 0, 3,
                                                                       114575, 81500, 114603,
                                                                       49694, 49739, 82025,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115359, 0, 3,
                                                                       114603, 81521, 114631,
                                                                       49739, 49784, 82088,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115443, 0, 3,
                                                                       114631, 81542, 114659,
                                                                       49784, 49829, 82151,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115527, 0, 3,
                                                                       114659, 81563, 114687,
                                                                       49829, 49874, 82214,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115611, 0, 3,
                                                                       114687, 81584, 114715,
                                                                       49874, 49919, 82277,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 115695, 0, 3,
                                                                       114715, 81605, 114743,
                                                                       49919, 49964, 82340,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 115779, 0, 3,
                                                                       114771, 81647, 114855,
                                                                       50054, 50144, 82403,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 115947, 0, 3,
                                                                       114855, 81710, 114939,
                                                                       50144, 50234, 82529,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116115, 0, 3,
                                                                       114939, 81773, 115023,
                                                                       50234, 50324, 82655,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116283, 0, 3,
                                                                       115023, 81836, 115107,
                                                                       50324, 50414, 82781,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116451, 0, 3,
                                                                       115107, 81899, 115191,
                                                                       50414, 50504, 82907,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116619, 0, 3,
                                                                       115191, 81962, 115275,
                                                                       50504, 50594, 83033,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116787, 0, 3,
                                                                       115275, 82025, 115359,
                                                                       50594, 50684, 83159,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 116955, 0, 3,
                                                                       115359, 82088, 115443,
                                                                       50684, 50774, 83285,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 117123, 0, 3,
                                                                       115443, 82151, 115527,
                                                                       50774, 50864, 83411,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 117291, 0, 3,
                                                                       115527, 82214, 115611,
                                                                       50864, 50954, 83537,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 117459, 0, 3,
                                                                       115611, 82277, 115695,
                                                                       50954, 51044, 83663,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 117627, 0, 3,
                                                                       115779, 82403, 115947,
                                                                       51224, 51374, 83789,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 117907, 0, 3,
                                                                       115947, 82529, 116115,
                                                                       51374, 51524, 83999,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 118187, 0, 3,
                                                                       116115, 82655, 116283,
                                                                       51524, 51674, 84209,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 118467, 0, 3,
                                                                       116283, 82781, 116451,
                                                                       51674, 51824, 84419,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 118747, 0, 3,
                                                                       116451, 82907, 116619,
                                                                       51824, 51974, 84629,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 119027, 0, 3,
                                                                       116619, 83033, 116787,
                                                                       51974, 52124, 84839,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 119307, 0, 3,
                                                                       116787, 83159, 116955,
                                                                       52124, 52274, 85049,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 119587, 0, 3,
                                                                       116955, 83285, 117123,
                                                                       52274, 52424, 85259,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 119867, 0, 3,
                                                                       117123, 83411, 117291,
                                                                       52424, 52574, 85469,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 120147, 0, 3,
                                                                       117291, 83537, 117459,
                                                                       52574, 52724, 85679,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 120427, 0, 3,
                                                                       117627, 83789, 117907,
                                                                       53024, 53249, 85889,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 120847, 0, 3,
                                                                       117907, 83999, 118187,
                                                                       53249, 53474, 86204,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 121267, 0, 3,
                                                                       118187, 84209, 118467,
                                                                       53474, 53699, 86519,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 121687, 0, 3,
                                                                       118467, 84419, 118747,
                                                                       53699, 53924, 86834,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 122107, 0, 3,
                                                                       118747, 84629, 119027,
                                                                       53924, 54149, 87149,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 122527, 0, 3,
                                                                       119027, 84839, 119307,
                                                                       54149, 54374, 87464,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 122947, 0, 3,
                                                                       119307, 85049, 119587,
                                                                       54374, 54599, 87779,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 123367, 0, 3,
                                                                       119587, 85259, 119867,
                                                                       54599, 54824, 88094,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 123787, 0, 3,
                                                                       119867, 85469, 120147,
                                                                       54824, 55049, 88409,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 124207, 0, 3,
                                                                       120427, 85889, 120847,
                                                                       55499, 55814, 88724,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 124795, 0, 3,
                                                                       120847, 86204, 121267,
                                                                       55814, 56129, 89165,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 125383, 0, 3,
                                                                       121267, 86519, 121687,
                                                                       56129, 56444, 89606,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 125971, 0, 3,
                                                                       121687, 86834, 122107,
                                                                       56444, 56759, 90047,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 126559, 0, 3,
                                                                       122107, 87149, 122527,
                                                                       56759, 57074, 90488,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 127147, 0, 3,
                                                                       122527, 87464, 122947,
                                                                       57074, 57389, 90929,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 127735, 0, 3,
                                                                       122947, 87779, 123367,
                                                                       57389, 57704, 91370,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 128323, 0, 3,
                                                                       123367, 88094, 123787,
                                                                       57704, 58019, 91811,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 128911, 0, 3,
                                                                       124207, 88724, 124795,
                                                                       58649, 59069, 92252,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 129695, 0, 3,
                                                                       124795, 89165, 125383,
                                                                       59069, 59489, 92840,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 130479, 0, 3,
                                                                       125383, 89606, 125971,
                                                                       59489, 59909, 93428,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 131263, 0, 3,
                                                                       125971, 90047, 126559,
                                                                       59909, 60329, 94016,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 132047, 0, 3,
                                                                       126559, 90488, 127147,
                                                                       60329, 60749, 94604,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 132831, 0, 3,
                                                                       127147, 90929, 127735,
                                                                       60749, 61169, 95192,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 133615, 0, 3,
                                                                       127735, 91370, 128323,
                                                                       61169, 61589, 95780,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 134399, 0, 3,
                                                                       128911, 92252, 129695,
                                                                       62429, 62969, 96368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 135407, 0, 3,
                                                                       129695, 92840, 130479,
                                                                       62969, 63509, 97124,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 136415, 0, 3,
                                                                       130479, 93428, 131263,
                                                                       63509, 64049, 97880,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 137423, 0, 3,
                                                                       131263, 94016, 132047,
                                                                       64049, 64589, 98636,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 138431, 0, 3,
                                                                       132047, 94604, 132831,
                                                                       64589, 65129, 99392,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 139439, 0, 3,
                                                                       132831, 95192, 133615,
                                                                       65129, 65669, 100148,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 140447, 0, 3,
                                                                       134399, 96368, 135407,
                                                                       66749, 67424, 100904,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 141707, 0, 3,
                                                                       135407, 97124, 136415,
                                                                       67424, 68099, 101849,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 142967, 0, 3,
                                                                       136415, 97880, 137423,
                                                                       68099, 68774, 102794,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 144227, 0, 3,
                                                                       137423, 98636, 138431,
                                                                       68774, 69449, 103739,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 145487, 0, 3,
                                                                       138431, 99392, 139439,
                                                                       69449, 70124, 104684,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 146747, 0, 3,
                                                                       140447, 100904, 141707,
                                                                       71474, 72299, 105629,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 148287, 0, 3,
                                                                       141707, 101849, 142967,
                                                                       72299, 73124, 106784,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 149827, 0, 3,
                                                                       142967, 102794, 144227,
                                                                       73124, 73949, 107939,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 151367, 0, 3,
                                                                       144227, 103739, 145487,
                                                                       73949, 74774, 109094,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 152907, 0, 3,
                                                                       146747, 105629, 148287,
                                                                       76424, 77414, 110249,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 154755, 0, 3,
                                                                       148287, 106784, 149827,
                                                                       77414, 78404, 111635,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 156603, 0, 3,
                                                                       149827, 107939, 151367,
                                                                       78404, 79394, 113021,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158451, 3, 81374,
                                                                       81395, 114463, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158487, 3, 81395,
                                                                       81416, 114491, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158523, 3, 81416,
                                                                       81437, 114519, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158559, 3, 81437,
                                                                       81458, 114547, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158595, 3, 81458,
                                                                       81479, 114575, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158631, 3, 81479,
                                                                       81500, 114603, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158667, 3, 81500,
                                                                       81521, 114631, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158703, 3, 81521,
                                                                       81542, 114659, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158739, 3, 81542,
                                                                       81563, 114687, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158775, 3, 81563,
                                                                       81584, 114715, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 158811, 3, 81584,
                                                                       81605, 114743, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 158847, 0, 3,
                                                                       158451, 114463, 158487,
                                                                       81647, 81710, 114939,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 158955, 0, 3,
                                                                       158487, 114491, 158523,
                                                                       81710, 81773, 115023,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159063, 0, 3,
                                                                       158523, 114519, 158559,
                                                                       81773, 81836, 115107,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159171, 0, 3,
                                                                       158559, 114547, 158595,
                                                                       81836, 81899, 115191,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159279, 0, 3,
                                                                       158595, 114575, 158631,
                                                                       81899, 81962, 115275,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159387, 0, 3,
                                                                       158631, 114603, 158667,
                                                                       81962, 82025, 115359,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159495, 0, 3,
                                                                       158667, 114631, 158703,
                                                                       82025, 82088, 115443,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159603, 0, 3,
                                                                       158703, 114659, 158739,
                                                                       82088, 82151, 115527,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159711, 0, 3,
                                                                       158739, 114687, 158775,
                                                                       82151, 82214, 115611,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 159819, 0, 3,
                                                                       158775, 114715, 158811,
                                                                       82214, 82277, 115695,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 159927, 0, 3,
                                                                       158847, 114939, 158955,
                                                                       82403, 82529, 116115,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 160143, 0, 3,
                                                                       158955, 115023, 159063,
                                                                       82529, 82655, 116283,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 160359, 0, 3,
                                                                       159063, 115107, 159171,
                                                                       82655, 82781, 116451,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 160575, 0, 3,
                                                                       159171, 115191, 159279,
                                                                       82781, 82907, 116619,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 160791, 0, 3,
                                                                       159279, 115275, 159387,
                                                                       82907, 83033, 116787,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 161007, 0, 3,
                                                                       159387, 115359, 159495,
                                                                       83033, 83159, 116955,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 161223, 0, 3,
                                                                       159495, 115443, 159603,
                                                                       83159, 83285, 117123,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 161439, 0, 3,
                                                                       159603, 115527, 159711,
                                                                       83285, 83411, 117291,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 161655, 0, 3,
                                                                       159711, 115611, 159819,
                                                                       83411, 83537, 117459,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 161871, 0, 3,
                                                                       159927, 116115, 160143,
                                                                       83789, 83999, 118187,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 162231, 0, 3,
                                                                       160143, 116283, 160359,
                                                                       83999, 84209, 118467,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 162591, 0, 3,
                                                                       160359, 116451, 160575,
                                                                       84209, 84419, 118747,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 162951, 0, 3,
                                                                       160575, 116619, 160791,
                                                                       84419, 84629, 119027,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 163311, 0, 3,
                                                                       160791, 116787, 161007,
                                                                       84629, 84839, 119307,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 163671, 0, 3,
                                                                       161007, 116955, 161223,
                                                                       84839, 85049, 119587,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 164031, 0, 3,
                                                                       161223, 117123, 161439,
                                                                       85049, 85259, 119867,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 164391, 0, 3,
                                                                       161439, 117291, 161655,
                                                                       85259, 85469, 120147,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 164751, 0, 3,
                                                                       161871, 118187, 162231,
                                                                       85889, 86204, 121267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 165291, 0, 3,
                                                                       162231, 118467, 162591,
                                                                       86204, 86519, 121687,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 165831, 0, 3,
                                                                       162591, 118747, 162951,
                                                                       86519, 86834, 122107,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 166371, 0, 3,
                                                                       162951, 119027, 163311,
                                                                       86834, 87149, 122527,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 166911, 0, 3,
                                                                       163311, 119307, 163671,
                                                                       87149, 87464, 122947,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 167451, 0, 3,
                                                                       163671, 119587, 164031,
                                                                       87464, 87779, 123367,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 167991, 0, 3,
                                                                       164031, 119867, 164391,
                                                                       87779, 88094, 123787,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 168531, 0, 3,
                                                                       164751, 121267, 165291,
                                                                       88724, 89165, 125383,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 169287, 0, 3,
                                                                       165291, 121687, 165831,
                                                                       89165, 89606, 125971,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 170043, 0, 3,
                                                                       165831, 122107, 166371,
                                                                       89606, 90047, 126559,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 170799, 0, 3,
                                                                       166371, 122527, 166911,
                                                                       90047, 90488, 127147,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 171555, 0, 3,
                                                                       166911, 122947, 167451,
                                                                       90488, 90929, 127735,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 172311, 0, 3,
                                                                       167451, 123367, 167991,
                                                                       90929, 91370, 128323,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 173067, 0, 3,
                                                                       168531, 125383, 169287,
                                                                       92252, 92840, 130479,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 174075, 0, 3,
                                                                       169287, 125971, 170043,
                                                                       92840, 93428, 131263,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 175083, 0, 3,
                                                                       170043, 126559, 170799,
                                                                       93428, 94016, 132047,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 176091, 0, 3,
                                                                       170799, 127147, 171555,
                                                                       94016, 94604, 132831,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 177099, 0, 3,
                                                                       171555, 127735, 172311,
                                                                       94604, 95192, 133615,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 178107, 0, 3,
                                                                       173067, 130479, 174075,
                                                                       96368, 97124, 136415,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 179403, 0, 3,
                                                                       174075, 131263, 175083,
                                                                       97124, 97880, 137423,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 180699, 0, 3,
                                                                       175083, 132047, 176091,
                                                                       97880, 98636, 138431,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 181995, 0, 3,
                                                                       176091, 132831, 177099,
                                                                       98636, 99392, 139439,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 183291, 0, 3,
                                                                       178107, 136415, 179403,
                                                                       100904, 101849, 142967,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 184911, 0, 3,
                                                                       179403, 137423, 180699,
                                                                       101849, 102794, 144227,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 186531, 0, 3,
                                                                       180699, 138431, 181995,
                                                                       102794, 103739, 145487,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 188151, 0, 3,
                                                                       183291, 142967, 184911,
                                                                       105629, 106784, 149827,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 190131, 0, 3,
                                                                       184911, 144227, 186531,
                                                                       106784, 107939, 151367,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 192111, 0, 3,
                                                                       188151, 149827, 190131,
                                                                       110249, 111635, 156603,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194487, 3, 114407,
                                                                       114435, 158451, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194532, 3, 114435,
                                                                       114463, 158487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194577, 3, 114463,
                                                                       114491, 158523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194622, 3, 114491,
                                                                       114519, 158559, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194667, 3, 114519,
                                                                       114547, 158595, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194712, 3, 114547,
                                                                       114575, 158631, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194757, 3, 114575,
                                                                       114603, 158667, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194802, 3, 114603,
                                                                       114631, 158703, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194847, 3, 114631,
                                                                       114659, 158739, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194892, 3, 114659,
                                                                       114687, 158775, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 194937, 3, 114687,
                                                                       114715, 158811, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 194982, 0, 3,
                                                                       194487, 158451, 194532,
                                                                       114771, 114855, 158847,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195117, 0, 3,
                                                                       194532, 158487, 194577,
                                                                       114855, 114939, 158955,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195252, 0, 3,
                                                                       194577, 158523, 194622,
                                                                       114939, 115023, 159063,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195387, 0, 3,
                                                                       194622, 158559, 194667,
                                                                       115023, 115107, 159171,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195522, 0, 3,
                                                                       194667, 158595, 194712,
                                                                       115107, 115191, 159279,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195657, 0, 3,
                                                                       194712, 158631, 194757,
                                                                       115191, 115275, 159387,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195792, 0, 3,
                                                                       194757, 158667, 194802,
                                                                       115275, 115359, 159495,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 195927, 0, 3,
                                                                       194802, 158703, 194847,
                                                                       115359, 115443, 159603,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 196062, 0, 3,
                                                                       194847, 158739, 194892,
                                                                       115443, 115527, 159711,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 196197, 0, 3,
                                                                       194892, 158775, 194937,
                                                                       115527, 115611, 159819,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 196332, 0, 3,
                                                                       194982, 158847, 195117,
                                                                       115779, 115947, 159927,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 196602, 0, 3,
                                                                       195117, 158955, 195252,
                                                                       115947, 116115, 160143,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 196872, 0, 3,
                                                                       195252, 159063, 195387,
                                                                       116115, 116283, 160359,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 197142, 0, 3,
                                                                       195387, 159171, 195522,
                                                                       116283, 116451, 160575,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 197412, 0, 3,
                                                                       195522, 159279, 195657,
                                                                       116451, 116619, 160791,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 197682, 0, 3,
                                                                       195657, 159387, 195792,
                                                                       116619, 116787, 161007,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 197952, 0, 3,
                                                                       195792, 159495, 195927,
                                                                       116787, 116955, 161223,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 198222, 0, 3,
                                                                       195927, 159603, 196062,
                                                                       116955, 117123, 161439,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 198492, 0, 3,
                                                                       196062, 159711, 196197,
                                                                       117123, 117291, 161655,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 198762, 0, 3,
                                                                       196332, 159927, 196602,
                                                                       117627, 117907, 161871,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 199212, 0, 3,
                                                                       196602, 160143, 196872,
                                                                       117907, 118187, 162231,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 199662, 0, 3,
                                                                       196872, 160359, 197142,
                                                                       118187, 118467, 162591,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 200112, 0, 3,
                                                                       197142, 160575, 197412,
                                                                       118467, 118747, 162951,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 200562, 0, 3,
                                                                       197412, 160791, 197682,
                                                                       118747, 119027, 163311,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 201012, 0, 3,
                                                                       197682, 161007, 197952,
                                                                       119027, 119307, 163671,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 201462, 0, 3,
                                                                       197952, 161223, 198222,
                                                                       119307, 119587, 164031,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 201912, 0, 3,
                                                                       198222, 161439, 198492,
                                                                       119587, 119867, 164391,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 202362, 0, 3,
                                                                       198762, 161871, 199212,
                                                                       120427, 120847, 164751,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 203037, 0, 3,
                                                                       199212, 162231, 199662,
                                                                       120847, 121267, 165291,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 203712, 0, 3,
                                                                       199662, 162591, 200112,
                                                                       121267, 121687, 165831,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 204387, 0, 3,
                                                                       200112, 162951, 200562,
                                                                       121687, 122107, 166371,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 205062, 0, 3,
                                                                       200562, 163311, 201012,
                                                                       122107, 122527, 166911,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 205737, 0, 3,
                                                                       201012, 163671, 201462,
                                                                       122527, 122947, 167451,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 206412, 0, 3,
                                                                       201462, 164031, 201912,
                                                                       122947, 123367, 167991,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 207087, 0, 3,
                                                                       202362, 164751, 203037,
                                                                       124207, 124795, 168531,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 208032, 0, 3,
                                                                       203037, 165291, 203712,
                                                                       124795, 125383, 169287,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 208977, 0, 3,
                                                                       203712, 165831, 204387,
                                                                       125383, 125971, 170043,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 209922, 0, 3,
                                                                       204387, 166371, 205062,
                                                                       125971, 126559, 170799,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 210867, 0, 3,
                                                                       205062, 166911, 205737,
                                                                       126559, 127147, 171555,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 211812, 0, 3,
                                                                       205737, 167451, 206412,
                                                                       127147, 127735, 172311,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 212757, 0, 3,
                                                                       207087, 168531, 208032,
                                                                       128911, 129695, 173067,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 214017, 0, 3,
                                                                       208032, 169287, 208977,
                                                                       129695, 130479, 174075,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 215277, 0, 3,
                                                                       208977, 170043, 209922,
                                                                       130479, 131263, 175083,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 216537, 0, 3,
                                                                       209922, 170799, 210867,
                                                                       131263, 132047, 176091,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 217797, 0, 3,
                                                                       210867, 171555, 211812,
                                                                       132047, 132831, 177099,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 219057, 0, 3,
                                                                       212757, 173067, 214017,
                                                                       134399, 135407, 178107,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 220677, 0, 3,
                                                                       214017, 174075, 215277,
                                                                       135407, 136415, 179403,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 222297, 0, 3,
                                                                       215277, 175083, 216537,
                                                                       136415, 137423, 180699,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 223917, 0, 3,
                                                                       216537, 176091, 217797,
                                                                       137423, 138431, 181995,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 225537, 0, 3,
                                                                       219057, 178107, 220677,
                                                                       140447, 141707, 183291,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 227562, 0, 3,
                                                                       220677, 179403, 222297,
                                                                       141707, 142967, 184911,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 229587, 0, 3,
                                                                       222297, 180699, 223917,
                                                                       142967, 144227, 186531,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 231612, 0, 3,
                                                                       225537, 183291, 227562,
                                                                       146747, 148287, 188151,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 234087, 0, 3,
                                                                       227562, 184911, 229587,
                                                                       148287, 149827, 190131,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 236562, 0, 3,
                                                                       231612, 188151, 234087,
                                                                       152907, 154755, 192111,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 239532, 207087, 945, ncols);

                    simdfunc::contract_primitives(buffer, 240834, 212757, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 242570, 219057, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 244802, 225537, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 247592, 231612, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 251002, 236562, 2970, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 240477, 239532, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 242094, 240834, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 244190, 242570, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 246827, 244802, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 250067, 247592, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 253972, 251002, 66, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 255094, 240477, 242094, 17,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 256165, 242094, 244190, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 257593, 244190, 246827, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 259429, 246827, 250067, 17,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 261724, 250067, 253972, 17,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 264529, 255094, 256165, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 266671, 256165, 257593, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 269527, 257593, 259429, 17,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 273199, 259429, 261724, 17,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 277789, 264529, 266671, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 281359, 266671, 269527, 17,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 286119, 269527, 273199, 17,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 292239, 277789, 281359, 17,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 297594, 281359, 286119, 17,
                                             nmax);

        simdtrf::compute_hrr_hh_out_of_first(buffer, coordinates, 304734, 292239, 297594, 17,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 312231, 304734, 21, 17, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 312231, 187, nmax);
    }

    for (size_t m = 0; m < 2057; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
