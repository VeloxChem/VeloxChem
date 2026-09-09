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


#include "SimdThreeCenterElectronRepulsionRecIII.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIH.hpp"
#include "SimdTransferII.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKH.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLG.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMF.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferND.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransferOP.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iii_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iii_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 316029, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2197 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 316029, 194747, 15176, dimensions);

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

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3296, 0, 3, 2152,
                                                                       2207, 2702, 2768, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3374, 0, 3, 2207,
                                                                       2262, 2768, 2834, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3452, 0, 3, 2262,
                                                                       2317, 2834, 2900, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3530, 0, 3, 2317,
                                                                       2372, 2900, 2966, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3608, 0, 3, 2372,
                                                                       2427, 2966, 3032, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3686, 0, 3, 2427,
                                                                       2482, 3032, 3098, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3764, 0, 3, 2482,
                                                                       2537, 3098, 3164, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3842, 0, 3, 2537,
                                                                       2592, 3164, 3230, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 3920, 0, 3, 2702,
                                                                       2768, 3296, 3374, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4011, 0, 3, 2768,
                                                                       2834, 3374, 3452, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4102, 0, 3, 2834,
                                                                       2900, 3452, 3530, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4193, 0, 3, 2900,
                                                                       2966, 3530, 3608, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4284, 0, 3, 2966,
                                                                       3032, 3608, 3686, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4375, 0, 3, 3032,
                                                                       3098, 3686, 3764, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4466, 0, 3, 3098,
                                                                       3164, 3764, 3842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4557, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4560, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4563, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4566, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4569, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4572, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4575, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4578, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4581, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4584, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4587, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4590, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4593, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4596, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4599, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4602, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4605, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4608, 3, 9, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4617, 3, 10, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4626, 3, 11, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4635, 3, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4644, 3, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4653, 3, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4662, 3, 15, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4671, 3, 16, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4680, 3, 17, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4689, 3, 18, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4698, 3, 19, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4707, 3, 20, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4716, 3, 21, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4725, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4734, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4743, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4752, 3, 32, 92,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4770, 3, 35, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4788, 3, 38, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4806, 3, 41, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4824, 3, 44, 116,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4842, 3, 47, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4860, 3, 50, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4878, 3, 53, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4896, 3, 56, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4914, 3, 59, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4932, 3, 62, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4950, 3, 65, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4968, 3, 68, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4986, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5004, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5022, 3, 92, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5052, 3, 98, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5082, 3, 104, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5112, 3, 110, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5142, 3, 116, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5172, 3, 122, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5202, 3, 128, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5232, 3, 134, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5262, 3, 140, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5292, 3, 146, 292,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5322, 3, 152, 302,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5352, 3, 158, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5382, 3, 164, 322,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5412, 3, 170, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5442, 3, 202, 372,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5487, 3, 212, 387,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5532, 3, 222, 402,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5577, 3, 232, 417,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5622, 3, 242, 432,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5667, 3, 252, 447,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5712, 3, 262, 462,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5757, 3, 272, 477,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5802, 3, 282, 492,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5847, 3, 292, 507,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5892, 3, 302, 522,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5937, 3, 312, 537,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5982, 3, 322, 552,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6027, 3, 372, 609,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6090, 3, 387, 630,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6153, 3, 402, 651,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6216, 3, 417, 672,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6279, 3, 432, 693,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6342, 3, 447, 714,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6405, 3, 462, 735,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6468, 3, 477, 756,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6531, 3, 492, 777,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6594, 3, 507, 798,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6657, 3, 522, 819,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6720, 3, 537, 840,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6783, 3, 609, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6867, 3, 630, 945,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6951, 3, 651, 973,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7035, 3, 672,
                                                                       1001, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7119, 3, 693,
                                                                       1029, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7203, 3, 714,
                                                                       1057, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7287, 3, 735,
                                                                       1085, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7371, 3, 756,
                                                                       1113, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7455, 3, 777,
                                                                       1141, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7539, 3, 798,
                                                                       1169, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7623, 3, 819,
                                                                       1197, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7707, 3, 917,
                                                                       1297, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7815, 3, 945,
                                                                       1333, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7923, 3, 973,
                                                                       1369, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8031, 3, 1001,
                                                                       1405, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8139, 3, 1029,
                                                                       1441, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8247, 3, 1057,
                                                                       1477, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8355, 3, 1085,
                                                                       1513, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8463, 3, 1113,
                                                                       1549, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8571, 3, 1141,
                                                                       1585, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8679, 3, 1169,
                                                                       1621, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8787, 3, 1297,
                                                                       1747, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8922, 3, 1333,
                                                                       1792, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9057, 3, 1369,
                                                                       1837, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9192, 3, 1405,
                                                                       1882, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9327, 3, 1441,
                                                                       1927, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9462, 3, 1477,
                                                                       1972, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9597, 3, 1513,
                                                                       2017, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9732, 3, 1549,
                                                                       2062, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9867, 3, 1585,
                                                                       2107, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10002, 3, 1747,
                                                                       2262, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10167, 3, 1792,
                                                                       2317, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10332, 3, 1837,
                                                                       2372, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10497, 3, 1882,
                                                                       2427, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10662, 3, 1927,
                                                                       2482, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10827, 3, 1972,
                                                                       2537, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10992, 3, 2017,
                                                                       2592, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11157, 3, 2062,
                                                                       2647, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11322, 3, 2262,
                                                                       2834, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11520, 3, 2317,
                                                                       2900, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11718, 3, 2372,
                                                                       2966, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 11916, 3, 2427,
                                                                       3032, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12114, 3, 2482,
                                                                       3098, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12312, 3, 2537,
                                                                       3164, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12510, 3, 2592,
                                                                       3230, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12708, 3, 2834,
                                                                       3452, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 12942, 3, 2900,
                                                                       3530, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13176, 3, 2966,
                                                                       3608, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13410, 3, 3032,
                                                                       3686, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13644, 3, 3098,
                                                                       3764, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 13878, 3, 3164,
                                                                       3842, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14112, 3, 3452,
                                                                       4102, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14385, 3, 3530,
                                                                       4193, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14658, 3, 3608,
                                                                       4284, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 14931, 3, 3686,
                                                                       4375, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 15204, 3, 3764,
                                                                       4466, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15477, 3, 7, 8,
                                                                       4557, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15483, 3, 8, 9,
                                                                       4560, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15489, 3, 9, 10,
                                                                       4563, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15495, 3, 10, 11,
                                                                       4566, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15501, 3, 11, 12,
                                                                       4569, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15507, 3, 12, 13,
                                                                       4572, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15513, 3, 13, 14,
                                                                       4575, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15519, 3, 14, 15,
                                                                       4578, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15525, 3, 15, 16,
                                                                       4581, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15531, 3, 16, 17,
                                                                       4584, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15537, 3, 17, 18,
                                                                       4587, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15543, 3, 18, 19,
                                                                       4590, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15549, 3, 19, 20,
                                                                       4593, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15555, 3, 20, 21,
                                                                       4596, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15561, 3, 21, 22,
                                                                       4599, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15567, 3, 22, 23,
                                                                       4602, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15573, 3, 23, 24,
                                                                       4605, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15579, 0, 3,
                                                                       15477, 4557, 15483, 4608,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15597, 0, 3,
                                                                       15483, 4560, 15489, 4617,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15615, 0, 3,
                                                                       15489, 4563, 15495, 4626,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15633, 0, 3,
                                                                       15495, 4566, 15501, 4635,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15651, 0, 3,
                                                                       15501, 4569, 15507, 4644,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15669, 0, 3,
                                                                       15507, 4572, 15513, 4653,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15687, 0, 3,
                                                                       15513, 4575, 15519, 4662,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15705, 0, 3,
                                                                       15519, 4578, 15525, 4671,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15723, 0, 3,
                                                                       15525, 4581, 15531, 4680,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15741, 0, 3,
                                                                       15531, 4584, 15537, 4689,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15759, 0, 3,
                                                                       15537, 4587, 15543, 4698,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15777, 0, 3,
                                                                       15543, 4590, 15549, 4707,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15795, 0, 3,
                                                                       15549, 4593, 15555, 4716,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15813, 0, 3,
                                                                       15555, 4596, 15561, 4725,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15831, 0, 3,
                                                                       15561, 4599, 15567, 4734,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 15849, 0, 3,
                                                                       15567, 4602, 15573, 4743,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15867, 0, 3,
                                                                       15579, 4608, 15597, 80,
                                                                       86, 4752, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15903, 0, 3,
                                                                       15597, 4617, 15615, 86,
                                                                       92, 4770, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15939, 0, 3,
                                                                       15615, 4626, 15633, 92,
                                                                       98, 4788, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 15975, 0, 3,
                                                                       15633, 4635, 15651, 98,
                                                                       104, 4806, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16011, 0, 3,
                                                                       15651, 4644, 15669, 104,
                                                                       110, 4824, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16047, 0, 3,
                                                                       15669, 4653, 15687, 110,
                                                                       116, 4842, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16083, 0, 3,
                                                                       15687, 4662, 15705, 116,
                                                                       122, 4860, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16119, 0, 3,
                                                                       15705, 4671, 15723, 122,
                                                                       128, 4878, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16155, 0, 3,
                                                                       15723, 4680, 15741, 128,
                                                                       134, 4896, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16191, 0, 3,
                                                                       15741, 4689, 15759, 134,
                                                                       140, 4914, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16227, 0, 3,
                                                                       15759, 4698, 15777, 140,
                                                                       146, 4932, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16263, 0, 3,
                                                                       15777, 4707, 15795, 146,
                                                                       152, 4950, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16299, 0, 3,
                                                                       15795, 4716, 15813, 152,
                                                                       158, 4968, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16335, 0, 3,
                                                                       15813, 4725, 15831, 158,
                                                                       164, 4986, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 16371, 0, 3,
                                                                       15831, 4734, 15849, 164,
                                                                       170, 5004, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16407, 0, 3,
                                                                       15867, 4752, 15903, 182,
                                                                       192, 5022, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16467, 0, 3,
                                                                       15903, 4770, 15939, 192,
                                                                       202, 5052, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16527, 0, 3,
                                                                       15939, 4788, 15975, 202,
                                                                       212, 5082, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16587, 0, 3,
                                                                       15975, 4806, 16011, 212,
                                                                       222, 5112, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16647, 0, 3,
                                                                       16011, 4824, 16047, 222,
                                                                       232, 5142, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16707, 0, 3,
                                                                       16047, 4842, 16083, 232,
                                                                       242, 5172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16767, 0, 3,
                                                                       16083, 4860, 16119, 242,
                                                                       252, 5202, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16827, 0, 3,
                                                                       16119, 4878, 16155, 252,
                                                                       262, 5232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16887, 0, 3,
                                                                       16155, 4896, 16191, 262,
                                                                       272, 5262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 16947, 0, 3,
                                                                       16191, 4914, 16227, 272,
                                                                       282, 5292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17007, 0, 3,
                                                                       16227, 4932, 16263, 282,
                                                                       292, 5322, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17067, 0, 3,
                                                                       16263, 4950, 16299, 292,
                                                                       302, 5352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17127, 0, 3,
                                                                       16299, 4968, 16335, 302,
                                                                       312, 5382, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 17187, 0, 3,
                                                                       16335, 4986, 16371, 312,
                                                                       322, 5412, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17247, 0, 3,
                                                                       16407, 5022, 16467, 342,
                                                                       357, 5442, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17337, 0, 3,
                                                                       16467, 5052, 16527, 357,
                                                                       372, 5487, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17427, 0, 3,
                                                                       16527, 5082, 16587, 372,
                                                                       387, 5532, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17517, 0, 3,
                                                                       16587, 5112, 16647, 387,
                                                                       402, 5577, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17607, 0, 3,
                                                                       16647, 5142, 16707, 402,
                                                                       417, 5622, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17697, 0, 3,
                                                                       16707, 5172, 16767, 417,
                                                                       432, 5667, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17787, 0, 3,
                                                                       16767, 5202, 16827, 432,
                                                                       447, 5712, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17877, 0, 3,
                                                                       16827, 5232, 16887, 447,
                                                                       462, 5757, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 17967, 0, 3,
                                                                       16887, 5262, 16947, 462,
                                                                       477, 5802, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18057, 0, 3,
                                                                       16947, 5292, 17007, 477,
                                                                       492, 5847, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18147, 0, 3,
                                                                       17007, 5322, 17067, 492,
                                                                       507, 5892, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18237, 0, 3,
                                                                       17067, 5352, 17127, 507,
                                                                       522, 5937, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 18327, 0, 3,
                                                                       17127, 5382, 17187, 522,
                                                                       537, 5982, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18417, 0, 3,
                                                                       17247, 5442, 17337, 567,
                                                                       588, 6027, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18543, 0, 3,
                                                                       17337, 5487, 17427, 588,
                                                                       609, 6090, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18669, 0, 3,
                                                                       17427, 5532, 17517, 609,
                                                                       630, 6153, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18795, 0, 3,
                                                                       17517, 5577, 17607, 630,
                                                                       651, 6216, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 18921, 0, 3,
                                                                       17607, 5622, 17697, 651,
                                                                       672, 6279, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19047, 0, 3,
                                                                       17697, 5667, 17787, 672,
                                                                       693, 6342, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19173, 0, 3,
                                                                       17787, 5712, 17877, 693,
                                                                       714, 6405, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19299, 0, 3,
                                                                       17877, 5757, 17967, 714,
                                                                       735, 6468, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19425, 0, 3,
                                                                       17967, 5802, 18057, 735,
                                                                       756, 6531, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19551, 0, 3,
                                                                       18057, 5847, 18147, 756,
                                                                       777, 6594, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19677, 0, 3,
                                                                       18147, 5892, 18237, 777,
                                                                       798, 6657, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 19803, 0, 3,
                                                                       18237, 5937, 18327, 798,
                                                                       819, 6720, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 19929, 0, 3,
                                                                       18417, 6027, 18543, 861,
                                                                       889, 6783, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20097, 0, 3,
                                                                       18543, 6090, 18669, 889,
                                                                       917, 6867, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20265, 0, 3,
                                                                       18669, 6153, 18795, 917,
                                                                       945, 6951, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20433, 0, 3,
                                                                       18795, 6216, 18921, 945,
                                                                       973, 7035, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20601, 0, 3,
                                                                       18921, 6279, 19047, 973,
                                                                       1001, 7119, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20769, 0, 3,
                                                                       19047, 6342, 19173, 1001,
                                                                       1029, 7203, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 20937, 0, 3,
                                                                       19173, 6405, 19299, 1029,
                                                                       1057, 7287, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21105, 0, 3,
                                                                       19299, 6468, 19425, 1057,
                                                                       1085, 7371, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21273, 0, 3,
                                                                       19425, 6531, 19551, 1085,
                                                                       1113, 7455, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21441, 0, 3,
                                                                       19551, 6594, 19677, 1113,
                                                                       1141, 7539, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 21609, 0, 3,
                                                                       19677, 6657, 19803, 1141,
                                                                       1169, 7623, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21777, 0, 3,
                                                                       19929, 6783, 20097, 1225,
                                                                       1261, 7707, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 21993, 0, 3,
                                                                       20097, 6867, 20265, 1261,
                                                                       1297, 7815, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22209, 0, 3,
                                                                       20265, 6951, 20433, 1297,
                                                                       1333, 7923, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22425, 0, 3,
                                                                       20433, 7035, 20601, 1333,
                                                                       1369, 8031, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22641, 0, 3,
                                                                       20601, 7119, 20769, 1369,
                                                                       1405, 8139, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 22857, 0, 3,
                                                                       20769, 7203, 20937, 1405,
                                                                       1441, 8247, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23073, 0, 3,
                                                                       20937, 7287, 21105, 1441,
                                                                       1477, 8355, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23289, 0, 3,
                                                                       21105, 7371, 21273, 1477,
                                                                       1513, 8463, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23505, 0, 3,
                                                                       21273, 7455, 21441, 1513,
                                                                       1549, 8571, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 23721, 0, 3,
                                                                       21441, 7539, 21609, 1549,
                                                                       1585, 8679, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23937, 0, 3,
                                                                       21777, 7707, 21993, 1657,
                                                                       1702, 8787, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24207, 0, 3,
                                                                       21993, 7815, 22209, 1702,
                                                                       1747, 8922, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24477, 0, 3,
                                                                       22209, 7923, 22425, 1747,
                                                                       1792, 9057, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 24747, 0, 3,
                                                                       22425, 8031, 22641, 1792,
                                                                       1837, 9192, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25017, 0, 3,
                                                                       22641, 8139, 22857, 1837,
                                                                       1882, 9327, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25287, 0, 3,
                                                                       22857, 8247, 23073, 1882,
                                                                       1927, 9462, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25557, 0, 3,
                                                                       23073, 8355, 23289, 1927,
                                                                       1972, 9597, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 25827, 0, 3,
                                                                       23289, 8463, 23505, 1972,
                                                                       2017, 9732, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26097, 0, 3,
                                                                       23505, 8571, 23721, 2017,
                                                                       2062, 9867, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 26367, 0, 3,
                                                                       23937, 8787, 24207, 2152,
                                                                       2207, 10002, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 26697, 0, 3,
                                                                       24207, 8922, 24477, 2207,
                                                                       2262, 10167, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 27027, 0, 3,
                                                                       24477, 9057, 24747, 2262,
                                                                       2317, 10332, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 27357, 0, 3,
                                                                       24747, 9192, 25017, 2317,
                                                                       2372, 10497, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 27687, 0, 3,
                                                                       25017, 9327, 25287, 2372,
                                                                       2427, 10662, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28017, 0, 3,
                                                                       25287, 9462, 25557, 2427,
                                                                       2482, 10827, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28347, 0, 3,
                                                                       25557, 9597, 25827, 2482,
                                                                       2537, 10992, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 28677, 0, 3,
                                                                       25827, 9732, 26097, 2537,
                                                                       2592, 11157, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 29007, 0, 3,
                                                                       26367, 10002, 26697, 2702,
                                                                       2768, 11322, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 29403, 0, 3,
                                                                       26697, 10167, 27027, 2768,
                                                                       2834, 11520, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 29799, 0, 3,
                                                                       27027, 10332, 27357, 2834,
                                                                       2900, 11718, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30195, 0, 3,
                                                                       27357, 10497, 27687, 2900,
                                                                       2966, 11916, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30591, 0, 3,
                                                                       27687, 10662, 28017, 2966,
                                                                       3032, 12114, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 30987, 0, 3,
                                                                       28017, 10827, 28347, 3032,
                                                                       3098, 12312, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31383, 0, 3,
                                                                       28347, 10992, 28677, 3098,
                                                                       3164, 12510, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 31779, 0, 3,
                                                                       29007, 11322, 29403, 3296,
                                                                       3374, 12708, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 32247, 0, 3,
                                                                       29403, 11520, 29799, 3374,
                                                                       3452, 12942, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 32715, 0, 3,
                                                                       29799, 11718, 30195, 3452,
                                                                       3530, 13176, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 33183, 0, 3,
                                                                       30195, 11916, 30591, 3530,
                                                                       3608, 13410, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 33651, 0, 3,
                                                                       30591, 12114, 30987, 3608,
                                                                       3686, 13644, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 34119, 0, 3,
                                                                       30987, 12312, 31383, 3686,
                                                                       3764, 13878, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 34587, 0, 3,
                                                                       31779, 12708, 32247, 3920,
                                                                       4011, 14112, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 35133, 0, 3,
                                                                       32247, 12942, 32715, 4011,
                                                                       4102, 14385, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 35679, 0, 3,
                                                                       32715, 13176, 33183, 4102,
                                                                       4193, 14658, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 36225, 0, 3,
                                                                       33183, 13410, 33651, 4193,
                                                                       4284, 14931, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 36771, 0, 3,
                                                                       33651, 13644, 34119, 4284,
                                                                       4375, 15204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37317, 3, 4557,
                                                                       4560, 15489, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37327, 3, 4560,
                                                                       4563, 15495, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37337, 3, 4563,
                                                                       4566, 15501, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37347, 3, 4566,
                                                                       4569, 15507, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37357, 3, 4569,
                                                                       4572, 15513, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37367, 3, 4572,
                                                                       4575, 15519, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37377, 3, 4575,
                                                                       4578, 15525, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37387, 3, 4578,
                                                                       4581, 15531, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37397, 3, 4581,
                                                                       4584, 15537, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37407, 3, 4584,
                                                                       4587, 15543, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37417, 3, 4587,
                                                                       4590, 15549, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37427, 3, 4590,
                                                                       4593, 15555, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37437, 3, 4593,
                                                                       4596, 15561, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37447, 3, 4596,
                                                                       4599, 15567, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 37457, 3, 4599,
                                                                       4602, 15573, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37467, 0, 3,
                                                                       37317, 15489, 37327, 4608,
                                                                       4617, 15615, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37497, 0, 3,
                                                                       37327, 15495, 37337, 4617,
                                                                       4626, 15633, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37527, 0, 3,
                                                                       37337, 15501, 37347, 4626,
                                                                       4635, 15651, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37557, 0, 3,
                                                                       37347, 15507, 37357, 4635,
                                                                       4644, 15669, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37587, 0, 3,
                                                                       37357, 15513, 37367, 4644,
                                                                       4653, 15687, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37617, 0, 3,
                                                                       37367, 15519, 37377, 4653,
                                                                       4662, 15705, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37647, 0, 3,
                                                                       37377, 15525, 37387, 4662,
                                                                       4671, 15723, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37677, 0, 3,
                                                                       37387, 15531, 37397, 4671,
                                                                       4680, 15741, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37707, 0, 3,
                                                                       37397, 15537, 37407, 4680,
                                                                       4689, 15759, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37737, 0, 3,
                                                                       37407, 15543, 37417, 4689,
                                                                       4698, 15777, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37767, 0, 3,
                                                                       37417, 15549, 37427, 4698,
                                                                       4707, 15795, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37797, 0, 3,
                                                                       37427, 15555, 37437, 4707,
                                                                       4716, 15813, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37827, 0, 3,
                                                                       37437, 15561, 37447, 4716,
                                                                       4725, 15831, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 37857, 0, 3,
                                                                       37447, 15567, 37457, 4725,
                                                                       4734, 15849, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37887, 0, 3,
                                                                       37467, 15615, 37497, 4752,
                                                                       4770, 15939, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 37947, 0, 3,
                                                                       37497, 15633, 37527, 4770,
                                                                       4788, 15975, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38007, 0, 3,
                                                                       37527, 15651, 37557, 4788,
                                                                       4806, 16011, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38067, 0, 3,
                                                                       37557, 15669, 37587, 4806,
                                                                       4824, 16047, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38127, 0, 3,
                                                                       37587, 15687, 37617, 4824,
                                                                       4842, 16083, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38187, 0, 3,
                                                                       37617, 15705, 37647, 4842,
                                                                       4860, 16119, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38247, 0, 3,
                                                                       37647, 15723, 37677, 4860,
                                                                       4878, 16155, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38307, 0, 3,
                                                                       37677, 15741, 37707, 4878,
                                                                       4896, 16191, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38367, 0, 3,
                                                                       37707, 15759, 37737, 4896,
                                                                       4914, 16227, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38427, 0, 3,
                                                                       37737, 15777, 37767, 4914,
                                                                       4932, 16263, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38487, 0, 3,
                                                                       37767, 15795, 37797, 4932,
                                                                       4950, 16299, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38547, 0, 3,
                                                                       37797, 15813, 37827, 4950,
                                                                       4968, 16335, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 38607, 0, 3,
                                                                       37827, 15831, 37857, 4968,
                                                                       4986, 16371, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38667, 0, 3,
                                                                       37887, 15939, 37947, 5022,
                                                                       5052, 16527, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38767, 0, 3,
                                                                       37947, 15975, 38007, 5052,
                                                                       5082, 16587, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38867, 0, 3,
                                                                       38007, 16011, 38067, 5082,
                                                                       5112, 16647, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 38967, 0, 3,
                                                                       38067, 16047, 38127, 5112,
                                                                       5142, 16707, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39067, 0, 3,
                                                                       38127, 16083, 38187, 5142,
                                                                       5172, 16767, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39167, 0, 3,
                                                                       38187, 16119, 38247, 5172,
                                                                       5202, 16827, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39267, 0, 3,
                                                                       38247, 16155, 38307, 5202,
                                                                       5232, 16887, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39367, 0, 3,
                                                                       38307, 16191, 38367, 5232,
                                                                       5262, 16947, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39467, 0, 3,
                                                                       38367, 16227, 38427, 5262,
                                                                       5292, 17007, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39567, 0, 3,
                                                                       38427, 16263, 38487, 5292,
                                                                       5322, 17067, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39667, 0, 3,
                                                                       38487, 16299, 38547, 5322,
                                                                       5352, 17127, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 39767, 0, 3,
                                                                       38547, 16335, 38607, 5352,
                                                                       5382, 17187, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 39867, 0, 3,
                                                                       38667, 16527, 38767, 5442,
                                                                       5487, 17427, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40017, 0, 3,
                                                                       38767, 16587, 38867, 5487,
                                                                       5532, 17517, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40167, 0, 3,
                                                                       38867, 16647, 38967, 5532,
                                                                       5577, 17607, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40317, 0, 3,
                                                                       38967, 16707, 39067, 5577,
                                                                       5622, 17697, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40467, 0, 3,
                                                                       39067, 16767, 39167, 5622,
                                                                       5667, 17787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40617, 0, 3,
                                                                       39167, 16827, 39267, 5667,
                                                                       5712, 17877, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40767, 0, 3,
                                                                       39267, 16887, 39367, 5712,
                                                                       5757, 17967, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 40917, 0, 3,
                                                                       39367, 16947, 39467, 5757,
                                                                       5802, 18057, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41067, 0, 3,
                                                                       39467, 17007, 39567, 5802,
                                                                       5847, 18147, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41217, 0, 3,
                                                                       39567, 17067, 39667, 5847,
                                                                       5892, 18237, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 41367, 0, 3,
                                                                       39667, 17127, 39767, 5892,
                                                                       5937, 18327, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41517, 0, 3,
                                                                       39867, 17427, 40017, 6027,
                                                                       6090, 18669, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41727, 0, 3,
                                                                       40017, 17517, 40167, 6090,
                                                                       6153, 18795, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 41937, 0, 3,
                                                                       40167, 17607, 40317, 6153,
                                                                       6216, 18921, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42147, 0, 3,
                                                                       40317, 17697, 40467, 6216,
                                                                       6279, 19047, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42357, 0, 3,
                                                                       40467, 17787, 40617, 6279,
                                                                       6342, 19173, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42567, 0, 3,
                                                                       40617, 17877, 40767, 6342,
                                                                       6405, 19299, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42777, 0, 3,
                                                                       40767, 17967, 40917, 6405,
                                                                       6468, 19425, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 42987, 0, 3,
                                                                       40917, 18057, 41067, 6468,
                                                                       6531, 19551, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43197, 0, 3,
                                                                       41067, 18147, 41217, 6531,
                                                                       6594, 19677, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 43407, 0, 3,
                                                                       41217, 18237, 41367, 6594,
                                                                       6657, 19803, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43617, 0, 3,
                                                                       41517, 18669, 41727, 6783,
                                                                       6867, 20265, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 43897, 0, 3,
                                                                       41727, 18795, 41937, 6867,
                                                                       6951, 20433, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44177, 0, 3,
                                                                       41937, 18921, 42147, 6951,
                                                                       7035, 20601, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44457, 0, 3,
                                                                       42147, 19047, 42357, 7035,
                                                                       7119, 20769, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 44737, 0, 3,
                                                                       42357, 19173, 42567, 7119,
                                                                       7203, 20937, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45017, 0, 3,
                                                                       42567, 19299, 42777, 7203,
                                                                       7287, 21105, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45297, 0, 3,
                                                                       42777, 19425, 42987, 7287,
                                                                       7371, 21273, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45577, 0, 3,
                                                                       42987, 19551, 43197, 7371,
                                                                       7455, 21441, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 45857, 0, 3,
                                                                       43197, 19677, 43407, 7455,
                                                                       7539, 21609, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46137, 0, 3,
                                                                       43617, 20265, 43897, 7707,
                                                                       7815, 22209, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46497, 0, 3,
                                                                       43897, 20433, 44177, 7815,
                                                                       7923, 22425, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 46857, 0, 3,
                                                                       44177, 20601, 44457, 7923,
                                                                       8031, 22641, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47217, 0, 3,
                                                                       44457, 20769, 44737, 8031,
                                                                       8139, 22857, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47577, 0, 3,
                                                                       44737, 20937, 45017, 8139,
                                                                       8247, 23073, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 47937, 0, 3,
                                                                       45017, 21105, 45297, 8247,
                                                                       8355, 23289, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48297, 0, 3,
                                                                       45297, 21273, 45577, 8355,
                                                                       8463, 23505, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 48657, 0, 3,
                                                                       45577, 21441, 45857, 8463,
                                                                       8571, 23721, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49017, 0, 3,
                                                                       46137, 22209, 46497, 8787,
                                                                       8922, 24477, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49467, 0, 3,
                                                                       46497, 22425, 46857, 8922,
                                                                       9057, 24747, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 49917, 0, 3,
                                                                       46857, 22641, 47217, 9057,
                                                                       9192, 25017, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 50367, 0, 3,
                                                                       47217, 22857, 47577, 9192,
                                                                       9327, 25287, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 50817, 0, 3,
                                                                       47577, 23073, 47937, 9327,
                                                                       9462, 25557, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 51267, 0, 3,
                                                                       47937, 23289, 48297, 9462,
                                                                       9597, 25827, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 51717, 0, 3,
                                                                       48297, 23505, 48657, 9597,
                                                                       9732, 26097, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 52167, 0, 3,
                                                                       49017, 24477, 49467,
                                                                       10002, 10167, 27027,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 52717, 0, 3,
                                                                       49467, 24747, 49917,
                                                                       10167, 10332, 27357,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 53267, 0, 3,
                                                                       49917, 25017, 50367,
                                                                       10332, 10497, 27687,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 53817, 0, 3,
                                                                       50367, 25287, 50817,
                                                                       10497, 10662, 28017,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 54367, 0, 3,
                                                                       50817, 25557, 51267,
                                                                       10662, 10827, 28347,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 54917, 0, 3,
                                                                       51267, 25827, 51717,
                                                                       10827, 10992, 28677,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 55467, 0, 3,
                                                                       52167, 27027, 52717,
                                                                       11322, 11520, 29799,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 56127, 0, 3,
                                                                       52717, 27357, 53267,
                                                                       11520, 11718, 30195,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 56787, 0, 3,
                                                                       53267, 27687, 53817,
                                                                       11718, 11916, 30591,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 57447, 0, 3,
                                                                       53817, 28017, 54367,
                                                                       11916, 12114, 30987,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 58107, 0, 3,
                                                                       54367, 28347, 54917,
                                                                       12114, 12312, 31383,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 58767, 0, 3,
                                                                       55467, 29799, 56127,
                                                                       12708, 12942, 32715,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 59547, 0, 3,
                                                                       56127, 30195, 56787,
                                                                       12942, 13176, 33183,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 60327, 0, 3,
                                                                       56787, 30591, 57447,
                                                                       13176, 13410, 33651,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 61107, 0, 3,
                                                                       57447, 30987, 58107,
                                                                       13410, 13644, 34119,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 61887, 0, 3,
                                                                       58767, 32715, 59547,
                                                                       14112, 14385, 35679,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 62797, 0, 3,
                                                                       59547, 33183, 60327,
                                                                       14385, 14658, 36225,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 63707, 0, 3,
                                                                       60327, 33651, 61107,
                                                                       14658, 14931, 36771,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64617, 3, 15477,
                                                                       15483, 37317, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64632, 3, 15483,
                                                                       15489, 37327, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64647, 3, 15489,
                                                                       15495, 37337, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64662, 3, 15495,
                                                                       15501, 37347, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64677, 3, 15501,
                                                                       15507, 37357, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64692, 3, 15507,
                                                                       15513, 37367, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64707, 3, 15513,
                                                                       15519, 37377, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64722, 3, 15519,
                                                                       15525, 37387, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64737, 3, 15525,
                                                                       15531, 37397, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64752, 3, 15531,
                                                                       15537, 37407, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64767, 3, 15537,
                                                                       15543, 37417, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64782, 3, 15543,
                                                                       15549, 37427, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64797, 3, 15549,
                                                                       15555, 37437, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64812, 3, 15555,
                                                                       15561, 37447, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64827, 3, 15561,
                                                                       15567, 37457, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 64842, 0, 3,
                                                                       64617, 37317, 64632,
                                                                       15579, 15597, 37467,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 64887, 0, 3,
                                                                       64632, 37327, 64647,
                                                                       15597, 15615, 37497,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 64932, 0, 3,
                                                                       64647, 37337, 64662,
                                                                       15615, 15633, 37527,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 64977, 0, 3,
                                                                       64662, 37347, 64677,
                                                                       15633, 15651, 37557,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65022, 0, 3,
                                                                       64677, 37357, 64692,
                                                                       15651, 15669, 37587,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65067, 0, 3,
                                                                       64692, 37367, 64707,
                                                                       15669, 15687, 37617,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65112, 0, 3,
                                                                       64707, 37377, 64722,
                                                                       15687, 15705, 37647,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65157, 0, 3,
                                                                       64722, 37387, 64737,
                                                                       15705, 15723, 37677,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65202, 0, 3,
                                                                       64737, 37397, 64752,
                                                                       15723, 15741, 37707,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65247, 0, 3,
                                                                       64752, 37407, 64767,
                                                                       15741, 15759, 37737,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65292, 0, 3,
                                                                       64767, 37417, 64782,
                                                                       15759, 15777, 37767,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65337, 0, 3,
                                                                       64782, 37427, 64797,
                                                                       15777, 15795, 37797,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65382, 0, 3,
                                                                       64797, 37437, 64812,
                                                                       15795, 15813, 37827,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 65427, 0, 3,
                                                                       64812, 37447, 64827,
                                                                       15813, 15831, 37857,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65472, 0, 3,
                                                                       64842, 37467, 64887,
                                                                       15867, 15903, 37887,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65562, 0, 3,
                                                                       64887, 37497, 64932,
                                                                       15903, 15939, 37947,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65652, 0, 3,
                                                                       64932, 37527, 64977,
                                                                       15939, 15975, 38007,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65742, 0, 3,
                                                                       64977, 37557, 65022,
                                                                       15975, 16011, 38067,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65832, 0, 3,
                                                                       65022, 37587, 65067,
                                                                       16011, 16047, 38127,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 65922, 0, 3,
                                                                       65067, 37617, 65112,
                                                                       16047, 16083, 38187,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66012, 0, 3,
                                                                       65112, 37647, 65157,
                                                                       16083, 16119, 38247,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66102, 0, 3,
                                                                       65157, 37677, 65202,
                                                                       16119, 16155, 38307,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66192, 0, 3,
                                                                       65202, 37707, 65247,
                                                                       16155, 16191, 38367,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66282, 0, 3,
                                                                       65247, 37737, 65292,
                                                                       16191, 16227, 38427,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66372, 0, 3,
                                                                       65292, 37767, 65337,
                                                                       16227, 16263, 38487,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66462, 0, 3,
                                                                       65337, 37797, 65382,
                                                                       16263, 16299, 38547,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 66552, 0, 3,
                                                                       65382, 37827, 65427,
                                                                       16299, 16335, 38607,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 66642, 0, 3,
                                                                       65472, 37887, 65562,
                                                                       16407, 16467, 38667,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 66792, 0, 3,
                                                                       65562, 37947, 65652,
                                                                       16467, 16527, 38767,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 66942, 0, 3,
                                                                       65652, 38007, 65742,
                                                                       16527, 16587, 38867,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67092, 0, 3,
                                                                       65742, 38067, 65832,
                                                                       16587, 16647, 38967,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67242, 0, 3,
                                                                       65832, 38127, 65922,
                                                                       16647, 16707, 39067,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67392, 0, 3,
                                                                       65922, 38187, 66012,
                                                                       16707, 16767, 39167,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67542, 0, 3,
                                                                       66012, 38247, 66102,
                                                                       16767, 16827, 39267,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67692, 0, 3,
                                                                       66102, 38307, 66192,
                                                                       16827, 16887, 39367,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67842, 0, 3,
                                                                       66192, 38367, 66282,
                                                                       16887, 16947, 39467,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 67992, 0, 3,
                                                                       66282, 38427, 66372,
                                                                       16947, 17007, 39567,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 68142, 0, 3,
                                                                       66372, 38487, 66462,
                                                                       17007, 17067, 39667,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 68292, 0, 3,
                                                                       66462, 38547, 66552,
                                                                       17067, 17127, 39767,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 68442, 0, 3,
                                                                       66642, 38667, 66792,
                                                                       17247, 17337, 39867,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 68667, 0, 3,
                                                                       66792, 38767, 66942,
                                                                       17337, 17427, 40017,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 68892, 0, 3,
                                                                       66942, 38867, 67092,
                                                                       17427, 17517, 40167,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 69117, 0, 3,
                                                                       67092, 38967, 67242,
                                                                       17517, 17607, 40317,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 69342, 0, 3,
                                                                       67242, 39067, 67392,
                                                                       17607, 17697, 40467,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 69567, 0, 3,
                                                                       67392, 39167, 67542,
                                                                       17697, 17787, 40617,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 69792, 0, 3,
                                                                       67542, 39267, 67692,
                                                                       17787, 17877, 40767,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 70017, 0, 3,
                                                                       67692, 39367, 67842,
                                                                       17877, 17967, 40917,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 70242, 0, 3,
                                                                       67842, 39467, 67992,
                                                                       17967, 18057, 41067,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 70467, 0, 3,
                                                                       67992, 39567, 68142,
                                                                       18057, 18147, 41217,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 70692, 0, 3,
                                                                       68142, 39667, 68292,
                                                                       18147, 18237, 41367,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 70917, 0, 3,
                                                                       68442, 39867, 68667,
                                                                       18417, 18543, 41517,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 71232, 0, 3,
                                                                       68667, 40017, 68892,
                                                                       18543, 18669, 41727,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 71547, 0, 3,
                                                                       68892, 40167, 69117,
                                                                       18669, 18795, 41937,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 71862, 0, 3,
                                                                       69117, 40317, 69342,
                                                                       18795, 18921, 42147,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 72177, 0, 3,
                                                                       69342, 40467, 69567,
                                                                       18921, 19047, 42357,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 72492, 0, 3,
                                                                       69567, 40617, 69792,
                                                                       19047, 19173, 42567,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 72807, 0, 3,
                                                                       69792, 40767, 70017,
                                                                       19173, 19299, 42777,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 73122, 0, 3,
                                                                       70017, 40917, 70242,
                                                                       19299, 19425, 42987,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 73437, 0, 3,
                                                                       70242, 41067, 70467,
                                                                       19425, 19551, 43197,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 73752, 0, 3,
                                                                       70467, 41217, 70692,
                                                                       19551, 19677, 43407,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 74067, 0, 3,
                                                                       70917, 41517, 71232,
                                                                       19929, 20097, 43617,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 74487, 0, 3,
                                                                       71232, 41727, 71547,
                                                                       20097, 20265, 43897,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 74907, 0, 3,
                                                                       71547, 41937, 71862,
                                                                       20265, 20433, 44177,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 75327, 0, 3,
                                                                       71862, 42147, 72177,
                                                                       20433, 20601, 44457,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 75747, 0, 3,
                                                                       72177, 42357, 72492,
                                                                       20601, 20769, 44737,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 76167, 0, 3,
                                                                       72492, 42567, 72807,
                                                                       20769, 20937, 45017,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 76587, 0, 3,
                                                                       72807, 42777, 73122,
                                                                       20937, 21105, 45297,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 77007, 0, 3,
                                                                       73122, 42987, 73437,
                                                                       21105, 21273, 45577,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 77427, 0, 3,
                                                                       73437, 43197, 73752,
                                                                       21273, 21441, 45857,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 77847, 0, 3,
                                                                       74067, 43617, 74487,
                                                                       21777, 21993, 46137,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 78387, 0, 3,
                                                                       74487, 43897, 74907,
                                                                       21993, 22209, 46497,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 78927, 0, 3,
                                                                       74907, 44177, 75327,
                                                                       22209, 22425, 46857,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 79467, 0, 3,
                                                                       75327, 44457, 75747,
                                                                       22425, 22641, 47217,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 80007, 0, 3,
                                                                       75747, 44737, 76167,
                                                                       22641, 22857, 47577,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 80547, 0, 3,
                                                                       76167, 45017, 76587,
                                                                       22857, 23073, 47937,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 81087, 0, 3,
                                                                       76587, 45297, 77007,
                                                                       23073, 23289, 48297,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 81627, 0, 3,
                                                                       77007, 45577, 77427,
                                                                       23289, 23505, 48657,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 82167, 0, 3,
                                                                       77847, 46137, 78387,
                                                                       23937, 24207, 49017,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 82842, 0, 3,
                                                                       78387, 46497, 78927,
                                                                       24207, 24477, 49467,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 83517, 0, 3,
                                                                       78927, 46857, 79467,
                                                                       24477, 24747, 49917,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 84192, 0, 3,
                                                                       79467, 47217, 80007,
                                                                       24747, 25017, 50367,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 84867, 0, 3,
                                                                       80007, 47577, 80547,
                                                                       25017, 25287, 50817,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 85542, 0, 3,
                                                                       80547, 47937, 81087,
                                                                       25287, 25557, 51267,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 86217, 0, 3,
                                                                       81087, 48297, 81627,
                                                                       25557, 25827, 51717,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 86892, 0, 3,
                                                                       82167, 49017, 82842,
                                                                       26367, 26697, 52167,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 87717, 0, 3,
                                                                       82842, 49467, 83517,
                                                                       26697, 27027, 52717,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 88542, 0, 3,
                                                                       83517, 49917, 84192,
                                                                       27027, 27357, 53267,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 89367, 0, 3,
                                                                       84192, 50367, 84867,
                                                                       27357, 27687, 53817,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 90192, 0, 3,
                                                                       84867, 50817, 85542,
                                                                       27687, 28017, 54367,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 91017, 0, 3,
                                                                       85542, 51267, 86217,
                                                                       28017, 28347, 54917,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 91842, 0, 3,
                                                                       86892, 52167, 87717,
                                                                       29007, 29403, 55467,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 92832, 0, 3,
                                                                       87717, 52717, 88542,
                                                                       29403, 29799, 56127,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 93822, 0, 3,
                                                                       88542, 53267, 89367,
                                                                       29799, 30195, 56787,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 94812, 0, 3,
                                                                       89367, 53817, 90192,
                                                                       30195, 30591, 57447,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 95802, 0, 3,
                                                                       90192, 54367, 91017,
                                                                       30591, 30987, 58107,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 96792, 0, 3,
                                                                       91842, 55467, 92832,
                                                                       31779, 32247, 58767,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 97962, 0, 3,
                                                                       92832, 56127, 93822,
                                                                       32247, 32715, 59547,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 99132, 0, 3,
                                                                       93822, 56787, 94812,
                                                                       32715, 33183, 60327,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 100302, 0, 3,
                                                                       94812, 57447, 95802,
                                                                       33183, 33651, 61107,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 101472, 0, 3,
                                                                       96792, 58767, 97962,
                                                                       34587, 35133, 61887,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 102837, 0, 3,
                                                                       97962, 59547, 99132,
                                                                       35133, 35679, 62797,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 104202, 0, 3,
                                                                       99132, 60327, 100302,
                                                                       35679, 36225, 63707,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105567, 3, 37317,
                                                                       37327, 64647, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105588, 3, 37327,
                                                                       37337, 64662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105609, 3, 37337,
                                                                       37347, 64677, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105630, 3, 37347,
                                                                       37357, 64692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105651, 3, 37357,
                                                                       37367, 64707, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105672, 3, 37367,
                                                                       37377, 64722, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105693, 3, 37377,
                                                                       37387, 64737, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105714, 3, 37387,
                                                                       37397, 64752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105735, 3, 37397,
                                                                       37407, 64767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105756, 3, 37407,
                                                                       37417, 64782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105777, 3, 37417,
                                                                       37427, 64797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105798, 3, 37427,
                                                                       37437, 64812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 105819, 3, 37437,
                                                                       37447, 64827, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 105840, 0, 3,
                                                                       105567, 64647, 105588,
                                                                       37467, 37497, 64932,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 105903, 0, 3,
                                                                       105588, 64662, 105609,
                                                                       37497, 37527, 64977,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 105966, 0, 3,
                                                                       105609, 64677, 105630,
                                                                       37527, 37557, 65022,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106029, 0, 3,
                                                                       105630, 64692, 105651,
                                                                       37557, 37587, 65067,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106092, 0, 3,
                                                                       105651, 64707, 105672,
                                                                       37587, 37617, 65112,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106155, 0, 3,
                                                                       105672, 64722, 105693,
                                                                       37617, 37647, 65157,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106218, 0, 3,
                                                                       105693, 64737, 105714,
                                                                       37647, 37677, 65202,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106281, 0, 3,
                                                                       105714, 64752, 105735,
                                                                       37677, 37707, 65247,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106344, 0, 3,
                                                                       105735, 64767, 105756,
                                                                       37707, 37737, 65292,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106407, 0, 3,
                                                                       105756, 64782, 105777,
                                                                       37737, 37767, 65337,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106470, 0, 3,
                                                                       105777, 64797, 105798,
                                                                       37767, 37797, 65382,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 106533, 0, 3,
                                                                       105798, 64812, 105819,
                                                                       37797, 37827, 65427,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106596, 0, 3,
                                                                       105840, 64932, 105903,
                                                                       37887, 37947, 65652,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106722, 0, 3,
                                                                       105903, 64977, 105966,
                                                                       37947, 38007, 65742,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106848, 0, 3,
                                                                       105966, 65022, 106029,
                                                                       38007, 38067, 65832,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 106974, 0, 3,
                                                                       106029, 65067, 106092,
                                                                       38067, 38127, 65922,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107100, 0, 3,
                                                                       106092, 65112, 106155,
                                                                       38127, 38187, 66012,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107226, 0, 3,
                                                                       106155, 65157, 106218,
                                                                       38187, 38247, 66102,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107352, 0, 3,
                                                                       106218, 65202, 106281,
                                                                       38247, 38307, 66192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107478, 0, 3,
                                                                       106281, 65247, 106344,
                                                                       38307, 38367, 66282,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107604, 0, 3,
                                                                       106344, 65292, 106407,
                                                                       38367, 38427, 66372,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107730, 0, 3,
                                                                       106407, 65337, 106470,
                                                                       38427, 38487, 66462,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 107856, 0, 3,
                                                                       106470, 65382, 106533,
                                                                       38487, 38547, 66552,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 107982, 0, 3,
                                                                       106596, 65652, 106722,
                                                                       38667, 38767, 66942,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108192, 0, 3,
                                                                       106722, 65742, 106848,
                                                                       38767, 38867, 67092,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108402, 0, 3,
                                                                       106848, 65832, 106974,
                                                                       38867, 38967, 67242,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108612, 0, 3,
                                                                       106974, 65922, 107100,
                                                                       38967, 39067, 67392,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 108822, 0, 3,
                                                                       107100, 66012, 107226,
                                                                       39067, 39167, 67542,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109032, 0, 3,
                                                                       107226, 66102, 107352,
                                                                       39167, 39267, 67692,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109242, 0, 3,
                                                                       107352, 66192, 107478,
                                                                       39267, 39367, 67842,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109452, 0, 3,
                                                                       107478, 66282, 107604,
                                                                       39367, 39467, 67992,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109662, 0, 3,
                                                                       107604, 66372, 107730,
                                                                       39467, 39567, 68142,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 109872, 0, 3,
                                                                       107730, 66462, 107856,
                                                                       39567, 39667, 68292,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110082, 0, 3,
                                                                       107982, 66942, 108192,
                                                                       39867, 40017, 68892,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110397, 0, 3,
                                                                       108192, 67092, 108402,
                                                                       40017, 40167, 69117,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 110712, 0, 3,
                                                                       108402, 67242, 108612,
                                                                       40167, 40317, 69342,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111027, 0, 3,
                                                                       108612, 67392, 108822,
                                                                       40317, 40467, 69567,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111342, 0, 3,
                                                                       108822, 67542, 109032,
                                                                       40467, 40617, 69792,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111657, 0, 3,
                                                                       109032, 67692, 109242,
                                                                       40617, 40767, 70017,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 111972, 0, 3,
                                                                       109242, 67842, 109452,
                                                                       40767, 40917, 70242,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 112287, 0, 3,
                                                                       109452, 67992, 109662,
                                                                       40917, 41067, 70467,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 112602, 0, 3,
                                                                       109662, 68142, 109872,
                                                                       41067, 41217, 70692,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 112917, 0, 3,
                                                                       110082, 68892, 110397,
                                                                       41517, 41727, 71547,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 113358, 0, 3,
                                                                       110397, 69117, 110712,
                                                                       41727, 41937, 71862,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 113799, 0, 3,
                                                                       110712, 69342, 111027,
                                                                       41937, 42147, 72177,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 114240, 0, 3,
                                                                       111027, 69567, 111342,
                                                                       42147, 42357, 72492,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 114681, 0, 3,
                                                                       111342, 69792, 111657,
                                                                       42357, 42567, 72807,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 115122, 0, 3,
                                                                       111657, 70017, 111972,
                                                                       42567, 42777, 73122,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 115563, 0, 3,
                                                                       111972, 70242, 112287,
                                                                       42777, 42987, 73437,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 116004, 0, 3,
                                                                       112287, 70467, 112602,
                                                                       42987, 43197, 73752,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 116445, 0, 3,
                                                                       112917, 71547, 113358,
                                                                       43617, 43897, 74907,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 117033, 0, 3,
                                                                       113358, 71862, 113799,
                                                                       43897, 44177, 75327,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 117621, 0, 3,
                                                                       113799, 72177, 114240,
                                                                       44177, 44457, 75747,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 118209, 0, 3,
                                                                       114240, 72492, 114681,
                                                                       44457, 44737, 76167,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 118797, 0, 3,
                                                                       114681, 72807, 115122,
                                                                       44737, 45017, 76587,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 119385, 0, 3,
                                                                       115122, 73122, 115563,
                                                                       45017, 45297, 77007,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 119973, 0, 3,
                                                                       115563, 73437, 116004,
                                                                       45297, 45577, 77427,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 120561, 0, 3,
                                                                       116445, 74907, 117033,
                                                                       46137, 46497, 78927,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 121317, 0, 3,
                                                                       117033, 75327, 117621,
                                                                       46497, 46857, 79467,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 122073, 0, 3,
                                                                       117621, 75747, 118209,
                                                                       46857, 47217, 80007,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 122829, 0, 3,
                                                                       118209, 76167, 118797,
                                                                       47217, 47577, 80547,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 123585, 0, 3,
                                                                       118797, 76587, 119385,
                                                                       47577, 47937, 81087,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 124341, 0, 3,
                                                                       119385, 77007, 119973,
                                                                       47937, 48297, 81627,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 125097, 0, 3,
                                                                       120561, 78927, 121317,
                                                                       49017, 49467, 83517,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 126042, 0, 3,
                                                                       121317, 79467, 122073,
                                                                       49467, 49917, 84192,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 126987, 0, 3,
                                                                       122073, 80007, 122829,
                                                                       49917, 50367, 84867,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 127932, 0, 3,
                                                                       122829, 80547, 123585,
                                                                       50367, 50817, 85542,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 128877, 0, 3,
                                                                       123585, 81087, 124341,
                                                                       50817, 51267, 86217,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 129822, 0, 3,
                                                                       125097, 83517, 126042,
                                                                       52167, 52717, 88542,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 130977, 0, 3,
                                                                       126042, 84192, 126987,
                                                                       52717, 53267, 89367,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 132132, 0, 3,
                                                                       126987, 84867, 127932,
                                                                       53267, 53817, 90192,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 133287, 0, 3,
                                                                       127932, 85542, 128877,
                                                                       53817, 54367, 91017,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 134442, 0, 3,
                                                                       129822, 88542, 130977,
                                                                       55467, 56127, 93822,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 135828, 0, 3,
                                                                       130977, 89367, 132132,
                                                                       56127, 56787, 94812,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 137214, 0, 3,
                                                                       132132, 90192, 133287,
                                                                       56787, 57447, 95802,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 138600, 0, 3,
                                                                       134442, 93822, 135828,
                                                                       58767, 59547, 99132,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 140238, 0, 3,
                                                                       135828, 94812, 137214,
                                                                       59547, 60327, 100302,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 141876, 0, 3,
                                                                       138600, 99132, 140238,
                                                                       61887, 62797, 104202,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143787, 3, 64617,
                                                                       64632, 105567, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143815, 3, 64632,
                                                                       64647, 105588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143843, 3, 64647,
                                                                       64662, 105609, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143871, 3, 64662,
                                                                       64677, 105630, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143899, 3, 64677,
                                                                       64692, 105651, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143927, 3, 64692,
                                                                       64707, 105672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143955, 3, 64707,
                                                                       64722, 105693, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 143983, 3, 64722,
                                                                       64737, 105714, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144011, 3, 64737,
                                                                       64752, 105735, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144039, 3, 64752,
                                                                       64767, 105756, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144067, 3, 64767,
                                                                       64782, 105777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144095, 3, 64782,
                                                                       64797, 105798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 144123, 3, 64797,
                                                                       64812, 105819, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144151, 0, 3,
                                                                       143787, 105567, 143815,
                                                                       64842, 64887, 105840,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144235, 0, 3,
                                                                       143815, 105588, 143843,
                                                                       64887, 64932, 105903,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144319, 0, 3,
                                                                       143843, 105609, 143871,
                                                                       64932, 64977, 105966,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144403, 0, 3,
                                                                       143871, 105630, 143899,
                                                                       64977, 65022, 106029,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144487, 0, 3,
                                                                       143899, 105651, 143927,
                                                                       65022, 65067, 106092,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144571, 0, 3,
                                                                       143927, 105672, 143955,
                                                                       65067, 65112, 106155,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144655, 0, 3,
                                                                       143955, 105693, 143983,
                                                                       65112, 65157, 106218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144739, 0, 3,
                                                                       143983, 105714, 144011,
                                                                       65157, 65202, 106281,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144823, 0, 3,
                                                                       144011, 105735, 144039,
                                                                       65202, 65247, 106344,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144907, 0, 3,
                                                                       144039, 105756, 144067,
                                                                       65247, 65292, 106407,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 144991, 0, 3,
                                                                       144067, 105777, 144095,
                                                                       65292, 65337, 106470,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 145075, 0, 3,
                                                                       144095, 105798, 144123,
                                                                       65337, 65382, 106533,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145159, 0, 3,
                                                                       144151, 105840, 144235,
                                                                       65472, 65562, 106596,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145327, 0, 3,
                                                                       144235, 105903, 144319,
                                                                       65562, 65652, 106722,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145495, 0, 3,
                                                                       144319, 105966, 144403,
                                                                       65652, 65742, 106848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145663, 0, 3,
                                                                       144403, 106029, 144487,
                                                                       65742, 65832, 106974,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145831, 0, 3,
                                                                       144487, 106092, 144571,
                                                                       65832, 65922, 107100,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 145999, 0, 3,
                                                                       144571, 106155, 144655,
                                                                       65922, 66012, 107226,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146167, 0, 3,
                                                                       144655, 106218, 144739,
                                                                       66012, 66102, 107352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146335, 0, 3,
                                                                       144739, 106281, 144823,
                                                                       66102, 66192, 107478,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146503, 0, 3,
                                                                       144823, 106344, 144907,
                                                                       66192, 66282, 107604,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146671, 0, 3,
                                                                       144907, 106407, 144991,
                                                                       66282, 66372, 107730,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 146839, 0, 3,
                                                                       144991, 106470, 145075,
                                                                       66372, 66462, 107856,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 147007, 0, 3,
                                                                       145159, 106596, 145327,
                                                                       66642, 66792, 107982,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 147287, 0, 3,
                                                                       145327, 106722, 145495,
                                                                       66792, 66942, 108192,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 147567, 0, 3,
                                                                       145495, 106848, 145663,
                                                                       66942, 67092, 108402,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 147847, 0, 3,
                                                                       145663, 106974, 145831,
                                                                       67092, 67242, 108612,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 148127, 0, 3,
                                                                       145831, 107100, 145999,
                                                                       67242, 67392, 108822,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 148407, 0, 3,
                                                                       145999, 107226, 146167,
                                                                       67392, 67542, 109032,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 148687, 0, 3,
                                                                       146167, 107352, 146335,
                                                                       67542, 67692, 109242,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 148967, 0, 3,
                                                                       146335, 107478, 146503,
                                                                       67692, 67842, 109452,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 149247, 0, 3,
                                                                       146503, 107604, 146671,
                                                                       67842, 67992, 109662,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 149527, 0, 3,
                                                                       146671, 107730, 146839,
                                                                       67992, 68142, 109872,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 149807, 0, 3,
                                                                       147007, 107982, 147287,
                                                                       68442, 68667, 110082,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 150227, 0, 3,
                                                                       147287, 108192, 147567,
                                                                       68667, 68892, 110397,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 150647, 0, 3,
                                                                       147567, 108402, 147847,
                                                                       68892, 69117, 110712,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 151067, 0, 3,
                                                                       147847, 108612, 148127,
                                                                       69117, 69342, 111027,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 151487, 0, 3,
                                                                       148127, 108822, 148407,
                                                                       69342, 69567, 111342,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 151907, 0, 3,
                                                                       148407, 109032, 148687,
                                                                       69567, 69792, 111657,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 152327, 0, 3,
                                                                       148687, 109242, 148967,
                                                                       69792, 70017, 111972,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 152747, 0, 3,
                                                                       148967, 109452, 149247,
                                                                       70017, 70242, 112287,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 153167, 0, 3,
                                                                       149247, 109662, 149527,
                                                                       70242, 70467, 112602,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 153587, 0, 3,
                                                                       149807, 110082, 150227,
                                                                       70917, 71232, 112917,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 154175, 0, 3,
                                                                       150227, 110397, 150647,
                                                                       71232, 71547, 113358,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 154763, 0, 3,
                                                                       150647, 110712, 151067,
                                                                       71547, 71862, 113799,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 155351, 0, 3,
                                                                       151067, 111027, 151487,
                                                                       71862, 72177, 114240,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 155939, 0, 3,
                                                                       151487, 111342, 151907,
                                                                       72177, 72492, 114681,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 156527, 0, 3,
                                                                       151907, 111657, 152327,
                                                                       72492, 72807, 115122,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 157115, 0, 3,
                                                                       152327, 111972, 152747,
                                                                       72807, 73122, 115563,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 157703, 0, 3,
                                                                       152747, 112287, 153167,
                                                                       73122, 73437, 116004,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 158291, 0, 3,
                                                                       153587, 112917, 154175,
                                                                       74067, 74487, 116445,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 159075, 0, 3,
                                                                       154175, 113358, 154763,
                                                                       74487, 74907, 117033,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 159859, 0, 3,
                                                                       154763, 113799, 155351,
                                                                       74907, 75327, 117621,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 160643, 0, 3,
                                                                       155351, 114240, 155939,
                                                                       75327, 75747, 118209,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 161427, 0, 3,
                                                                       155939, 114681, 156527,
                                                                       75747, 76167, 118797,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 162211, 0, 3,
                                                                       156527, 115122, 157115,
                                                                       76167, 76587, 119385,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 162995, 0, 3,
                                                                       157115, 115563, 157703,
                                                                       76587, 77007, 119973,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 163779, 0, 3,
                                                                       158291, 116445, 159075,
                                                                       77847, 78387, 120561,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 164787, 0, 3,
                                                                       159075, 117033, 159859,
                                                                       78387, 78927, 121317,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 165795, 0, 3,
                                                                       159859, 117621, 160643,
                                                                       78927, 79467, 122073,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 166803, 0, 3,
                                                                       160643, 118209, 161427,
                                                                       79467, 80007, 122829,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 167811, 0, 3,
                                                                       161427, 118797, 162211,
                                                                       80007, 80547, 123585,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 168819, 0, 3,
                                                                       162211, 119385, 162995,
                                                                       80547, 81087, 124341,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 169827, 0, 3,
                                                                       163779, 120561, 164787,
                                                                       82167, 82842, 125097,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 171087, 0, 3,
                                                                       164787, 121317, 165795,
                                                                       82842, 83517, 126042,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 172347, 0, 3,
                                                                       165795, 122073, 166803,
                                                                       83517, 84192, 126987,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 173607, 0, 3,
                                                                       166803, 122829, 167811,
                                                                       84192, 84867, 127932,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 174867, 0, 3,
                                                                       167811, 123585, 168819,
                                                                       84867, 85542, 128877,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 176127, 0, 3,
                                                                       169827, 125097, 171087,
                                                                       86892, 87717, 129822,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 177667, 0, 3,
                                                                       171087, 126042, 172347,
                                                                       87717, 88542, 130977,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 179207, 0, 3,
                                                                       172347, 126987, 173607,
                                                                       88542, 89367, 132132,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 180747, 0, 3,
                                                                       173607, 127932, 174867,
                                                                       89367, 90192, 133287,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 182287, 0, 3,
                                                                       176127, 129822, 177667,
                                                                       91842, 92832, 134442,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 184135, 0, 3,
                                                                       177667, 130977, 179207,
                                                                       92832, 93822, 135828,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 185983, 0, 3,
                                                                       179207, 132132, 180747,
                                                                       93822, 94812, 137214,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 187831, 0, 3,
                                                                       182287, 134442, 184135,
                                                                       96792, 97962, 138600,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 190015, 0, 3,
                                                                       184135, 135828, 185983,
                                                                       97962, 99132, 140238,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 192199, 0, 3,
                                                                       187831, 138600, 190015,
                                                                       101472, 102837, 141876,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 194747, 158291, 784, ncols);

                    simdfunc::contract_primitives(buffer, 195895, 163779, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 197371, 169827, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 199216, 176127, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 201471, 182287, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 204177, 187831, 2184, ncols);

                    simdfunc::contract_primitives(buffer, 207375, 192199, 2548, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 195531, 194747, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 196903, 195895, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 198631, 197371, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 200756, 199216, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 203319, 201471, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 206361, 204177, 78, 1, nmax);

        simdtrf::transform_i_inner(buffer, 209923, 207375, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 211106, 195531, 196903, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 212198, 196903, 198631, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 213602, 198631, 200756, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 215357, 200756, 203319, 13,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 217502, 203319, 206361, 13,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 220076, 206361, 209923, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 223118, 211106, 212198, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 225302, 212198, 213602, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 228110, 213602, 215357, 13,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 231620, 215357, 217502, 13,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 235910, 217502, 220076, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 241058, 223118, 225302, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 244698, 225302, 228110, 13,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 249378, 228110, 231620, 13,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 255228, 231620, 235910, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 262378, 241058, 244698, 13,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 267838, 244698, 249378, 13,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 274858, 249378, 255228, 13,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 283633, 262378, 267838, 13,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 291277, 267838, 274858, 13,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 301105, 283633, 291277, 13,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 311297, 301105, 28, 13, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 311297, 169, nmax);
    }

    for (size_t m = 0; m < 2197; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
