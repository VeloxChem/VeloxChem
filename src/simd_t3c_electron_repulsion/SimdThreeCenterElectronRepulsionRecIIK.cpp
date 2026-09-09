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


#include "SimdThreeCenterElectronRepulsionRecIIK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
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
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iik_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 432621, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2535 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 432621, 291207, 18984, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17, 18, 19}, ncols, fj, mu, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4557, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4560, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4563, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4566, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4569, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4572, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4575, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4578, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4581, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4584, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4587, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4590, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4593, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4596, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4599, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4602, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4605, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4608, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4611, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4614, 3, 7, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4623, 3, 8, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4632, 3, 9, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4641, 3, 10, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4650, 3, 11, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4659, 3, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4668, 3, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4677, 3, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4686, 3, 15, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4695, 3, 16, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4704, 3, 17, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4713, 3, 18, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4722, 3, 19, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4731, 3, 20, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4740, 3, 21, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4749, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4758, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 4767, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4776, 3, 26, 80,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4794, 3, 29, 86,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4812, 3, 32, 92,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4830, 3, 35, 98,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4848, 3, 38, 104,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4866, 3, 41, 110,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4884, 3, 44, 116,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4902, 3, 47, 122,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4920, 3, 50, 128,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4938, 3, 53, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4956, 3, 56, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4974, 3, 59, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 4992, 3, 62, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5010, 3, 65, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5028, 3, 68, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5046, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5064, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5082, 3, 80, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5112, 3, 86, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5142, 3, 92, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5172, 3, 98, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5202, 3, 104, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5232, 3, 110, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5262, 3, 116, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5292, 3, 122, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5322, 3, 128, 262,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5352, 3, 134, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5382, 3, 140, 282,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5412, 3, 146, 292,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5442, 3, 152, 302,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5472, 3, 158, 312,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5502, 3, 164, 322,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5532, 3, 170, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5562, 3, 182, 342,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5607, 3, 192, 357,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5652, 3, 202, 372,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5697, 3, 212, 387,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5742, 3, 222, 402,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5787, 3, 232, 417,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5832, 3, 242, 432,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5877, 3, 252, 447,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5922, 3, 262, 462,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5967, 3, 272, 477,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6012, 3, 282, 492,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6057, 3, 292, 507,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6102, 3, 302, 522,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6147, 3, 312, 537,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6192, 3, 322, 552,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6237, 3, 342, 567,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6300, 3, 357, 588,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6363, 3, 372, 609,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6426, 3, 387, 630,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6489, 3, 402, 651,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6552, 3, 417, 672,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6615, 3, 432, 693,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6678, 3, 447, 714,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6741, 3, 462, 735,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6804, 3, 477, 756,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6867, 3, 492, 777,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6930, 3, 507, 798,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6993, 3, 522, 819,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7056, 3, 537, 840,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7119, 3, 567, 861,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7203, 3, 588, 889,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7287, 3, 609, 917,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7371, 3, 630, 945,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7455, 3, 651, 973,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7539, 3, 672,
                                                                       1001, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7623, 3, 693,
                                                                       1029, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7707, 3, 714,
                                                                       1057, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7791, 3, 735,
                                                                       1085, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7875, 3, 756,
                                                                       1113, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7959, 3, 777,
                                                                       1141, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8043, 3, 798,
                                                                       1169, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8127, 3, 819,
                                                                       1197, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8211, 3, 861,
                                                                       1225, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8319, 3, 889,
                                                                       1261, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8427, 3, 917,
                                                                       1297, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8535, 3, 945,
                                                                       1333, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8643, 3, 973,
                                                                       1369, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8751, 3, 1001,
                                                                       1405, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8859, 3, 1029,
                                                                       1441, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8967, 3, 1057,
                                                                       1477, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9075, 3, 1085,
                                                                       1513, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9183, 3, 1113,
                                                                       1549, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9291, 3, 1141,
                                                                       1585, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9399, 3, 1169,
                                                                       1621, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9507, 3, 1225,
                                                                       1657, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9642, 3, 1261,
                                                                       1702, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9777, 3, 1297,
                                                                       1747, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9912, 3, 1333,
                                                                       1792, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10047, 3, 1369,
                                                                       1837, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10182, 3, 1405,
                                                                       1882, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10317, 3, 1441,
                                                                       1927, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10452, 3, 1477,
                                                                       1972, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10587, 3, 1513,
                                                                       2017, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10722, 3, 1549,
                                                                       2062, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10857, 3, 1585,
                                                                       2107, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 10992, 3, 1657,
                                                                       2152, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11157, 3, 1702,
                                                                       2207, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11322, 3, 1747,
                                                                       2262, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11487, 3, 1792,
                                                                       2317, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11652, 3, 1837,
                                                                       2372, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11817, 3, 1882,
                                                                       2427, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11982, 3, 1927,
                                                                       2482, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12147, 3, 1972,
                                                                       2537, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12312, 3, 2017,
                                                                       2592, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12477, 3, 2062,
                                                                       2647, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12642, 3, 2152,
                                                                       2702, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 12840, 3, 2207,
                                                                       2768, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13038, 3, 2262,
                                                                       2834, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13236, 3, 2317,
                                                                       2900, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13434, 3, 2372,
                                                                       2966, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13632, 3, 2427,
                                                                       3032, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13830, 3, 2482,
                                                                       3098, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14028, 3, 2537,
                                                                       3164, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14226, 3, 2592,
                                                                       3230, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 14424, 3, 2702,
                                                                       3296, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 14658, 3, 2768,
                                                                       3374, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 14892, 3, 2834,
                                                                       3452, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15126, 3, 2900,
                                                                       3530, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15360, 3, 2966,
                                                                       3608, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15594, 3, 3032,
                                                                       3686, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15828, 3, 3098,
                                                                       3764, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16062, 3, 3164,
                                                                       3842, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 16296, 3, 3296,
                                                                       3920, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 16569, 3, 3374,
                                                                       4011, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 16842, 3, 3452,
                                                                       4102, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17115, 3, 3530,
                                                                       4193, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17388, 3, 3608,
                                                                       4284, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17661, 3, 3686,
                                                                       4375, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17934, 3, 3764,
                                                                       4466, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18207, 3, 7, 8,
                                                                       4563, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18213, 3, 8, 9,
                                                                       4566, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18219, 3, 9, 10,
                                                                       4569, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18225, 3, 10, 11,
                                                                       4572, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18231, 3, 11, 12,
                                                                       4575, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18237, 3, 12, 13,
                                                                       4578, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18243, 3, 13, 14,
                                                                       4581, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18249, 3, 14, 15,
                                                                       4584, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18255, 3, 15, 16,
                                                                       4587, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18261, 3, 16, 17,
                                                                       4590, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18267, 3, 17, 18,
                                                                       4593, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18273, 3, 18, 19,
                                                                       4596, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18279, 3, 19, 20,
                                                                       4599, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18285, 3, 20, 21,
                                                                       4602, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18291, 3, 21, 22,
                                                                       4605, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18297, 3, 22, 23,
                                                                       4608, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 18303, 3, 23, 24,
                                                                       4611, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18309, 0, 3,
                                                                       18207, 4563, 18213, 4632,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18327, 0, 3,
                                                                       18213, 4566, 18219, 4641,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18345, 0, 3,
                                                                       18219, 4569, 18225, 4650,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18363, 0, 3,
                                                                       18225, 4572, 18231, 4659,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18381, 0, 3,
                                                                       18231, 4575, 18237, 4668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18399, 0, 3,
                                                                       18237, 4578, 18243, 4677,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18417, 0, 3,
                                                                       18243, 4581, 18249, 4686,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18435, 0, 3,
                                                                       18249, 4584, 18255, 4695,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18453, 0, 3,
                                                                       18255, 4587, 18261, 4704,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18471, 0, 3,
                                                                       18261, 4590, 18267, 4713,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18489, 0, 3,
                                                                       18267, 4593, 18273, 4722,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18507, 0, 3,
                                                                       18273, 4596, 18279, 4731,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18525, 0, 3,
                                                                       18279, 4599, 18285, 4740,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18543, 0, 3,
                                                                       18285, 4602, 18291, 4749,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18561, 0, 3,
                                                                       18291, 4605, 18297, 4758,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 18579, 0, 3,
                                                                       18297, 4608, 18303, 4767,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18597, 0, 3,
                                                                       18309, 4632, 18327, 80,
                                                                       86, 4812, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18633, 0, 3,
                                                                       18327, 4641, 18345, 86,
                                                                       92, 4830, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18669, 0, 3,
                                                                       18345, 4650, 18363, 92,
                                                                       98, 4848, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18705, 0, 3,
                                                                       18363, 4659, 18381, 98,
                                                                       104, 4866, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18741, 0, 3,
                                                                       18381, 4668, 18399, 104,
                                                                       110, 4884, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18777, 0, 3,
                                                                       18399, 4677, 18417, 110,
                                                                       116, 4902, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18813, 0, 3,
                                                                       18417, 4686, 18435, 116,
                                                                       122, 4920, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18849, 0, 3,
                                                                       18435, 4695, 18453, 122,
                                                                       128, 4938, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18885, 0, 3,
                                                                       18453, 4704, 18471, 128,
                                                                       134, 4956, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18921, 0, 3,
                                                                       18471, 4713, 18489, 134,
                                                                       140, 4974, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18957, 0, 3,
                                                                       18489, 4722, 18507, 140,
                                                                       146, 4992, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 18993, 0, 3,
                                                                       18507, 4731, 18525, 146,
                                                                       152, 5010, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19029, 0, 3,
                                                                       18525, 4740, 18543, 152,
                                                                       158, 5028, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19065, 0, 3,
                                                                       18543, 4749, 18561, 158,
                                                                       164, 5046, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19101, 0, 3,
                                                                       18561, 4758, 18579, 164,
                                                                       170, 5064, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19137, 0, 3,
                                                                       18597, 4812, 18633, 182,
                                                                       192, 5142, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19197, 0, 3,
                                                                       18633, 4830, 18669, 192,
                                                                       202, 5172, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19257, 0, 3,
                                                                       18669, 4848, 18705, 202,
                                                                       212, 5202, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19317, 0, 3,
                                                                       18705, 4866, 18741, 212,
                                                                       222, 5232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19377, 0, 3,
                                                                       18741, 4884, 18777, 222,
                                                                       232, 5262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19437, 0, 3,
                                                                       18777, 4902, 18813, 232,
                                                                       242, 5292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19497, 0, 3,
                                                                       18813, 4920, 18849, 242,
                                                                       252, 5322, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19557, 0, 3,
                                                                       18849, 4938, 18885, 252,
                                                                       262, 5352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19617, 0, 3,
                                                                       18885, 4956, 18921, 262,
                                                                       272, 5382, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19677, 0, 3,
                                                                       18921, 4974, 18957, 272,
                                                                       282, 5412, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19737, 0, 3,
                                                                       18957, 4992, 18993, 282,
                                                                       292, 5442, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19797, 0, 3,
                                                                       18993, 5010, 19029, 292,
                                                                       302, 5472, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19857, 0, 3,
                                                                       19029, 5028, 19065, 302,
                                                                       312, 5502, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 19917, 0, 3,
                                                                       19065, 5046, 19101, 312,
                                                                       322, 5532, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 19977, 0, 3,
                                                                       19137, 5142, 19197, 342,
                                                                       357, 5652, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20067, 0, 3,
                                                                       19197, 5172, 19257, 357,
                                                                       372, 5697, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20157, 0, 3,
                                                                       19257, 5202, 19317, 372,
                                                                       387, 5742, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20247, 0, 3,
                                                                       19317, 5232, 19377, 387,
                                                                       402, 5787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20337, 0, 3,
                                                                       19377, 5262, 19437, 402,
                                                                       417, 5832, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20427, 0, 3,
                                                                       19437, 5292, 19497, 417,
                                                                       432, 5877, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20517, 0, 3,
                                                                       19497, 5322, 19557, 432,
                                                                       447, 5922, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20607, 0, 3,
                                                                       19557, 5352, 19617, 447,
                                                                       462, 5967, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20697, 0, 3,
                                                                       19617, 5382, 19677, 462,
                                                                       477, 6012, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20787, 0, 3,
                                                                       19677, 5412, 19737, 477,
                                                                       492, 6057, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20877, 0, 3,
                                                                       19737, 5442, 19797, 492,
                                                                       507, 6102, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 20967, 0, 3,
                                                                       19797, 5472, 19857, 507,
                                                                       522, 6147, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21057, 0, 3,
                                                                       19857, 5502, 19917, 522,
                                                                       537, 6192, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21147, 0, 3,
                                                                       19977, 5652, 20067, 567,
                                                                       588, 6363, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21273, 0, 3,
                                                                       20067, 5697, 20157, 588,
                                                                       609, 6426, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21399, 0, 3,
                                                                       20157, 5742, 20247, 609,
                                                                       630, 6489, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21525, 0, 3,
                                                                       20247, 5787, 20337, 630,
                                                                       651, 6552, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21651, 0, 3,
                                                                       20337, 5832, 20427, 651,
                                                                       672, 6615, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21777, 0, 3,
                                                                       20427, 5877, 20517, 672,
                                                                       693, 6678, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 21903, 0, 3,
                                                                       20517, 5922, 20607, 693,
                                                                       714, 6741, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22029, 0, 3,
                                                                       20607, 5967, 20697, 714,
                                                                       735, 6804, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22155, 0, 3,
                                                                       20697, 6012, 20787, 735,
                                                                       756, 6867, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22281, 0, 3,
                                                                       20787, 6057, 20877, 756,
                                                                       777, 6930, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22407, 0, 3,
                                                                       20877, 6102, 20967, 777,
                                                                       798, 6993, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22533, 0, 3,
                                                                       20967, 6147, 21057, 798,
                                                                       819, 7056, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22659, 0, 3,
                                                                       21147, 6363, 21273, 861,
                                                                       889, 7287, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22827, 0, 3,
                                                                       21273, 6426, 21399, 889,
                                                                       917, 7371, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 22995, 0, 3,
                                                                       21399, 6489, 21525, 917,
                                                                       945, 7455, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23163, 0, 3,
                                                                       21525, 6552, 21651, 945,
                                                                       973, 7539, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23331, 0, 3,
                                                                       21651, 6615, 21777, 973,
                                                                       1001, 7623, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23499, 0, 3,
                                                                       21777, 6678, 21903, 1001,
                                                                       1029, 7707, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23667, 0, 3,
                                                                       21903, 6741, 22029, 1029,
                                                                       1057, 7791, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 23835, 0, 3,
                                                                       22029, 6804, 22155, 1057,
                                                                       1085, 7875, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24003, 0, 3,
                                                                       22155, 6867, 22281, 1085,
                                                                       1113, 7959, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24171, 0, 3,
                                                                       22281, 6930, 22407, 1113,
                                                                       1141, 8043, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24339, 0, 3,
                                                                       22407, 6993, 22533, 1141,
                                                                       1169, 8127, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24507, 0, 3,
                                                                       22659, 7287, 22827, 1225,
                                                                       1261, 8427, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24723, 0, 3,
                                                                       22827, 7371, 22995, 1261,
                                                                       1297, 8535, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 24939, 0, 3,
                                                                       22995, 7455, 23163, 1297,
                                                                       1333, 8643, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25155, 0, 3,
                                                                       23163, 7539, 23331, 1333,
                                                                       1369, 8751, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25371, 0, 3,
                                                                       23331, 7623, 23499, 1369,
                                                                       1405, 8859, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25587, 0, 3,
                                                                       23499, 7707, 23667, 1405,
                                                                       1441, 8967, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 25803, 0, 3,
                                                                       23667, 7791, 23835, 1441,
                                                                       1477, 9075, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26019, 0, 3,
                                                                       23835, 7875, 24003, 1477,
                                                                       1513, 9183, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26235, 0, 3,
                                                                       24003, 7959, 24171, 1513,
                                                                       1549, 9291, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26451, 0, 3,
                                                                       24171, 8043, 24339, 1549,
                                                                       1585, 9399, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26667, 0, 3,
                                                                       24507, 8427, 24723, 1657,
                                                                       1702, 9777, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 26937, 0, 3,
                                                                       24723, 8535, 24939, 1702,
                                                                       1747, 9912, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27207, 0, 3,
                                                                       24939, 8643, 25155, 1747,
                                                                       1792, 10047, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27477, 0, 3,
                                                                       25155, 8751, 25371, 1792,
                                                                       1837, 10182, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 27747, 0, 3,
                                                                       25371, 8859, 25587, 1837,
                                                                       1882, 10317, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28017, 0, 3,
                                                                       25587, 8967, 25803, 1882,
                                                                       1927, 10452, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28287, 0, 3,
                                                                       25803, 9075, 26019, 1927,
                                                                       1972, 10587, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28557, 0, 3,
                                                                       26019, 9183, 26235, 1972,
                                                                       2017, 10722, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 28827, 0, 3,
                                                                       26235, 9291, 26451, 2017,
                                                                       2062, 10857, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29097, 0, 3,
                                                                       26667, 9777, 26937, 2152,
                                                                       2207, 11322, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29427, 0, 3,
                                                                       26937, 9912, 27207, 2207,
                                                                       2262, 11487, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 29757, 0, 3,
                                                                       27207, 10047, 27477, 2262,
                                                                       2317, 11652, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30087, 0, 3,
                                                                       27477, 10182, 27747, 2317,
                                                                       2372, 11817, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30417, 0, 3,
                                                                       27747, 10317, 28017, 2372,
                                                                       2427, 11982, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 30747, 0, 3,
                                                                       28017, 10452, 28287, 2427,
                                                                       2482, 12147, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31077, 0, 3,
                                                                       28287, 10587, 28557, 2482,
                                                                       2537, 12312, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31407, 0, 3,
                                                                       28557, 10722, 28827, 2537,
                                                                       2592, 12477, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 31737, 0, 3,
                                                                       29097, 11322, 29427, 2702,
                                                                       2768, 13038, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 32133, 0, 3,
                                                                       29427, 11487, 29757, 2768,
                                                                       2834, 13236, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 32529, 0, 3,
                                                                       29757, 11652, 30087, 2834,
                                                                       2900, 13434, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 32925, 0, 3,
                                                                       30087, 11817, 30417, 2900,
                                                                       2966, 13632, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 33321, 0, 3,
                                                                       30417, 11982, 30747, 2966,
                                                                       3032, 13830, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 33717, 0, 3,
                                                                       30747, 12147, 31077, 3032,
                                                                       3098, 14028, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 34113, 0, 3,
                                                                       31077, 12312, 31407, 3098,
                                                                       3164, 14226, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 34509, 0, 3,
                                                                       31737, 13038, 32133, 3296,
                                                                       3374, 14892, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 34977, 0, 3,
                                                                       32133, 13236, 32529, 3374,
                                                                       3452, 15126, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 35445, 0, 3,
                                                                       32529, 13434, 32925, 3452,
                                                                       3530, 15360, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 35913, 0, 3,
                                                                       32925, 13632, 33321, 3530,
                                                                       3608, 15594, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 36381, 0, 3,
                                                                       33321, 13830, 33717, 3608,
                                                                       3686, 15828, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 36849, 0, 3,
                                                                       33717, 14028, 34113, 3686,
                                                                       3764, 16062, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 37317, 0, 3,
                                                                       34509, 14892, 34977, 3920,
                                                                       4011, 16842, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 37863, 0, 3,
                                                                       34977, 15126, 35445, 4011,
                                                                       4102, 17115, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 38409, 0, 3,
                                                                       35445, 15360, 35913, 4102,
                                                                       4193, 17388, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 38955, 0, 3,
                                                                       35913, 15594, 36381, 4193,
                                                                       4284, 17661, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 39501, 0, 3,
                                                                       36381, 15828, 36849, 4284,
                                                                       4375, 17934, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40047, 3, 4557,
                                                                       4560, 18207, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40057, 3, 4560,
                                                                       4563, 18213, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40067, 3, 4563,
                                                                       4566, 18219, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40077, 3, 4566,
                                                                       4569, 18225, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40087, 3, 4569,
                                                                       4572, 18231, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40097, 3, 4572,
                                                                       4575, 18237, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40107, 3, 4575,
                                                                       4578, 18243, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40117, 3, 4578,
                                                                       4581, 18249, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40127, 3, 4581,
                                                                       4584, 18255, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40137, 3, 4584,
                                                                       4587, 18261, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40147, 3, 4587,
                                                                       4590, 18267, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40157, 3, 4590,
                                                                       4593, 18273, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40167, 3, 4593,
                                                                       4596, 18279, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40177, 3, 4596,
                                                                       4599, 18285, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40187, 3, 4599,
                                                                       4602, 18291, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40197, 3, 4602,
                                                                       4605, 18297, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 40207, 3, 4605,
                                                                       4608, 18303, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40217, 0, 3,
                                                                       40047, 18207, 40057, 4614,
                                                                       4623, 18309, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40247, 0, 3,
                                                                       40057, 18213, 40067, 4623,
                                                                       4632, 18327, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40277, 0, 3,
                                                                       40067, 18219, 40077, 4632,
                                                                       4641, 18345, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40307, 0, 3,
                                                                       40077, 18225, 40087, 4641,
                                                                       4650, 18363, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40337, 0, 3,
                                                                       40087, 18231, 40097, 4650,
                                                                       4659, 18381, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40367, 0, 3,
                                                                       40097, 18237, 40107, 4659,
                                                                       4668, 18399, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40397, 0, 3,
                                                                       40107, 18243, 40117, 4668,
                                                                       4677, 18417, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40427, 0, 3,
                                                                       40117, 18249, 40127, 4677,
                                                                       4686, 18435, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40457, 0, 3,
                                                                       40127, 18255, 40137, 4686,
                                                                       4695, 18453, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40487, 0, 3,
                                                                       40137, 18261, 40147, 4695,
                                                                       4704, 18471, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40517, 0, 3,
                                                                       40147, 18267, 40157, 4704,
                                                                       4713, 18489, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40547, 0, 3,
                                                                       40157, 18273, 40167, 4713,
                                                                       4722, 18507, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40577, 0, 3,
                                                                       40167, 18279, 40177, 4722,
                                                                       4731, 18525, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40607, 0, 3,
                                                                       40177, 18285, 40187, 4731,
                                                                       4740, 18543, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40637, 0, 3,
                                                                       40187, 18291, 40197, 4740,
                                                                       4749, 18561, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 40667, 0, 3,
                                                                       40197, 18297, 40207, 4749,
                                                                       4758, 18579, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40697, 0, 3,
                                                                       40217, 18309, 40247, 4776,
                                                                       4794, 18597, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40757, 0, 3,
                                                                       40247, 18327, 40277, 4794,
                                                                       4812, 18633, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40817, 0, 3,
                                                                       40277, 18345, 40307, 4812,
                                                                       4830, 18669, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40877, 0, 3,
                                                                       40307, 18363, 40337, 4830,
                                                                       4848, 18705, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40937, 0, 3,
                                                                       40337, 18381, 40367, 4848,
                                                                       4866, 18741, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 40997, 0, 3,
                                                                       40367, 18399, 40397, 4866,
                                                                       4884, 18777, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41057, 0, 3,
                                                                       40397, 18417, 40427, 4884,
                                                                       4902, 18813, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41117, 0, 3,
                                                                       40427, 18435, 40457, 4902,
                                                                       4920, 18849, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41177, 0, 3,
                                                                       40457, 18453, 40487, 4920,
                                                                       4938, 18885, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41237, 0, 3,
                                                                       40487, 18471, 40517, 4938,
                                                                       4956, 18921, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41297, 0, 3,
                                                                       40517, 18489, 40547, 4956,
                                                                       4974, 18957, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41357, 0, 3,
                                                                       40547, 18507, 40577, 4974,
                                                                       4992, 18993, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41417, 0, 3,
                                                                       40577, 18525, 40607, 4992,
                                                                       5010, 19029, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41477, 0, 3,
                                                                       40607, 18543, 40637, 5010,
                                                                       5028, 19065, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 41537, 0, 3,
                                                                       40637, 18561, 40667, 5028,
                                                                       5046, 19101, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41597, 0, 3,
                                                                       40697, 18597, 40757, 5082,
                                                                       5112, 19137, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41697, 0, 3,
                                                                       40757, 18633, 40817, 5112,
                                                                       5142, 19197, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41797, 0, 3,
                                                                       40817, 18669, 40877, 5142,
                                                                       5172, 19257, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41897, 0, 3,
                                                                       40877, 18705, 40937, 5172,
                                                                       5202, 19317, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 41997, 0, 3,
                                                                       40937, 18741, 40997, 5202,
                                                                       5232, 19377, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42097, 0, 3,
                                                                       40997, 18777, 41057, 5232,
                                                                       5262, 19437, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42197, 0, 3,
                                                                       41057, 18813, 41117, 5262,
                                                                       5292, 19497, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42297, 0, 3,
                                                                       41117, 18849, 41177, 5292,
                                                                       5322, 19557, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42397, 0, 3,
                                                                       41177, 18885, 41237, 5322,
                                                                       5352, 19617, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42497, 0, 3,
                                                                       41237, 18921, 41297, 5352,
                                                                       5382, 19677, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42597, 0, 3,
                                                                       41297, 18957, 41357, 5382,
                                                                       5412, 19737, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42697, 0, 3,
                                                                       41357, 18993, 41417, 5412,
                                                                       5442, 19797, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42797, 0, 3,
                                                                       41417, 19029, 41477, 5442,
                                                                       5472, 19857, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 42897, 0, 3,
                                                                       41477, 19065, 41537, 5472,
                                                                       5502, 19917, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 42997, 0, 3,
                                                                       41597, 19137, 41697, 5562,
                                                                       5607, 19977, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43147, 0, 3,
                                                                       41697, 19197, 41797, 5607,
                                                                       5652, 20067, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43297, 0, 3,
                                                                       41797, 19257, 41897, 5652,
                                                                       5697, 20157, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43447, 0, 3,
                                                                       41897, 19317, 41997, 5697,
                                                                       5742, 20247, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43597, 0, 3,
                                                                       41997, 19377, 42097, 5742,
                                                                       5787, 20337, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43747, 0, 3,
                                                                       42097, 19437, 42197, 5787,
                                                                       5832, 20427, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 43897, 0, 3,
                                                                       42197, 19497, 42297, 5832,
                                                                       5877, 20517, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44047, 0, 3,
                                                                       42297, 19557, 42397, 5877,
                                                                       5922, 20607, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44197, 0, 3,
                                                                       42397, 19617, 42497, 5922,
                                                                       5967, 20697, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44347, 0, 3,
                                                                       42497, 19677, 42597, 5967,
                                                                       6012, 20787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44497, 0, 3,
                                                                       42597, 19737, 42697, 6012,
                                                                       6057, 20877, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44647, 0, 3,
                                                                       42697, 19797, 42797, 6057,
                                                                       6102, 20967, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 44797, 0, 3,
                                                                       42797, 19857, 42897, 6102,
                                                                       6147, 21057, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 44947, 0, 3,
                                                                       42997, 19977, 43147, 6237,
                                                                       6300, 21147, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45157, 0, 3,
                                                                       43147, 20067, 43297, 6300,
                                                                       6363, 21273, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45367, 0, 3,
                                                                       43297, 20157, 43447, 6363,
                                                                       6426, 21399, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45577, 0, 3,
                                                                       43447, 20247, 43597, 6426,
                                                                       6489, 21525, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45787, 0, 3,
                                                                       43597, 20337, 43747, 6489,
                                                                       6552, 21651, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 45997, 0, 3,
                                                                       43747, 20427, 43897, 6552,
                                                                       6615, 21777, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46207, 0, 3,
                                                                       43897, 20517, 44047, 6615,
                                                                       6678, 21903, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46417, 0, 3,
                                                                       44047, 20607, 44197, 6678,
                                                                       6741, 22029, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46627, 0, 3,
                                                                       44197, 20697, 44347, 6741,
                                                                       6804, 22155, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 46837, 0, 3,
                                                                       44347, 20787, 44497, 6804,
                                                                       6867, 22281, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 47047, 0, 3,
                                                                       44497, 20877, 44647, 6867,
                                                                       6930, 22407, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 47257, 0, 3,
                                                                       44647, 20967, 44797, 6930,
                                                                       6993, 22533, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47467, 0, 3,
                                                                       44947, 21147, 45157, 7119,
                                                                       7203, 22659, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 47747, 0, 3,
                                                                       45157, 21273, 45367, 7203,
                                                                       7287, 22827, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48027, 0, 3,
                                                                       45367, 21399, 45577, 7287,
                                                                       7371, 22995, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48307, 0, 3,
                                                                       45577, 21525, 45787, 7371,
                                                                       7455, 23163, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48587, 0, 3,
                                                                       45787, 21651, 45997, 7455,
                                                                       7539, 23331, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 48867, 0, 3,
                                                                       45997, 21777, 46207, 7539,
                                                                       7623, 23499, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49147, 0, 3,
                                                                       46207, 21903, 46417, 7623,
                                                                       7707, 23667, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49427, 0, 3,
                                                                       46417, 22029, 46627, 7707,
                                                                       7791, 23835, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49707, 0, 3,
                                                                       46627, 22155, 46837, 7791,
                                                                       7875, 24003, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 49987, 0, 3,
                                                                       46837, 22281, 47047, 7875,
                                                                       7959, 24171, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 50267, 0, 3,
                                                                       47047, 22407, 47257, 7959,
                                                                       8043, 24339, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 50547, 0, 3,
                                                                       47467, 22659, 47747, 8211,
                                                                       8319, 24507, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 50907, 0, 3,
                                                                       47747, 22827, 48027, 8319,
                                                                       8427, 24723, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51267, 0, 3,
                                                                       48027, 22995, 48307, 8427,
                                                                       8535, 24939, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51627, 0, 3,
                                                                       48307, 23163, 48587, 8535,
                                                                       8643, 25155, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 51987, 0, 3,
                                                                       48587, 23331, 48867, 8643,
                                                                       8751, 25371, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52347, 0, 3,
                                                                       48867, 23499, 49147, 8751,
                                                                       8859, 25587, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 52707, 0, 3,
                                                                       49147, 23667, 49427, 8859,
                                                                       8967, 25803, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53067, 0, 3,
                                                                       49427, 23835, 49707, 8967,
                                                                       9075, 26019, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53427, 0, 3,
                                                                       49707, 24003, 49987, 9075,
                                                                       9183, 26235, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 53787, 0, 3,
                                                                       49987, 24171, 50267, 9183,
                                                                       9291, 26451, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 54147, 0, 3,
                                                                       50547, 24507, 50907, 9507,
                                                                       9642, 26667, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 54597, 0, 3,
                                                                       50907, 24723, 51267, 9642,
                                                                       9777, 26937, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55047, 0, 3,
                                                                       51267, 24939, 51627, 9777,
                                                                       9912, 27207, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55497, 0, 3,
                                                                       51627, 25155, 51987, 9912,
                                                                       10047, 27477, ncols,
                                                                       gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 55947, 0, 3,
                                                                       51987, 25371, 52347,
                                                                       10047, 10182, 27747,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 56397, 0, 3,
                                                                       52347, 25587, 52707,
                                                                       10182, 10317, 28017,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 56847, 0, 3,
                                                                       52707, 25803, 53067,
                                                                       10317, 10452, 28287,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 57297, 0, 3,
                                                                       53067, 26019, 53427,
                                                                       10452, 10587, 28557,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 57747, 0, 3,
                                                                       53427, 26235, 53787,
                                                                       10587, 10722, 28827,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 58197, 0, 3,
                                                                       54147, 26667, 54597,
                                                                       10992, 11157, 29097,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 58747, 0, 3,
                                                                       54597, 26937, 55047,
                                                                       11157, 11322, 29427,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 59297, 0, 3,
                                                                       55047, 27207, 55497,
                                                                       11322, 11487, 29757,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 59847, 0, 3,
                                                                       55497, 27477, 55947,
                                                                       11487, 11652, 30087,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 60397, 0, 3,
                                                                       55947, 27747, 56397,
                                                                       11652, 11817, 30417,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 60947, 0, 3,
                                                                       56397, 28017, 56847,
                                                                       11817, 11982, 30747,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 61497, 0, 3,
                                                                       56847, 28287, 57297,
                                                                       11982, 12147, 31077,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 62047, 0, 3,
                                                                       57297, 28557, 57747,
                                                                       12147, 12312, 31407,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 62597, 0, 3,
                                                                       58197, 29097, 58747,
                                                                       12642, 12840, 31737,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 63257, 0, 3,
                                                                       58747, 29427, 59297,
                                                                       12840, 13038, 32133,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 63917, 0, 3,
                                                                       59297, 29757, 59847,
                                                                       13038, 13236, 32529,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 64577, 0, 3,
                                                                       59847, 30087, 60397,
                                                                       13236, 13434, 32925,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 65237, 0, 3,
                                                                       60397, 30417, 60947,
                                                                       13434, 13632, 33321,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 65897, 0, 3,
                                                                       60947, 30747, 61497,
                                                                       13632, 13830, 33717,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 66557, 0, 3,
                                                                       61497, 31077, 62047,
                                                                       13830, 14028, 34113,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 67217, 0, 3,
                                                                       62597, 31737, 63257,
                                                                       14424, 14658, 34509,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 67997, 0, 3,
                                                                       63257, 32133, 63917,
                                                                       14658, 14892, 34977,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 68777, 0, 3,
                                                                       63917, 32529, 64577,
                                                                       14892, 15126, 35445,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 69557, 0, 3,
                                                                       64577, 32925, 65237,
                                                                       15126, 15360, 35913,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 70337, 0, 3,
                                                                       65237, 33321, 65897,
                                                                       15360, 15594, 36381,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 71117, 0, 3,
                                                                       65897, 33717, 66557,
                                                                       15594, 15828, 36849,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 71897, 0, 3,
                                                                       67217, 34509, 67997,
                                                                       16296, 16569, 37317,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 72807, 0, 3,
                                                                       67997, 34977, 68777,
                                                                       16569, 16842, 37863,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 73717, 0, 3,
                                                                       68777, 35445, 69557,
                                                                       16842, 17115, 38409,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 74627, 0, 3,
                                                                       69557, 35913, 70337,
                                                                       17115, 17388, 38955,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 75537, 0, 3,
                                                                       70337, 36381, 71117,
                                                                       17388, 17661, 39501,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76447, 3, 18207,
                                                                       18213, 40067, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76462, 3, 18213,
                                                                       18219, 40077, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76477, 3, 18219,
                                                                       18225, 40087, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76492, 3, 18225,
                                                                       18231, 40097, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76507, 3, 18231,
                                                                       18237, 40107, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76522, 3, 18237,
                                                                       18243, 40117, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76537, 3, 18243,
                                                                       18249, 40127, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76552, 3, 18249,
                                                                       18255, 40137, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76567, 3, 18255,
                                                                       18261, 40147, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76582, 3, 18261,
                                                                       18267, 40157, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76597, 3, 18267,
                                                                       18273, 40167, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76612, 3, 18273,
                                                                       18279, 40177, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76627, 3, 18279,
                                                                       18285, 40187, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76642, 3, 18285,
                                                                       18291, 40197, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 76657, 3, 18291,
                                                                       18297, 40207, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76672, 0, 3,
                                                                       76447, 40067, 76462,
                                                                       18309, 18327, 40277,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76717, 0, 3,
                                                                       76462, 40077, 76477,
                                                                       18327, 18345, 40307,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76762, 0, 3,
                                                                       76477, 40087, 76492,
                                                                       18345, 18363, 40337,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76807, 0, 3,
                                                                       76492, 40097, 76507,
                                                                       18363, 18381, 40367,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76852, 0, 3,
                                                                       76507, 40107, 76522,
                                                                       18381, 18399, 40397,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76897, 0, 3,
                                                                       76522, 40117, 76537,
                                                                       18399, 18417, 40427,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76942, 0, 3,
                                                                       76537, 40127, 76552,
                                                                       18417, 18435, 40457,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 76987, 0, 3,
                                                                       76552, 40137, 76567,
                                                                       18435, 18453, 40487,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77032, 0, 3,
                                                                       76567, 40147, 76582,
                                                                       18453, 18471, 40517,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77077, 0, 3,
                                                                       76582, 40157, 76597,
                                                                       18471, 18489, 40547,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77122, 0, 3,
                                                                       76597, 40167, 76612,
                                                                       18489, 18507, 40577,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77167, 0, 3,
                                                                       76612, 40177, 76627,
                                                                       18507, 18525, 40607,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77212, 0, 3,
                                                                       76627, 40187, 76642,
                                                                       18525, 18543, 40637,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 77257, 0, 3,
                                                                       76642, 40197, 76657,
                                                                       18543, 18561, 40667,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77302, 0, 3,
                                                                       76672, 40277, 76717,
                                                                       18597, 18633, 40817,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77392, 0, 3,
                                                                       76717, 40307, 76762,
                                                                       18633, 18669, 40877,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77482, 0, 3,
                                                                       76762, 40337, 76807,
                                                                       18669, 18705, 40937,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77572, 0, 3,
                                                                       76807, 40367, 76852,
                                                                       18705, 18741, 40997,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77662, 0, 3,
                                                                       76852, 40397, 76897,
                                                                       18741, 18777, 41057,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77752, 0, 3,
                                                                       76897, 40427, 76942,
                                                                       18777, 18813, 41117,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77842, 0, 3,
                                                                       76942, 40457, 76987,
                                                                       18813, 18849, 41177,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 77932, 0, 3,
                                                                       76987, 40487, 77032,
                                                                       18849, 18885, 41237,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78022, 0, 3,
                                                                       77032, 40517, 77077,
                                                                       18885, 18921, 41297,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78112, 0, 3,
                                                                       77077, 40547, 77122,
                                                                       18921, 18957, 41357,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78202, 0, 3,
                                                                       77122, 40577, 77167,
                                                                       18957, 18993, 41417,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78292, 0, 3,
                                                                       77167, 40607, 77212,
                                                                       18993, 19029, 41477,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 78382, 0, 3,
                                                                       77212, 40637, 77257,
                                                                       19029, 19065, 41537,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78472, 0, 3,
                                                                       77302, 40817, 77392,
                                                                       19137, 19197, 41797,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78622, 0, 3,
                                                                       77392, 40877, 77482,
                                                                       19197, 19257, 41897,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78772, 0, 3,
                                                                       77482, 40937, 77572,
                                                                       19257, 19317, 41997,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 78922, 0, 3,
                                                                       77572, 40997, 77662,
                                                                       19317, 19377, 42097,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79072, 0, 3,
                                                                       77662, 41057, 77752,
                                                                       19377, 19437, 42197,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79222, 0, 3,
                                                                       77752, 41117, 77842,
                                                                       19437, 19497, 42297,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79372, 0, 3,
                                                                       77842, 41177, 77932,
                                                                       19497, 19557, 42397,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79522, 0, 3,
                                                                       77932, 41237, 78022,
                                                                       19557, 19617, 42497,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79672, 0, 3,
                                                                       78022, 41297, 78112,
                                                                       19617, 19677, 42597,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79822, 0, 3,
                                                                       78112, 41357, 78202,
                                                                       19677, 19737, 42697,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 79972, 0, 3,
                                                                       78202, 41417, 78292,
                                                                       19737, 19797, 42797,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 80122, 0, 3,
                                                                       78292, 41477, 78382,
                                                                       19797, 19857, 42897,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80272, 0, 3,
                                                                       78472, 41797, 78622,
                                                                       19977, 20067, 43297,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80497, 0, 3,
                                                                       78622, 41897, 78772,
                                                                       20067, 20157, 43447,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80722, 0, 3,
                                                                       78772, 41997, 78922,
                                                                       20157, 20247, 43597,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 80947, 0, 3,
                                                                       78922, 42097, 79072,
                                                                       20247, 20337, 43747,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81172, 0, 3,
                                                                       79072, 42197, 79222,
                                                                       20337, 20427, 43897,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81397, 0, 3,
                                                                       79222, 42297, 79372,
                                                                       20427, 20517, 44047,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81622, 0, 3,
                                                                       79372, 42397, 79522,
                                                                       20517, 20607, 44197,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 81847, 0, 3,
                                                                       79522, 42497, 79672,
                                                                       20607, 20697, 44347,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82072, 0, 3,
                                                                       79672, 42597, 79822,
                                                                       20697, 20787, 44497,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82297, 0, 3,
                                                                       79822, 42697, 79972,
                                                                       20787, 20877, 44647,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 82522, 0, 3,
                                                                       79972, 42797, 80122,
                                                                       20877, 20967, 44797,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 82747, 0, 3,
                                                                       80272, 43297, 80497,
                                                                       21147, 21273, 45367,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83062, 0, 3,
                                                                       80497, 43447, 80722,
                                                                       21273, 21399, 45577,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83377, 0, 3,
                                                                       80722, 43597, 80947,
                                                                       21399, 21525, 45787,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 83692, 0, 3,
                                                                       80947, 43747, 81172,
                                                                       21525, 21651, 45997,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84007, 0, 3,
                                                                       81172, 43897, 81397,
                                                                       21651, 21777, 46207,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84322, 0, 3,
                                                                       81397, 44047, 81622,
                                                                       21777, 21903, 46417,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84637, 0, 3,
                                                                       81622, 44197, 81847,
                                                                       21903, 22029, 46627,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 84952, 0, 3,
                                                                       81847, 44347, 82072,
                                                                       22029, 22155, 46837,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 85267, 0, 3,
                                                                       82072, 44497, 82297,
                                                                       22155, 22281, 47047,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 85582, 0, 3,
                                                                       82297, 44647, 82522,
                                                                       22281, 22407, 47257,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 85897, 0, 3,
                                                                       82747, 45367, 83062,
                                                                       22659, 22827, 48027,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 86317, 0, 3,
                                                                       83062, 45577, 83377,
                                                                       22827, 22995, 48307,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 86737, 0, 3,
                                                                       83377, 45787, 83692,
                                                                       22995, 23163, 48587,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87157, 0, 3,
                                                                       83692, 45997, 84007,
                                                                       23163, 23331, 48867,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87577, 0, 3,
                                                                       84007, 46207, 84322,
                                                                       23331, 23499, 49147,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 87997, 0, 3,
                                                                       84322, 46417, 84637,
                                                                       23499, 23667, 49427,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 88417, 0, 3,
                                                                       84637, 46627, 84952,
                                                                       23667, 23835, 49707,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 88837, 0, 3,
                                                                       84952, 46837, 85267,
                                                                       23835, 24003, 49987,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 89257, 0, 3,
                                                                       85267, 47047, 85582,
                                                                       24003, 24171, 50267,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 89677, 0, 3,
                                                                       85897, 48027, 86317,
                                                                       24507, 24723, 51267,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 90217, 0, 3,
                                                                       86317, 48307, 86737,
                                                                       24723, 24939, 51627,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 90757, 0, 3,
                                                                       86737, 48587, 87157,
                                                                       24939, 25155, 51987,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 91297, 0, 3,
                                                                       87157, 48867, 87577,
                                                                       25155, 25371, 52347,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 91837, 0, 3,
                                                                       87577, 49147, 87997,
                                                                       25371, 25587, 52707,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 92377, 0, 3,
                                                                       87997, 49427, 88417,
                                                                       25587, 25803, 53067,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 92917, 0, 3,
                                                                       88417, 49707, 88837,
                                                                       25803, 26019, 53427,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 93457, 0, 3,
                                                                       88837, 49987, 89257,
                                                                       26019, 26235, 53787,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 93997, 0, 3,
                                                                       89677, 51267, 90217,
                                                                       26667, 26937, 55047,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 94672, 0, 3,
                                                                       90217, 51627, 90757,
                                                                       26937, 27207, 55497,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 95347, 0, 3,
                                                                       90757, 51987, 91297,
                                                                       27207, 27477, 55947,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 96022, 0, 3,
                                                                       91297, 52347, 91837,
                                                                       27477, 27747, 56397,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 96697, 0, 3,
                                                                       91837, 52707, 92377,
                                                                       27747, 28017, 56847,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 97372, 0, 3,
                                                                       92377, 53067, 92917,
                                                                       28017, 28287, 57297,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 98047, 0, 3,
                                                                       92917, 53427, 93457,
                                                                       28287, 28557, 57747,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 98722, 0, 3,
                                                                       93997, 55047, 94672,
                                                                       29097, 29427, 59297,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 99547, 0, 3,
                                                                       94672, 55497, 95347,
                                                                       29427, 29757, 59847,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 100372, 0, 3,
                                                                       95347, 55947, 96022,
                                                                       29757, 30087, 60397,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 101197, 0, 3,
                                                                       96022, 56397, 96697,
                                                                       30087, 30417, 60947,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 102022, 0, 3,
                                                                       96697, 56847, 97372,
                                                                       30417, 30747, 61497,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 102847, 0, 3,
                                                                       97372, 57297, 98047,
                                                                       30747, 31077, 62047,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 103672, 0, 3,
                                                                       98722, 59297, 99547,
                                                                       31737, 32133, 63917,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 104662, 0, 3,
                                                                       99547, 59847, 100372,
                                                                       32133, 32529, 64577,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 105652, 0, 3,
                                                                       100372, 60397, 101197,
                                                                       32529, 32925, 65237,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 106642, 0, 3,
                                                                       101197, 60947, 102022,
                                                                       32925, 33321, 65897,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 107632, 0, 3,
                                                                       102022, 61497, 102847,
                                                                       33321, 33717, 66557,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 108622, 0, 3,
                                                                       103672, 63917, 104662,
                                                                       34509, 34977, 68777,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 109792, 0, 3,
                                                                       104662, 64577, 105652,
                                                                       34977, 35445, 69557,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 110962, 0, 3,
                                                                       105652, 65237, 106642,
                                                                       35445, 35913, 70337,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 112132, 0, 3,
                                                                       106642, 65897, 107632,
                                                                       35913, 36381, 71117,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 113302, 0, 3,
                                                                       108622, 68777, 109792,
                                                                       37317, 37863, 73717,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 114667, 0, 3,
                                                                       109792, 69557, 110962,
                                                                       37863, 38409, 74627,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 116032, 0, 3,
                                                                       110962, 70337, 112132,
                                                                       38409, 38955, 75537,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117397, 3, 40047,
                                                                       40057, 76447, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117418, 3, 40057,
                                                                       40067, 76462, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117439, 3, 40067,
                                                                       40077, 76477, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117460, 3, 40077,
                                                                       40087, 76492, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117481, 3, 40087,
                                                                       40097, 76507, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117502, 3, 40097,
                                                                       40107, 76522, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117523, 3, 40107,
                                                                       40117, 76537, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117544, 3, 40117,
                                                                       40127, 76552, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117565, 3, 40127,
                                                                       40137, 76567, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117586, 3, 40137,
                                                                       40147, 76582, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117607, 3, 40147,
                                                                       40157, 76597, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117628, 3, 40157,
                                                                       40167, 76612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117649, 3, 40167,
                                                                       40177, 76627, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117670, 3, 40177,
                                                                       40187, 76642, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 117691, 3, 40187,
                                                                       40197, 76657, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117712, 0, 3,
                                                                       117397, 76447, 117418,
                                                                       40217, 40247, 76672,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117775, 0, 3,
                                                                       117418, 76462, 117439,
                                                                       40247, 40277, 76717,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117838, 0, 3,
                                                                       117439, 76477, 117460,
                                                                       40277, 40307, 76762,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117901, 0, 3,
                                                                       117460, 76492, 117481,
                                                                       40307, 40337, 76807,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 117964, 0, 3,
                                                                       117481, 76507, 117502,
                                                                       40337, 40367, 76852,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118027, 0, 3,
                                                                       117502, 76522, 117523,
                                                                       40367, 40397, 76897,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118090, 0, 3,
                                                                       117523, 76537, 117544,
                                                                       40397, 40427, 76942,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118153, 0, 3,
                                                                       117544, 76552, 117565,
                                                                       40427, 40457, 76987,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118216, 0, 3,
                                                                       117565, 76567, 117586,
                                                                       40457, 40487, 77032,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118279, 0, 3,
                                                                       117586, 76582, 117607,
                                                                       40487, 40517, 77077,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118342, 0, 3,
                                                                       117607, 76597, 117628,
                                                                       40517, 40547, 77122,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118405, 0, 3,
                                                                       117628, 76612, 117649,
                                                                       40547, 40577, 77167,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118468, 0, 3,
                                                                       117649, 76627, 117670,
                                                                       40577, 40607, 77212,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 118531, 0, 3,
                                                                       117670, 76642, 117691,
                                                                       40607, 40637, 77257,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 118594, 0, 3,
                                                                       117712, 76672, 117775,
                                                                       40697, 40757, 77302,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 118720, 0, 3,
                                                                       117775, 76717, 117838,
                                                                       40757, 40817, 77392,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 118846, 0, 3,
                                                                       117838, 76762, 117901,
                                                                       40817, 40877, 77482,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 118972, 0, 3,
                                                                       117901, 76807, 117964,
                                                                       40877, 40937, 77572,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119098, 0, 3,
                                                                       117964, 76852, 118027,
                                                                       40937, 40997, 77662,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119224, 0, 3,
                                                                       118027, 76897, 118090,
                                                                       40997, 41057, 77752,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119350, 0, 3,
                                                                       118090, 76942, 118153,
                                                                       41057, 41117, 77842,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119476, 0, 3,
                                                                       118153, 76987, 118216,
                                                                       41117, 41177, 77932,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119602, 0, 3,
                                                                       118216, 77032, 118279,
                                                                       41177, 41237, 78022,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119728, 0, 3,
                                                                       118279, 77077, 118342,
                                                                       41237, 41297, 78112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119854, 0, 3,
                                                                       118342, 77122, 118405,
                                                                       41297, 41357, 78202,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 119980, 0, 3,
                                                                       118405, 77167, 118468,
                                                                       41357, 41417, 78292,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 120106, 0, 3,
                                                                       118468, 77212, 118531,
                                                                       41417, 41477, 78382,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 120232, 0, 3,
                                                                       118594, 77302, 118720,
                                                                       41597, 41697, 78472,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 120442, 0, 3,
                                                                       118720, 77392, 118846,
                                                                       41697, 41797, 78622,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 120652, 0, 3,
                                                                       118846, 77482, 118972,
                                                                       41797, 41897, 78772,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 120862, 0, 3,
                                                                       118972, 77572, 119098,
                                                                       41897, 41997, 78922,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121072, 0, 3,
                                                                       119098, 77662, 119224,
                                                                       41997, 42097, 79072,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121282, 0, 3,
                                                                       119224, 77752, 119350,
                                                                       42097, 42197, 79222,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121492, 0, 3,
                                                                       119350, 77842, 119476,
                                                                       42197, 42297, 79372,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121702, 0, 3,
                                                                       119476, 77932, 119602,
                                                                       42297, 42397, 79522,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 121912, 0, 3,
                                                                       119602, 78022, 119728,
                                                                       42397, 42497, 79672,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 122122, 0, 3,
                                                                       119728, 78112, 119854,
                                                                       42497, 42597, 79822,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 122332, 0, 3,
                                                                       119854, 78202, 119980,
                                                                       42597, 42697, 79972,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 122542, 0, 3,
                                                                       119980, 78292, 120106,
                                                                       42697, 42797, 80122,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 122752, 0, 3,
                                                                       120232, 78472, 120442,
                                                                       42997, 43147, 80272,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 123067, 0, 3,
                                                                       120442, 78622, 120652,
                                                                       43147, 43297, 80497,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 123382, 0, 3,
                                                                       120652, 78772, 120862,
                                                                       43297, 43447, 80722,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 123697, 0, 3,
                                                                       120862, 78922, 121072,
                                                                       43447, 43597, 80947,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 124012, 0, 3,
                                                                       121072, 79072, 121282,
                                                                       43597, 43747, 81172,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 124327, 0, 3,
                                                                       121282, 79222, 121492,
                                                                       43747, 43897, 81397,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 124642, 0, 3,
                                                                       121492, 79372, 121702,
                                                                       43897, 44047, 81622,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 124957, 0, 3,
                                                                       121702, 79522, 121912,
                                                                       44047, 44197, 81847,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 125272, 0, 3,
                                                                       121912, 79672, 122122,
                                                                       44197, 44347, 82072,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 125587, 0, 3,
                                                                       122122, 79822, 122332,
                                                                       44347, 44497, 82297,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 125902, 0, 3,
                                                                       122332, 79972, 122542,
                                                                       44497, 44647, 82522,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 126217, 0, 3,
                                                                       122752, 80272, 123067,
                                                                       44947, 45157, 82747,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 126658, 0, 3,
                                                                       123067, 80497, 123382,
                                                                       45157, 45367, 83062,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 127099, 0, 3,
                                                                       123382, 80722, 123697,
                                                                       45367, 45577, 83377,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 127540, 0, 3,
                                                                       123697, 80947, 124012,
                                                                       45577, 45787, 83692,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 127981, 0, 3,
                                                                       124012, 81172, 124327,
                                                                       45787, 45997, 84007,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 128422, 0, 3,
                                                                       124327, 81397, 124642,
                                                                       45997, 46207, 84322,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 128863, 0, 3,
                                                                       124642, 81622, 124957,
                                                                       46207, 46417, 84637,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 129304, 0, 3,
                                                                       124957, 81847, 125272,
                                                                       46417, 46627, 84952,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 129745, 0, 3,
                                                                       125272, 82072, 125587,
                                                                       46627, 46837, 85267,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 130186, 0, 3,
                                                                       125587, 82297, 125902,
                                                                       46837, 47047, 85582,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 130627, 0, 3,
                                                                       126217, 82747, 126658,
                                                                       47467, 47747, 85897,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 131215, 0, 3,
                                                                       126658, 83062, 127099,
                                                                       47747, 48027, 86317,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 131803, 0, 3,
                                                                       127099, 83377, 127540,
                                                                       48027, 48307, 86737,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 132391, 0, 3,
                                                                       127540, 83692, 127981,
                                                                       48307, 48587, 87157,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 132979, 0, 3,
                                                                       127981, 84007, 128422,
                                                                       48587, 48867, 87577,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 133567, 0, 3,
                                                                       128422, 84322, 128863,
                                                                       48867, 49147, 87997,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 134155, 0, 3,
                                                                       128863, 84637, 129304,
                                                                       49147, 49427, 88417,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 134743, 0, 3,
                                                                       129304, 84952, 129745,
                                                                       49427, 49707, 88837,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 135331, 0, 3,
                                                                       129745, 85267, 130186,
                                                                       49707, 49987, 89257,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 135919, 0, 3,
                                                                       130627, 85897, 131215,
                                                                       50547, 50907, 89677,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 136675, 0, 3,
                                                                       131215, 86317, 131803,
                                                                       50907, 51267, 90217,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 137431, 0, 3,
                                                                       131803, 86737, 132391,
                                                                       51267, 51627, 90757,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 138187, 0, 3,
                                                                       132391, 87157, 132979,
                                                                       51627, 51987, 91297,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 138943, 0, 3,
                                                                       132979, 87577, 133567,
                                                                       51987, 52347, 91837,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 139699, 0, 3,
                                                                       133567, 87997, 134155,
                                                                       52347, 52707, 92377,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 140455, 0, 3,
                                                                       134155, 88417, 134743,
                                                                       52707, 53067, 92917,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 141211, 0, 3,
                                                                       134743, 88837, 135331,
                                                                       53067, 53427, 93457,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 141967, 0, 3,
                                                                       135919, 89677, 136675,
                                                                       54147, 54597, 93997,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 142912, 0, 3,
                                                                       136675, 90217, 137431,
                                                                       54597, 55047, 94672,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 143857, 0, 3,
                                                                       137431, 90757, 138187,
                                                                       55047, 55497, 95347,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 144802, 0, 3,
                                                                       138187, 91297, 138943,
                                                                       55497, 55947, 96022,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 145747, 0, 3,
                                                                       138943, 91837, 139699,
                                                                       55947, 56397, 96697,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 146692, 0, 3,
                                                                       139699, 92377, 140455,
                                                                       56397, 56847, 97372,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 147637, 0, 3,
                                                                       140455, 92917, 141211,
                                                                       56847, 57297, 98047,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 148582, 0, 3,
                                                                       141967, 93997, 142912,
                                                                       58197, 58747, 98722,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 149737, 0, 3,
                                                                       142912, 94672, 143857,
                                                                       58747, 59297, 99547,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 150892, 0, 3,
                                                                       143857, 95347, 144802,
                                                                       59297, 59847, 100372,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 152047, 0, 3,
                                                                       144802, 96022, 145747,
                                                                       59847, 60397, 101197,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 153202, 0, 3,
                                                                       145747, 96697, 146692,
                                                                       60397, 60947, 102022,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 154357, 0, 3,
                                                                       146692, 97372, 147637,
                                                                       60947, 61497, 102847,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 155512, 0, 3,
                                                                       148582, 98722, 149737,
                                                                       62597, 63257, 103672,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 156898, 0, 3,
                                                                       149737, 99547, 150892,
                                                                       63257, 63917, 104662,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 158284, 0, 3,
                                                                       150892, 100372, 152047,
                                                                       63917, 64577, 105652,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 159670, 0, 3,
                                                                       152047, 101197, 153202,
                                                                       64577, 65237, 106642,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 161056, 0, 3,
                                                                       153202, 102022, 154357,
                                                                       65237, 65897, 107632,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 162442, 0, 3,
                                                                       155512, 103672, 156898,
                                                                       67217, 67997, 108622,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 164080, 0, 3,
                                                                       156898, 104662, 158284,
                                                                       67997, 68777, 109792,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 165718, 0, 3,
                                                                       158284, 105652, 159670,
                                                                       68777, 69557, 110962,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 167356, 0, 3,
                                                                       159670, 106642, 161056,
                                                                       69557, 70337, 112132,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 168994, 0, 3,
                                                                       162442, 108622, 164080,
                                                                       71897, 72807, 113302,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 170905, 0, 3,
                                                                       164080, 109792, 165718,
                                                                       72807, 73717, 114667,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 172816, 0, 3,
                                                                       165718, 110962, 167356,
                                                                       73717, 74627, 116032,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174727, 3, 76447,
                                                                       76462, 117439, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174755, 3, 76462,
                                                                       76477, 117460, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174783, 3, 76477,
                                                                       76492, 117481, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174811, 3, 76492,
                                                                       76507, 117502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174839, 3, 76507,
                                                                       76522, 117523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174867, 3, 76522,
                                                                       76537, 117544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174895, 3, 76537,
                                                                       76552, 117565, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174923, 3, 76552,
                                                                       76567, 117586, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174951, 3, 76567,
                                                                       76582, 117607, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 174979, 3, 76582,
                                                                       76597, 117628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 175007, 3, 76597,
                                                                       76612, 117649, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 175035, 3, 76612,
                                                                       76627, 117670, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 175063, 3, 76627,
                                                                       76642, 117691, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175091, 0, 3,
                                                                       174727, 117439, 174755,
                                                                       76672, 76717, 117838,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175175, 0, 3,
                                                                       174755, 117460, 174783,
                                                                       76717, 76762, 117901,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175259, 0, 3,
                                                                       174783, 117481, 174811,
                                                                       76762, 76807, 117964,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175343, 0, 3,
                                                                       174811, 117502, 174839,
                                                                       76807, 76852, 118027,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175427, 0, 3,
                                                                       174839, 117523, 174867,
                                                                       76852, 76897, 118090,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175511, 0, 3,
                                                                       174867, 117544, 174895,
                                                                       76897, 76942, 118153,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175595, 0, 3,
                                                                       174895, 117565, 174923,
                                                                       76942, 76987, 118216,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175679, 0, 3,
                                                                       174923, 117586, 174951,
                                                                       76987, 77032, 118279,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175763, 0, 3,
                                                                       174951, 117607, 174979,
                                                                       77032, 77077, 118342,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175847, 0, 3,
                                                                       174979, 117628, 175007,
                                                                       77077, 77122, 118405,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 175931, 0, 3,
                                                                       175007, 117649, 175035,
                                                                       77122, 77167, 118468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 176015, 0, 3,
                                                                       175035, 117670, 175063,
                                                                       77167, 77212, 118531,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176099, 0, 3,
                                                                       175091, 117838, 175175,
                                                                       77302, 77392, 118846,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176267, 0, 3,
                                                                       175175, 117901, 175259,
                                                                       77392, 77482, 118972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176435, 0, 3,
                                                                       175259, 117964, 175343,
                                                                       77482, 77572, 119098,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176603, 0, 3,
                                                                       175343, 118027, 175427,
                                                                       77572, 77662, 119224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176771, 0, 3,
                                                                       175427, 118090, 175511,
                                                                       77662, 77752, 119350,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 176939, 0, 3,
                                                                       175511, 118153, 175595,
                                                                       77752, 77842, 119476,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177107, 0, 3,
                                                                       175595, 118216, 175679,
                                                                       77842, 77932, 119602,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177275, 0, 3,
                                                                       175679, 118279, 175763,
                                                                       77932, 78022, 119728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177443, 0, 3,
                                                                       175763, 118342, 175847,
                                                                       78022, 78112, 119854,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177611, 0, 3,
                                                                       175847, 118405, 175931,
                                                                       78112, 78202, 119980,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 177779, 0, 3,
                                                                       175931, 118468, 176015,
                                                                       78202, 78292, 120106,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 177947, 0, 3,
                                                                       176099, 118846, 176267,
                                                                       78472, 78622, 120652,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178227, 0, 3,
                                                                       176267, 118972, 176435,
                                                                       78622, 78772, 120862,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178507, 0, 3,
                                                                       176435, 119098, 176603,
                                                                       78772, 78922, 121072,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 178787, 0, 3,
                                                                       176603, 119224, 176771,
                                                                       78922, 79072, 121282,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179067, 0, 3,
                                                                       176771, 119350, 176939,
                                                                       79072, 79222, 121492,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179347, 0, 3,
                                                                       176939, 119476, 177107,
                                                                       79222, 79372, 121702,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179627, 0, 3,
                                                                       177107, 119602, 177275,
                                                                       79372, 79522, 121912,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 179907, 0, 3,
                                                                       177275, 119728, 177443,
                                                                       79522, 79672, 122122,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 180187, 0, 3,
                                                                       177443, 119854, 177611,
                                                                       79672, 79822, 122332,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 180467, 0, 3,
                                                                       177611, 119980, 177779,
                                                                       79822, 79972, 122542,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 180747, 0, 3,
                                                                       177947, 120652, 178227,
                                                                       80272, 80497, 123382,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 181167, 0, 3,
                                                                       178227, 120862, 178507,
                                                                       80497, 80722, 123697,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 181587, 0, 3,
                                                                       178507, 121072, 178787,
                                                                       80722, 80947, 124012,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182007, 0, 3,
                                                                       178787, 121282, 179067,
                                                                       80947, 81172, 124327,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182427, 0, 3,
                                                                       179067, 121492, 179347,
                                                                       81172, 81397, 124642,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 182847, 0, 3,
                                                                       179347, 121702, 179627,
                                                                       81397, 81622, 124957,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 183267, 0, 3,
                                                                       179627, 121912, 179907,
                                                                       81622, 81847, 125272,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 183687, 0, 3,
                                                                       179907, 122122, 180187,
                                                                       81847, 82072, 125587,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 184107, 0, 3,
                                                                       180187, 122332, 180467,
                                                                       82072, 82297, 125902,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 184527, 0, 3,
                                                                       180747, 123382, 181167,
                                                                       82747, 83062, 127099,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 185115, 0, 3,
                                                                       181167, 123697, 181587,
                                                                       83062, 83377, 127540,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 185703, 0, 3,
                                                                       181587, 124012, 182007,
                                                                       83377, 83692, 127981,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 186291, 0, 3,
                                                                       182007, 124327, 182427,
                                                                       83692, 84007, 128422,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 186879, 0, 3,
                                                                       182427, 124642, 182847,
                                                                       84007, 84322, 128863,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 187467, 0, 3,
                                                                       182847, 124957, 183267,
                                                                       84322, 84637, 129304,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 188055, 0, 3,
                                                                       183267, 125272, 183687,
                                                                       84637, 84952, 129745,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 188643, 0, 3,
                                                                       183687, 125587, 184107,
                                                                       84952, 85267, 130186,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 189231, 0, 3,
                                                                       184527, 127099, 185115,
                                                                       85897, 86317, 131803,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 190015, 0, 3,
                                                                       185115, 127540, 185703,
                                                                       86317, 86737, 132391,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 190799, 0, 3,
                                                                       185703, 127981, 186291,
                                                                       86737, 87157, 132979,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 191583, 0, 3,
                                                                       186291, 128422, 186879,
                                                                       87157, 87577, 133567,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 192367, 0, 3,
                                                                       186879, 128863, 187467,
                                                                       87577, 87997, 134155,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 193151, 0, 3,
                                                                       187467, 129304, 188055,
                                                                       87997, 88417, 134743,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 193935, 0, 3,
                                                                       188055, 129745, 188643,
                                                                       88417, 88837, 135331,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 194719, 0, 3,
                                                                       189231, 131803, 190015,
                                                                       89677, 90217, 137431,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 195727, 0, 3,
                                                                       190015, 132391, 190799,
                                                                       90217, 90757, 138187,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 196735, 0, 3,
                                                                       190799, 132979, 191583,
                                                                       90757, 91297, 138943,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 197743, 0, 3,
                                                                       191583, 133567, 192367,
                                                                       91297, 91837, 139699,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 198751, 0, 3,
                                                                       192367, 134155, 193151,
                                                                       91837, 92377, 140455,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 199759, 0, 3,
                                                                       193151, 134743, 193935,
                                                                       92377, 92917, 141211,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 200767, 0, 3,
                                                                       194719, 137431, 195727,
                                                                       93997, 94672, 143857,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 202027, 0, 3,
                                                                       195727, 138187, 196735,
                                                                       94672, 95347, 144802,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 203287, 0, 3,
                                                                       196735, 138943, 197743,
                                                                       95347, 96022, 145747,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 204547, 0, 3,
                                                                       197743, 139699, 198751,
                                                                       96022, 96697, 146692,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 205807, 0, 3,
                                                                       198751, 140455, 199759,
                                                                       96697, 97372, 147637,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 207067, 0, 3,
                                                                       200767, 143857, 202027,
                                                                       98722, 99547, 150892,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 208607, 0, 3,
                                                                       202027, 144802, 203287,
                                                                       99547, 100372, 152047,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 210147, 0, 3,
                                                                       203287, 145747, 204547,
                                                                       100372, 101197, 153202,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 211687, 0, 3,
                                                                       204547, 146692, 205807,
                                                                       101197, 102022, 154357,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 213227, 0, 3,
                                                                       207067, 150892, 208607,
                                                                       103672, 104662, 158284,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 215075, 0, 3,
                                                                       208607, 152047, 210147,
                                                                       104662, 105652, 159670,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 216923, 0, 3,
                                                                       210147, 153202, 211687,
                                                                       105652, 106642, 161056,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 218771, 0, 3,
                                                                       213227, 158284, 215075,
                                                                       108622, 109792, 165718,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 220955, 0, 3,
                                                                       215075, 159670, 216923,
                                                                       109792, 110962, 167356,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 223139, 0, 3,
                                                                       218771, 165718, 220955,
                                                                       113302, 114667, 172816,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225687, 3, 117397,
                                                                       117418, 174727, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225723, 3, 117418,
                                                                       117439, 174755, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225759, 3, 117439,
                                                                       117460, 174783, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225795, 3, 117460,
                                                                       117481, 174811, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225831, 3, 117481,
                                                                       117502, 174839, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225867, 3, 117502,
                                                                       117523, 174867, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225903, 3, 117523,
                                                                       117544, 174895, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225939, 3, 117544,
                                                                       117565, 174923, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 225975, 3, 117565,
                                                                       117586, 174951, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 226011, 3, 117586,
                                                                       117607, 174979, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 226047, 3, 117607,
                                                                       117628, 175007, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 226083, 3, 117628,
                                                                       117649, 175035, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 226119, 3, 117649,
                                                                       117670, 175063, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226155, 0, 3,
                                                                       225687, 174727, 225723,
                                                                       117712, 117775, 175091,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226263, 0, 3,
                                                                       225723, 174755, 225759,
                                                                       117775, 117838, 175175,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226371, 0, 3,
                                                                       225759, 174783, 225795,
                                                                       117838, 117901, 175259,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226479, 0, 3,
                                                                       225795, 174811, 225831,
                                                                       117901, 117964, 175343,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226587, 0, 3,
                                                                       225831, 174839, 225867,
                                                                       117964, 118027, 175427,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226695, 0, 3,
                                                                       225867, 174867, 225903,
                                                                       118027, 118090, 175511,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226803, 0, 3,
                                                                       225903, 174895, 225939,
                                                                       118090, 118153, 175595,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 226911, 0, 3,
                                                                       225939, 174923, 225975,
                                                                       118153, 118216, 175679,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 227019, 0, 3,
                                                                       225975, 174951, 226011,
                                                                       118216, 118279, 175763,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 227127, 0, 3,
                                                                       226011, 174979, 226047,
                                                                       118279, 118342, 175847,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 227235, 0, 3,
                                                                       226047, 175007, 226083,
                                                                       118342, 118405, 175931,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 227343, 0, 3,
                                                                       226083, 175035, 226119,
                                                                       118405, 118468, 176015,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 227451, 0, 3,
                                                                       226155, 175091, 226263,
                                                                       118594, 118720, 176099,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 227667, 0, 3,
                                                                       226263, 175175, 226371,
                                                                       118720, 118846, 176267,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 227883, 0, 3,
                                                                       226371, 175259, 226479,
                                                                       118846, 118972, 176435,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228099, 0, 3,
                                                                       226479, 175343, 226587,
                                                                       118972, 119098, 176603,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228315, 0, 3,
                                                                       226587, 175427, 226695,
                                                                       119098, 119224, 176771,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228531, 0, 3,
                                                                       226695, 175511, 226803,
                                                                       119224, 119350, 176939,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228747, 0, 3,
                                                                       226803, 175595, 226911,
                                                                       119350, 119476, 177107,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 228963, 0, 3,
                                                                       226911, 175679, 227019,
                                                                       119476, 119602, 177275,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 229179, 0, 3,
                                                                       227019, 175763, 227127,
                                                                       119602, 119728, 177443,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 229395, 0, 3,
                                                                       227127, 175847, 227235,
                                                                       119728, 119854, 177611,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 229611, 0, 3,
                                                                       227235, 175931, 227343,
                                                                       119854, 119980, 177779,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 229827, 0, 3,
                                                                       227451, 176099, 227667,
                                                                       120232, 120442, 177947,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 230187, 0, 3,
                                                                       227667, 176267, 227883,
                                                                       120442, 120652, 178227,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 230547, 0, 3,
                                                                       227883, 176435, 228099,
                                                                       120652, 120862, 178507,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 230907, 0, 3,
                                                                       228099, 176603, 228315,
                                                                       120862, 121072, 178787,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 231267, 0, 3,
                                                                       228315, 176771, 228531,
                                                                       121072, 121282, 179067,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 231627, 0, 3,
                                                                       228531, 176939, 228747,
                                                                       121282, 121492, 179347,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 231987, 0, 3,
                                                                       228747, 177107, 228963,
                                                                       121492, 121702, 179627,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 232347, 0, 3,
                                                                       228963, 177275, 229179,
                                                                       121702, 121912, 179907,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 232707, 0, 3,
                                                                       229179, 177443, 229395,
                                                                       121912, 122122, 180187,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 233067, 0, 3,
                                                                       229395, 177611, 229611,
                                                                       122122, 122332, 180467,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 233427, 0, 3,
                                                                       229827, 177947, 230187,
                                                                       122752, 123067, 180747,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 233967, 0, 3,
                                                                       230187, 178227, 230547,
                                                                       123067, 123382, 181167,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 234507, 0, 3,
                                                                       230547, 178507, 230907,
                                                                       123382, 123697, 181587,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 235047, 0, 3,
                                                                       230907, 178787, 231267,
                                                                       123697, 124012, 182007,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 235587, 0, 3,
                                                                       231267, 179067, 231627,
                                                                       124012, 124327, 182427,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 236127, 0, 3,
                                                                       231627, 179347, 231987,
                                                                       124327, 124642, 182847,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 236667, 0, 3,
                                                                       231987, 179627, 232347,
                                                                       124642, 124957, 183267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 237207, 0, 3,
                                                                       232347, 179907, 232707,
                                                                       124957, 125272, 183687,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 237747, 0, 3,
                                                                       232707, 180187, 233067,
                                                                       125272, 125587, 184107,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 238287, 0, 3,
                                                                       233427, 180747, 233967,
                                                                       126217, 126658, 184527,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 239043, 0, 3,
                                                                       233967, 181167, 234507,
                                                                       126658, 127099, 185115,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 239799, 0, 3,
                                                                       234507, 181587, 235047,
                                                                       127099, 127540, 185703,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 240555, 0, 3,
                                                                       235047, 182007, 235587,
                                                                       127540, 127981, 186291,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 241311, 0, 3,
                                                                       235587, 182427, 236127,
                                                                       127981, 128422, 186879,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 242067, 0, 3,
                                                                       236127, 182847, 236667,
                                                                       128422, 128863, 187467,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 242823, 0, 3,
                                                                       236667, 183267, 237207,
                                                                       128863, 129304, 188055,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 243579, 0, 3,
                                                                       237207, 183687, 237747,
                                                                       129304, 129745, 188643,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 244335, 0, 3,
                                                                       238287, 184527, 239043,
                                                                       130627, 131215, 189231,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 245343, 0, 3,
                                                                       239043, 185115, 239799,
                                                                       131215, 131803, 190015,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 246351, 0, 3,
                                                                       239799, 185703, 240555,
                                                                       131803, 132391, 190799,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 247359, 0, 3,
                                                                       240555, 186291, 241311,
                                                                       132391, 132979, 191583,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 248367, 0, 3,
                                                                       241311, 186879, 242067,
                                                                       132979, 133567, 192367,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 249375, 0, 3,
                                                                       242067, 187467, 242823,
                                                                       133567, 134155, 193151,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 250383, 0, 3,
                                                                       242823, 188055, 243579,
                                                                       134155, 134743, 193935,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 251391, 0, 3,
                                                                       244335, 189231, 245343,
                                                                       135919, 136675, 194719,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 252687, 0, 3,
                                                                       245343, 190015, 246351,
                                                                       136675, 137431, 195727,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 253983, 0, 3,
                                                                       246351, 190799, 247359,
                                                                       137431, 138187, 196735,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 255279, 0, 3,
                                                                       247359, 191583, 248367,
                                                                       138187, 138943, 197743,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 256575, 0, 3,
                                                                       248367, 192367, 249375,
                                                                       138943, 139699, 198751,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 257871, 0, 3,
                                                                       249375, 193151, 250383,
                                                                       139699, 140455, 199759,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 259167, 0, 3,
                                                                       251391, 194719, 252687,
                                                                       141967, 142912, 200767,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 260787, 0, 3,
                                                                       252687, 195727, 253983,
                                                                       142912, 143857, 202027,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 262407, 0, 3,
                                                                       253983, 196735, 255279,
                                                                       143857, 144802, 203287,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 264027, 0, 3,
                                                                       255279, 197743, 256575,
                                                                       144802, 145747, 204547,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 265647, 0, 3,
                                                                       256575, 198751, 257871,
                                                                       145747, 146692, 205807,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 267267, 0, 3,
                                                                       259167, 200767, 260787,
                                                                       148582, 149737, 207067,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 269247, 0, 3,
                                                                       260787, 202027, 262407,
                                                                       149737, 150892, 208607,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 271227, 0, 3,
                                                                       262407, 203287, 264027,
                                                                       150892, 152047, 210147,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 273207, 0, 3,
                                                                       264027, 204547, 265647,
                                                                       152047, 153202, 211687,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 275187, 0, 3,
                                                                       267267, 207067, 269247,
                                                                       155512, 156898, 213227,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 277563, 0, 3,
                                                                       269247, 208607, 271227,
                                                                       156898, 158284, 215075,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 279939, 0, 3,
                                                                       271227, 210147, 273207,
                                                                       158284, 159670, 216923,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 282315, 0, 3,
                                                                       275187, 213227, 277563,
                                                                       162442, 164080, 218771,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 285123, 0, 3,
                                                                       277563, 215075, 279939,
                                                                       164080, 165718, 220955,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsk_three_center_electron_repulsion_0(buffer, 287931, 0, 3,
                                                                       282315, 218771, 285123,
                                                                       168994, 170905, 223139,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 291207, 244335, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 292635, 251391, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 294471, 259167, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 296766, 267267, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 299571, 275187, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 302937, 282315, 2808, ncols);

                    simdfunc::contract_primitives(buffer, 306915, 287931, 3276, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 292215, 291207, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 293931, 292635, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 296091, 294471, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 298746, 296766, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 301947, 299571, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 305745, 302937, 78, 1, nmax);

        simdtrf::transform_k_inner(buffer, 310191, 306915, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 311556, 292215, 293931, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 312816, 293931, 296091, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 314436, 296091, 298746, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 316461, 298746, 301947, 15,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 318936, 301947, 305745, 15,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 321906, 305745, 310191, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 325416, 311556, 312816, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 327936, 312816, 314436, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 331176, 314436, 316461, 15,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 335226, 316461, 318936, 15,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 340176, 318936, 321906, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 346116, 325416, 327936, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 350316, 327936, 331176, 15,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 355716, 331176, 335226, 15,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 362466, 335226, 340176, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 370716, 346116, 350316, 15,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 377016, 350316, 355716, 15,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 385116, 355716, 362466, 15,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 395241, 370716, 377016, 15,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 404061, 377016, 385116, 15,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 415401, 395241, 404061, 15,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 427161, 415401, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 427161, 195, nmax);
    }

    for (size_t m = 0; m < 2535; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
