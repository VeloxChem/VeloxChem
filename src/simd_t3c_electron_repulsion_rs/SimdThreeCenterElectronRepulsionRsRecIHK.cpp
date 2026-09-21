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


#include "SimdThreeCenterElectronRepulsionRsRecIHK.hpp"

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
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ihk_three_center_electron_repulsion(double               *values,
                                               const size_t          npairs,
                                               const size_t          natoms,
                                               const CBasisFunction &a_function,
                                               const CBasisFunction &b_function,
                                               const CBasisFunction &c_function,
                                               const CSimdMatrix    &coordinates,
                                               const CSimdMatrix    &c_coordinates,
                                               CSimdMatrix          &buffer,
                                               const double          omega,
                                               const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_ihk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 609484, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 4290 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    simdfunc::compute_pair_exponents(a_function, b_function, coordinates, nmax);

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 609484, 444088, 30246, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

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

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 6, 3, {1, 2, 3,
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                                            15, 16, 17, 18}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 25, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17, 18}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 125, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 131, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 137, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 41, 42,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 42, 43,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 7, 8,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 8, 9,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 9, 10,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 164, 0, 3, 10, 11,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 170, 0, 3, 11, 12,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 176, 0, 3, 12, 13,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 13, 14,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 14, 15,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 194, 0, 3, 15, 16,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 200, 0, 3, 16, 17,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 206, 0, 3, 17, 18,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 18, 19,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 19, 20,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 224, 0, 3, 20, 21,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 230, 0, 3, 21, 22,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 236, 0, 3, 22, 23,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 26, 27,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 27, 28,
                                                                       98, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 254, 0, 3, 28, 29,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 260, 0, 3, 29, 30,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 266, 0, 3, 30, 31,
                                                                       107, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 31, 32,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 32, 33,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 284, 0, 3, 33, 34,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 290, 0, 3, 34, 35,
                                                                       119, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 296, 0, 3, 35, 36,
                                                                       122, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 36, 37,
                                                                       125, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 37, 38,
                                                                       128, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 314, 0, 3, 38, 39,
                                                                       131, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 320, 0, 3, 39, 40,
                                                                       134, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 326, 0, 3, 40, 41,
                                                                       137, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 332, 0, 3, 41, 42,
                                                                       140, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 44, 47,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 47, 50,
                                                                       152, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 50, 53,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 53, 56,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 56, 59,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 59, 62,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 62, 65,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 65, 68,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 68, 71,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 71, 74,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 74, 77,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 77, 80,
                                                                       212, 218, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 80, 83,
                                                                       218, 224, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 83, 86,
                                                                       224, 230, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 86, 89,
                                                                       230, 236, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 95, 98,
                                                                       242, 248, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 98,
                                                                       101, 248, 254, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 101,
                                                                       104, 254, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 104,
                                                                       107, 260, 266, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 107,
                                                                       110, 266, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 110,
                                                                       113, 272, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 113,
                                                                       116, 278, 284, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 116,
                                                                       119, 284, 290, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 119,
                                                                       122, 290, 296, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 122,
                                                                       125, 296, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 125,
                                                                       128, 302, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 598, 0, 3, 128,
                                                                       131, 308, 314, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 608, 0, 3, 131,
                                                                       134, 314, 320, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 618, 0, 3, 134,
                                                                       137, 320, 326, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 628, 0, 3, 137,
                                                                       140, 326, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 146,
                                                                       152, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 653, 0, 3, 152,
                                                                       158, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 668, 0, 3, 158,
                                                                       164, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 683, 0, 3, 164,
                                                                       170, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 698, 0, 3, 170,
                                                                       176, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 713, 0, 3, 176,
                                                                       182, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 728, 0, 3, 182,
                                                                       188, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 188,
                                                                       194, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 758, 0, 3, 194,
                                                                       200, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 773, 0, 3, 200,
                                                                       206, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 788, 0, 3, 206,
                                                                       212, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 803, 0, 3, 212,
                                                                       218, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 818, 0, 3, 218,
                                                                       224, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 224,
                                                                       230, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 848, 0, 3, 242,
                                                                       248, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 863, 0, 3, 248,
                                                                       254, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 878, 0, 3, 254,
                                                                       260, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 893, 0, 3, 260,
                                                                       266, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 908, 0, 3, 266,
                                                                       272, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 923, 0, 3, 272,
                                                                       278, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 938, 0, 3, 278,
                                                                       284, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 953, 0, 3, 284,
                                                                       290, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 968, 0, 3, 290,
                                                                       296, 568, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 983, 0, 3, 296,
                                                                       302, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 998, 0, 3, 302,
                                                                       308, 588, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1013, 0, 3, 308,
                                                                       314, 598, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1028, 0, 3, 314,
                                                                       320, 608, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1043, 0, 3, 320,
                                                                       326, 618, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 338,
                                                                       348, 638, 653, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1079, 0, 3, 348,
                                                                       358, 653, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 358,
                                                                       368, 668, 683, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1121, 0, 3, 368,
                                                                       378, 683, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 378,
                                                                       388, 698, 713, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 388,
                                                                       398, 713, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 398,
                                                                       408, 728, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1205, 0, 3, 408,
                                                                       418, 743, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 418,
                                                                       428, 758, 773, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1247, 0, 3, 428,
                                                                       438, 773, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 438,
                                                                       448, 788, 803, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1289, 0, 3, 448,
                                                                       458, 803, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 458,
                                                                       468, 818, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1331, 0, 3, 488,
                                                                       498, 848, 863, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 498,
                                                                       508, 863, 878, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1373, 0, 3, 508,
                                                                       518, 878, 893, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1394, 0, 3, 518,
                                                                       528, 893, 908, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1415, 0, 3, 528,
                                                                       538, 908, 923, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 538,
                                                                       548, 923, 938, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1457, 0, 3, 548,
                                                                       558, 938, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 558,
                                                                       568, 953, 968, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1499, 0, 3, 568,
                                                                       578, 968, 983, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 578,
                                                                       588, 983, 998, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1541, 0, 3, 588,
                                                                       598, 998, 1013, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1562, 0, 3, 598,
                                                                       608, 1013, 1028, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1583, 0, 3, 608,
                                                                       618, 1028, 1043, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 638,
                                                                       653, 1058, 1079, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 653,
                                                                       668, 1079, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 668,
                                                                       683, 1100, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 683,
                                                                       698, 1121, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 698,
                                                                       713, 1142, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 713,
                                                                       728, 1163, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 728,
                                                                       743, 1184, 1205, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 743,
                                                                       758, 1205, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 758,
                                                                       773, 1226, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 773,
                                                                       788, 1247, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 788,
                                                                       803, 1268, 1289, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 803,
                                                                       818, 1289, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 848,
                                                                       863, 1331, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 863,
                                                                       878, 1352, 1373, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 878,
                                                                       893, 1373, 1394, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 893,
                                                                       908, 1394, 1415, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 908,
                                                                       923, 1415, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 923,
                                                                       938, 1436, 1457, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 938,
                                                                       953, 1457, 1478, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2136, 0, 3, 953,
                                                                       968, 1478, 1499, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 968,
                                                                       983, 1499, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 983,
                                                                       998, 1520, 1541, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2220, 0, 3, 998,
                                                                       1013, 1541, 1562, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 1013,
                                                                       1028, 1562, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1058,
                                                                       1079, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 1079,
                                                                       1100, 1632, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1100,
                                                                       1121, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2384, 0, 3, 1121,
                                                                       1142, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2420, 0, 3, 1142,
                                                                       1163, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2456, 0, 3, 1163,
                                                                       1184, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2492, 0, 3, 1184,
                                                                       1205, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1205,
                                                                       1226, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2564, 0, 3, 1226,
                                                                       1247, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2600, 0, 3, 1247,
                                                                       1268, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2636, 0, 3, 1268,
                                                                       1289, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2672, 0, 3, 1331,
                                                                       1352, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2708, 0, 3, 1352,
                                                                       1373, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2744, 0, 3, 1373,
                                                                       1394, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2780, 0, 3, 1394,
                                                                       1415, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2816, 0, 3, 1415,
                                                                       1436, 2052, 2080, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2852, 0, 3, 1436,
                                                                       1457, 2080, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2888, 0, 3, 1457,
                                                                       1478, 2108, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2924, 0, 3, 1478,
                                                                       1499, 2136, 2164, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2960, 0, 3, 1499,
                                                                       1520, 2164, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2996, 0, 3, 1520,
                                                                       1541, 2192, 2220, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3032, 0, 3, 1541,
                                                                       1562, 2220, 2248, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3068, 0, 3, 1604,
                                                                       1632, 2276, 2312, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3113, 0, 3, 1632,
                                                                       1660, 2312, 2348, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3158, 0, 3, 1660,
                                                                       1688, 2348, 2384, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3203, 0, 3, 1688,
                                                                       1716, 2384, 2420, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3248, 0, 3, 1716,
                                                                       1744, 2420, 2456, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3293, 0, 3, 1744,
                                                                       1772, 2456, 2492, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3338, 0, 3, 1772,
                                                                       1800, 2492, 2528, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3383, 0, 3, 1800,
                                                                       1828, 2528, 2564, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3428, 0, 3, 1828,
                                                                       1856, 2564, 2600, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 1856,
                                                                       1884, 2600, 2636, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3518, 0, 3, 1940,
                                                                       1968, 2672, 2708, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3563, 0, 3, 1968,
                                                                       1996, 2708, 2744, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3608, 0, 3, 1996,
                                                                       2024, 2744, 2780, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3653, 0, 3, 2024,
                                                                       2052, 2780, 2816, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3698, 0, 3, 2052,
                                                                       2080, 2816, 2852, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3743, 0, 3, 2080,
                                                                       2108, 2852, 2888, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3788, 0, 3, 2108,
                                                                       2136, 2888, 2924, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3833, 0, 3, 2136,
                                                                       2164, 2924, 2960, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3878, 0, 3, 2164,
                                                                       2192, 2960, 2996, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3923, 0, 3, 2192,
                                                                       2220, 2996, 3032, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2276,
                                                                       2312, 3068, 3113, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 2312,
                                                                       2348, 3113, 3158, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 2348,
                                                                       2384, 3158, 3203, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4133, 0, 3, 2384,
                                                                       2420, 3203, 3248, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 2420,
                                                                       2456, 3248, 3293, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4243, 0, 3, 2456,
                                                                       2492, 3293, 3338, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 2492,
                                                                       2528, 3338, 3383, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4353, 0, 3, 2528,
                                                                       2564, 3383, 3428, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2564,
                                                                       2600, 3428, 3473, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4463, 0, 3, 2672,
                                                                       2708, 3518, 3563, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4518, 0, 3, 2708,
                                                                       2744, 3563, 3608, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4573, 0, 3, 2744,
                                                                       2780, 3608, 3653, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 2780,
                                                                       2816, 3653, 3698, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4683, 0, 3, 2816,
                                                                       2852, 3698, 3743, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4738, 0, 3, 2852,
                                                                       2888, 3743, 3788, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4793, 0, 3, 2888,
                                                                       2924, 3788, 3833, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 2924,
                                                                       2960, 3833, 3878, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4903, 0, 3, 2960,
                                                                       2996, 3878, 3923, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 3068,
                                                                       3113, 3968, 4023, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5024, 0, 3, 3113,
                                                                       3158, 4023, 4078, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5090, 0, 3, 3158,
                                                                       3203, 4078, 4133, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5156, 0, 3, 3203,
                                                                       3248, 4133, 4188, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5222, 0, 3, 3248,
                                                                       3293, 4188, 4243, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 3293,
                                                                       3338, 4243, 4298, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5354, 0, 3, 3338,
                                                                       3383, 4298, 4353, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 3383,
                                                                       3428, 4353, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5486, 0, 3, 3518,
                                                                       3563, 4463, 4518, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5552, 0, 3, 3563,
                                                                       3608, 4518, 4573, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 3608,
                                                                       3653, 4573, 4628, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5684, 0, 3, 3653,
                                                                       3698, 4628, 4683, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5750, 0, 3, 3698,
                                                                       3743, 4683, 4738, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5816, 0, 3, 3743,
                                                                       3788, 4738, 4793, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5882, 0, 3, 3788,
                                                                       3833, 4793, 4848, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 3833,
                                                                       3878, 4848, 4903, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6014, 0, 3, 3968,
                                                                       4023, 4958, 5024, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6092, 0, 3, 4023,
                                                                       4078, 5024, 5090, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6170, 0, 3, 4078,
                                                                       4133, 5090, 5156, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 4133,
                                                                       4188, 5156, 5222, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6326, 0, 3, 4188,
                                                                       4243, 5222, 5288, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6404, 0, 3, 4243,
                                                                       4298, 5288, 5354, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6482, 0, 3, 4298,
                                                                       4353, 5354, 5420, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6560, 0, 3, 4463,
                                                                       4518, 5486, 5552, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6638, 0, 3, 4518,
                                                                       4573, 5552, 5618, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 4573,
                                                                       4628, 5618, 5684, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6794, 0, 3, 4628,
                                                                       4683, 5684, 5750, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6872, 0, 3, 4683,
                                                                       4738, 5750, 5816, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6950, 0, 3, 4738,
                                                                       4793, 5816, 5882, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7028, 0, 3, 4793,
                                                                       4848, 5882, 5948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7106, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7109, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7112, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7115, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7118, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7121, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7124, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7127, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7130, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7133, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7136, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7139, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7142, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7145, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7148, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7151, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7154, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7157, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7160, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7163, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7166, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7169, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7172, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7175, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7178, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7181, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7184, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7187, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7190, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7193, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7196, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7199, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7202, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7205, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7208, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7211, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7214, 3, 7, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7223, 3, 8, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7232, 3, 9, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7241, 3, 10, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7250, 3, 11, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7259, 3, 12, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7268, 3, 13, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7277, 3, 14, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7286, 3, 15, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7295, 3, 16, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7304, 3, 17, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7313, 3, 18, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7322, 3, 19, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7331, 3, 20, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7340, 3, 21, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7349, 3, 22, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7358, 3, 23, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7367, 3, 26, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7376, 3, 27, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7385, 3, 28, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7394, 3, 29, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7403, 3, 30, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7412, 3, 31, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7421, 3, 32, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7430, 3, 33, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7439, 3, 34, 119,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7448, 3, 35, 122,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7457, 3, 36, 125,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7466, 3, 37, 128,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7475, 3, 38, 131,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7484, 3, 39, 134,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7493, 3, 40, 137,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7502, 3, 41, 140,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 7511, 3, 42, 143,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7520, 3, 44, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7538, 3, 47, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7556, 3, 50, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7574, 3, 53, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7592, 3, 56, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7610, 3, 59, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7628, 3, 62, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7646, 3, 65, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7664, 3, 68, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7682, 3, 71, 200,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7700, 3, 74, 206,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7718, 3, 77, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7736, 3, 80, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7754, 3, 83, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7772, 3, 86, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7790, 3, 89, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7808, 3, 95, 242,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7826, 3, 98, 248,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7844, 3, 101, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7862, 3, 104, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7880, 3, 107, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7898, 3, 110, 272,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7916, 3, 113, 278,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7934, 3, 116, 284,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7952, 3, 119, 290,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7970, 3, 122, 296,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 7988, 3, 125, 302,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8006, 3, 128, 308,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8024, 3, 131, 314,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8042, 3, 134, 320,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8060, 3, 137, 326,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 8078, 3, 140, 332,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8096, 3, 146, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8126, 3, 152, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8156, 3, 158, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8186, 3, 164, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8216, 3, 170, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8246, 3, 176, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8276, 3, 182, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8306, 3, 188, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8336, 3, 194, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8366, 3, 200, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8396, 3, 206, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8426, 3, 212, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8456, 3, 218, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8486, 3, 224, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8516, 3, 230, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8546, 3, 242, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8576, 3, 248, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8606, 3, 254, 508,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8636, 3, 260, 518,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8666, 3, 266, 528,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8696, 3, 272, 538,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8726, 3, 278, 548,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8756, 3, 284, 558,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8786, 3, 290, 568,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8816, 3, 296, 578,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8846, 3, 302, 588,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8876, 3, 308, 598,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8906, 3, 314, 608,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8936, 3, 320, 618,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 8966, 3, 326, 628,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 8996, 3, 338, 638,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9041, 3, 348, 653,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9086, 3, 358, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9131, 3, 368, 683,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9176, 3, 378, 698,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9221, 3, 388, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9266, 3, 398, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9311, 3, 408, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9356, 3, 418, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9401, 3, 428, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9446, 3, 438, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9491, 3, 448, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9536, 3, 458, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9581, 3, 468, 833,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9626, 3, 488, 848,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9671, 3, 498, 863,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9716, 3, 508, 878,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9761, 3, 518, 893,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9806, 3, 528, 908,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9851, 3, 538, 923,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9896, 3, 548, 938,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9941, 3, 558, 953,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 9986, 3, 568, 968,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10031, 3, 578,
                                                                       983, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10076, 3, 588,
                                                                       998, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10121, 3, 598,
                                                                       1013, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10166, 3, 608,
                                                                       1028, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10211, 3, 618,
                                                                       1043, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10256, 3, 638,
                                                                       1058, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10319, 3, 653,
                                                                       1079, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10382, 3, 668,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10445, 3, 683,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10508, 3, 698,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10571, 3, 713,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10634, 3, 728,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10697, 3, 743,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10760, 3, 758,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10823, 3, 773,
                                                                       1247, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10886, 3, 788,
                                                                       1268, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 10949, 3, 803,
                                                                       1289, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11012, 3, 818,
                                                                       1310, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11075, 3, 848,
                                                                       1331, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11138, 3, 863,
                                                                       1352, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11201, 3, 878,
                                                                       1373, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11264, 3, 893,
                                                                       1394, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11327, 3, 908,
                                                                       1415, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11390, 3, 923,
                                                                       1436, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11453, 3, 938,
                                                                       1457, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11516, 3, 953,
                                                                       1478, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11579, 3, 968,
                                                                       1499, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11642, 3, 983,
                                                                       1520, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11705, 3, 998,
                                                                       1541, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11768, 3, 1013,
                                                                       1562, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 11831, 3, 1028,
                                                                       1583, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11894, 3, 1058,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 11978, 3, 1079,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12062, 3, 1100,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12146, 3, 1121,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12230, 3, 1142,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12314, 3, 1163,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12398, 3, 1184,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12482, 3, 1205,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12566, 3, 1226,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12650, 3, 1247,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12734, 3, 1268,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12818, 3, 1289,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12902, 3, 1331,
                                                                       1940, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 12986, 3, 1352,
                                                                       1968, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13070, 3, 1373,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13154, 3, 1394,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13238, 3, 1415,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13322, 3, 1436,
                                                                       2080, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13406, 3, 1457,
                                                                       2108, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13490, 3, 1478,
                                                                       2136, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13574, 3, 1499,
                                                                       2164, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13658, 3, 1520,
                                                                       2192, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13742, 3, 1541,
                                                                       2220, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13826, 3, 1562,
                                                                       2248, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 13910, 3, 1604,
                                                                       2276, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14018, 3, 1632,
                                                                       2312, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14126, 3, 1660,
                                                                       2348, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14234, 3, 1688,
                                                                       2384, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14342, 3, 1716,
                                                                       2420, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14450, 3, 1744,
                                                                       2456, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14558, 3, 1772,
                                                                       2492, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14666, 3, 1800,
                                                                       2528, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14774, 3, 1828,
                                                                       2564, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14882, 3, 1856,
                                                                       2600, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 14990, 3, 1884,
                                                                       2636, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15098, 3, 1940,
                                                                       2672, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15206, 3, 1968,
                                                                       2708, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15314, 3, 1996,
                                                                       2744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15422, 3, 2024,
                                                                       2780, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15530, 3, 2052,
                                                                       2816, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15638, 3, 2080,
                                                                       2852, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15746, 3, 2108,
                                                                       2888, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15854, 3, 2136,
                                                                       2924, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15962, 3, 2164,
                                                                       2960, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16070, 3, 2192,
                                                                       2996, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16178, 3, 2220,
                                                                       3032, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16286, 3, 2276,
                                                                       3068, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16421, 3, 2312,
                                                                       3113, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16556, 3, 2348,
                                                                       3158, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16691, 3, 2384,
                                                                       3203, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16826, 3, 2420,
                                                                       3248, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 16961, 3, 2456,
                                                                       3293, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17096, 3, 2492,
                                                                       3338, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17231, 3, 2528,
                                                                       3383, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17366, 3, 2564,
                                                                       3428, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17501, 3, 2600,
                                                                       3473, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17636, 3, 2672,
                                                                       3518, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17771, 3, 2708,
                                                                       3563, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17906, 3, 2744,
                                                                       3608, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18041, 3, 2780,
                                                                       3653, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18176, 3, 2816,
                                                                       3698, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18311, 3, 2852,
                                                                       3743, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18446, 3, 2888,
                                                                       3788, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18581, 3, 2924,
                                                                       3833, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18716, 3, 2960,
                                                                       3878, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18851, 3, 2996,
                                                                       3923, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 18986, 3, 3068,
                                                                       3968, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19151, 3, 3113,
                                                                       4023, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19316, 3, 3158,
                                                                       4078, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19481, 3, 3203,
                                                                       4133, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19646, 3, 3248,
                                                                       4188, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19811, 3, 3293,
                                                                       4243, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19976, 3, 3338,
                                                                       4298, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20141, 3, 3383,
                                                                       4353, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20306, 3, 3428,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20471, 3, 3518,
                                                                       4463, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20636, 3, 3563,
                                                                       4518, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20801, 3, 3608,
                                                                       4573, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20966, 3, 3653,
                                                                       4628, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21131, 3, 3698,
                                                                       4683, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21296, 3, 3743,
                                                                       4738, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21461, 3, 3788,
                                                                       4793, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21626, 3, 3833,
                                                                       4848, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21791, 3, 3878,
                                                                       4903, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 21956, 3, 3968,
                                                                       4958, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22154, 3, 4023,
                                                                       5024, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22352, 3, 4078,
                                                                       5090, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22550, 3, 4133,
                                                                       5156, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22748, 3, 4188,
                                                                       5222, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22946, 3, 4243,
                                                                       5288, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23144, 3, 4298,
                                                                       5354, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23342, 3, 4353,
                                                                       5420, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23540, 3, 4463,
                                                                       5486, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23738, 3, 4518,
                                                                       5552, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23936, 3, 4573,
                                                                       5618, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24134, 3, 4628,
                                                                       5684, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24332, 3, 4683,
                                                                       5750, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24530, 3, 4738,
                                                                       5816, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24728, 3, 4793,
                                                                       5882, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24926, 3, 4848,
                                                                       5948, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25124, 3, 4958,
                                                                       6014, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25358, 3, 5024,
                                                                       6092, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25592, 3, 5090,
                                                                       6170, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25826, 3, 5156,
                                                                       6248, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26060, 3, 5222,
                                                                       6326, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26294, 3, 5288,
                                                                       6404, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26528, 3, 5354,
                                                                       6482, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26762, 3, 5486,
                                                                       6560, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26996, 3, 5552,
                                                                       6638, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27230, 3, 5618,
                                                                       6716, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27464, 3, 5684,
                                                                       6794, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27698, 3, 5750,
                                                                       6872, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27932, 3, 5816,
                                                                       6950, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 28166, 3, 5882,
                                                                       7028, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28400, 3, 7, 8,
                                                                       7112, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28406, 3, 8, 9,
                                                                       7115, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28412, 3, 9, 10,
                                                                       7118, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28418, 3, 10, 11,
                                                                       7121, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28424, 3, 11, 12,
                                                                       7124, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28430, 3, 12, 13,
                                                                       7127, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28436, 3, 13, 14,
                                                                       7130, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28442, 3, 14, 15,
                                                                       7133, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28448, 3, 15, 16,
                                                                       7136, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28454, 3, 16, 17,
                                                                       7139, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28460, 3, 17, 18,
                                                                       7142, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28466, 3, 18, 19,
                                                                       7145, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28472, 3, 19, 20,
                                                                       7148, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28478, 3, 20, 21,
                                                                       7151, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28484, 3, 21, 22,
                                                                       7154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28490, 3, 22, 23,
                                                                       7157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28496, 3, 26, 27,
                                                                       7166, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28502, 3, 27, 28,
                                                                       7169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28508, 3, 28, 29,
                                                                       7172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28514, 3, 29, 30,
                                                                       7175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28520, 3, 30, 31,
                                                                       7178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28526, 3, 31, 32,
                                                                       7181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28532, 3, 32, 33,
                                                                       7184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28538, 3, 33, 34,
                                                                       7187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28544, 3, 34, 35,
                                                                       7190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28550, 3, 35, 36,
                                                                       7193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28556, 3, 36, 37,
                                                                       7196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28562, 3, 37, 38,
                                                                       7199, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28568, 3, 38, 39,
                                                                       7202, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28574, 3, 39, 40,
                                                                       7205, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28580, 3, 40, 41,
                                                                       7208, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28586, 3, 41, 42,
                                                                       7211, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28592, 0, 3,
                                                                       28400, 7112, 28406, 7232,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28610, 0, 3,
                                                                       28406, 7115, 28412, 7241,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28628, 0, 3,
                                                                       28412, 7118, 28418, 7250,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28646, 0, 3,
                                                                       28418, 7121, 28424, 7259,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28664, 0, 3,
                                                                       28424, 7124, 28430, 7268,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28682, 0, 3,
                                                                       28430, 7127, 28436, 7277,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28700, 0, 3,
                                                                       28436, 7130, 28442, 7286,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28718, 0, 3,
                                                                       28442, 7133, 28448, 7295,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28736, 0, 3,
                                                                       28448, 7136, 28454, 7304,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28754, 0, 3,
                                                                       28454, 7139, 28460, 7313,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28772, 0, 3,
                                                                       28460, 7142, 28466, 7322,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28790, 0, 3,
                                                                       28466, 7145, 28472, 7331,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28808, 0, 3,
                                                                       28472, 7148, 28478, 7340,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28826, 0, 3,
                                                                       28478, 7151, 28484, 7349,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28844, 0, 3,
                                                                       28484, 7154, 28490, 7358,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28862, 0, 3,
                                                                       28496, 7166, 28502, 7385,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28880, 0, 3,
                                                                       28502, 7169, 28508, 7394,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28898, 0, 3,
                                                                       28508, 7172, 28514, 7403,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28916, 0, 3,
                                                                       28514, 7175, 28520, 7412,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28934, 0, 3,
                                                                       28520, 7178, 28526, 7421,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28952, 0, 3,
                                                                       28526, 7181, 28532, 7430,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28970, 0, 3,
                                                                       28532, 7184, 28538, 7439,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 28988, 0, 3,
                                                                       28538, 7187, 28544, 7448,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29006, 0, 3,
                                                                       28544, 7190, 28550, 7457,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29024, 0, 3,
                                                                       28550, 7193, 28556, 7466,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29042, 0, 3,
                                                                       28556, 7196, 28562, 7475,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29060, 0, 3,
                                                                       28562, 7199, 28568, 7484,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29078, 0, 3,
                                                                       28568, 7202, 28574, 7493,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29096, 0, 3,
                                                                       28574, 7205, 28580, 7502,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 29114, 0, 3,
                                                                       28580, 7208, 28586, 7511,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29132, 0, 3,
                                                                       28592, 7232, 28610, 146,
                                                                       152, 7556, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29168, 0, 3,
                                                                       28610, 7241, 28628, 152,
                                                                       158, 7574, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29204, 0, 3,
                                                                       28628, 7250, 28646, 158,
                                                                       164, 7592, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29240, 0, 3,
                                                                       28646, 7259, 28664, 164,
                                                                       170, 7610, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29276, 0, 3,
                                                                       28664, 7268, 28682, 170,
                                                                       176, 7628, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29312, 0, 3,
                                                                       28682, 7277, 28700, 176,
                                                                       182, 7646, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29348, 0, 3,
                                                                       28700, 7286, 28718, 182,
                                                                       188, 7664, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29384, 0, 3,
                                                                       28718, 7295, 28736, 188,
                                                                       194, 7682, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29420, 0, 3,
                                                                       28736, 7304, 28754, 194,
                                                                       200, 7700, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29456, 0, 3,
                                                                       28754, 7313, 28772, 200,
                                                                       206, 7718, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29492, 0, 3,
                                                                       28772, 7322, 28790, 206,
                                                                       212, 7736, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29528, 0, 3,
                                                                       28790, 7331, 28808, 212,
                                                                       218, 7754, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29564, 0, 3,
                                                                       28808, 7340, 28826, 218,
                                                                       224, 7772, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29600, 0, 3,
                                                                       28826, 7349, 28844, 224,
                                                                       230, 7790, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29636, 0, 3,
                                                                       28862, 7385, 28880, 242,
                                                                       248, 7844, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29672, 0, 3,
                                                                       28880, 7394, 28898, 248,
                                                                       254, 7862, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29708, 0, 3,
                                                                       28898, 7403, 28916, 254,
                                                                       260, 7880, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29744, 0, 3,
                                                                       28916, 7412, 28934, 260,
                                                                       266, 7898, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29780, 0, 3,
                                                                       28934, 7421, 28952, 266,
                                                                       272, 7916, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29816, 0, 3,
                                                                       28952, 7430, 28970, 272,
                                                                       278, 7934, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29852, 0, 3,
                                                                       28970, 7439, 28988, 278,
                                                                       284, 7952, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29888, 0, 3,
                                                                       28988, 7448, 29006, 284,
                                                                       290, 7970, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29924, 0, 3,
                                                                       29006, 7457, 29024, 290,
                                                                       296, 7988, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29960, 0, 3,
                                                                       29024, 7466, 29042, 296,
                                                                       302, 8006, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 29996, 0, 3,
                                                                       29042, 7475, 29060, 302,
                                                                       308, 8024, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30032, 0, 3,
                                                                       29060, 7484, 29078, 308,
                                                                       314, 8042, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30068, 0, 3,
                                                                       29078, 7493, 29096, 314,
                                                                       320, 8060, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 30104, 0, 3,
                                                                       29096, 7502, 29114, 320,
                                                                       326, 8078, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30140, 0, 3,
                                                                       29132, 7556, 29168, 338,
                                                                       348, 8156, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30200, 0, 3,
                                                                       29168, 7574, 29204, 348,
                                                                       358, 8186, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30260, 0, 3,
                                                                       29204, 7592, 29240, 358,
                                                                       368, 8216, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30320, 0, 3,
                                                                       29240, 7610, 29276, 368,
                                                                       378, 8246, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30380, 0, 3,
                                                                       29276, 7628, 29312, 378,
                                                                       388, 8276, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30440, 0, 3,
                                                                       29312, 7646, 29348, 388,
                                                                       398, 8306, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30500, 0, 3,
                                                                       29348, 7664, 29384, 398,
                                                                       408, 8336, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30560, 0, 3,
                                                                       29384, 7682, 29420, 408,
                                                                       418, 8366, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30620, 0, 3,
                                                                       29420, 7700, 29456, 418,
                                                                       428, 8396, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30680, 0, 3,
                                                                       29456, 7718, 29492, 428,
                                                                       438, 8426, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30740, 0, 3,
                                                                       29492, 7736, 29528, 438,
                                                                       448, 8456, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30800, 0, 3,
                                                                       29528, 7754, 29564, 448,
                                                                       458, 8486, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30860, 0, 3,
                                                                       29564, 7772, 29600, 458,
                                                                       468, 8516, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30920, 0, 3,
                                                                       29636, 7844, 29672, 488,
                                                                       498, 8606, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 30980, 0, 3,
                                                                       29672, 7862, 29708, 498,
                                                                       508, 8636, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31040, 0, 3,
                                                                       29708, 7880, 29744, 508,
                                                                       518, 8666, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31100, 0, 3,
                                                                       29744, 7898, 29780, 518,
                                                                       528, 8696, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31160, 0, 3,
                                                                       29780, 7916, 29816, 528,
                                                                       538, 8726, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31220, 0, 3,
                                                                       29816, 7934, 29852, 538,
                                                                       548, 8756, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31280, 0, 3,
                                                                       29852, 7952, 29888, 548,
                                                                       558, 8786, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31340, 0, 3,
                                                                       29888, 7970, 29924, 558,
                                                                       568, 8816, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31400, 0, 3,
                                                                       29924, 7988, 29960, 568,
                                                                       578, 8846, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31460, 0, 3,
                                                                       29960, 8006, 29996, 578,
                                                                       588, 8876, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31520, 0, 3,
                                                                       29996, 8024, 30032, 588,
                                                                       598, 8906, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31580, 0, 3,
                                                                       30032, 8042, 30068, 598,
                                                                       608, 8936, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 31640, 0, 3,
                                                                       30068, 8060, 30104, 608,
                                                                       618, 8966, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 31700, 0, 3,
                                                                       30140, 8156, 30200, 638,
                                                                       653, 9086, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 31790, 0, 3,
                                                                       30200, 8186, 30260, 653,
                                                                       668, 9131, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 31880, 0, 3,
                                                                       30260, 8216, 30320, 668,
                                                                       683, 9176, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 31970, 0, 3,
                                                                       30320, 8246, 30380, 683,
                                                                       698, 9221, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32060, 0, 3,
                                                                       30380, 8276, 30440, 698,
                                                                       713, 9266, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32150, 0, 3,
                                                                       30440, 8306, 30500, 713,
                                                                       728, 9311, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32240, 0, 3,
                                                                       30500, 8336, 30560, 728,
                                                                       743, 9356, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32330, 0, 3,
                                                                       30560, 8366, 30620, 743,
                                                                       758, 9401, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32420, 0, 3,
                                                                       30620, 8396, 30680, 758,
                                                                       773, 9446, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32510, 0, 3,
                                                                       30680, 8426, 30740, 773,
                                                                       788, 9491, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32600, 0, 3,
                                                                       30740, 8456, 30800, 788,
                                                                       803, 9536, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32690, 0, 3,
                                                                       30800, 8486, 30860, 803,
                                                                       818, 9581, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32780, 0, 3,
                                                                       30920, 8606, 30980, 848,
                                                                       863, 9716, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32870, 0, 3,
                                                                       30980, 8636, 31040, 863,
                                                                       878, 9761, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 32960, 0, 3,
                                                                       31040, 8666, 31100, 878,
                                                                       893, 9806, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33050, 0, 3,
                                                                       31100, 8696, 31160, 893,
                                                                       908, 9851, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33140, 0, 3,
                                                                       31160, 8726, 31220, 908,
                                                                       923, 9896, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33230, 0, 3,
                                                                       31220, 8756, 31280, 923,
                                                                       938, 9941, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33320, 0, 3,
                                                                       31280, 8786, 31340, 938,
                                                                       953, 9986, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33410, 0, 3,
                                                                       31340, 8816, 31400, 953,
                                                                       968, 10031, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33500, 0, 3,
                                                                       31400, 8846, 31460, 968,
                                                                       983, 10076, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33590, 0, 3,
                                                                       31460, 8876, 31520, 983,
                                                                       998, 10121, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33680, 0, 3,
                                                                       31520, 8906, 31580, 998,
                                                                       1013, 10166, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 33770, 0, 3,
                                                                       31580, 8936, 31640, 1013,
                                                                       1028, 10211, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 33860, 0, 3,
                                                                       31700, 9086, 31790, 1058,
                                                                       1079, 10382, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 33986, 0, 3,
                                                                       31790, 9131, 31880, 1079,
                                                                       1100, 10445, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34112, 0, 3,
                                                                       31880, 9176, 31970, 1100,
                                                                       1121, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34238, 0, 3,
                                                                       31970, 9221, 32060, 1121,
                                                                       1142, 10571, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34364, 0, 3,
                                                                       32060, 9266, 32150, 1142,
                                                                       1163, 10634, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34490, 0, 3,
                                                                       32150, 9311, 32240, 1163,
                                                                       1184, 10697, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34616, 0, 3,
                                                                       32240, 9356, 32330, 1184,
                                                                       1205, 10760, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34742, 0, 3,
                                                                       32330, 9401, 32420, 1205,
                                                                       1226, 10823, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34868, 0, 3,
                                                                       32420, 9446, 32510, 1226,
                                                                       1247, 10886, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 34994, 0, 3,
                                                                       32510, 9491, 32600, 1247,
                                                                       1268, 10949, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35120, 0, 3,
                                                                       32600, 9536, 32690, 1268,
                                                                       1289, 11012, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35246, 0, 3,
                                                                       32780, 9716, 32870, 1331,
                                                                       1352, 11201, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35372, 0, 3,
                                                                       32870, 9761, 32960, 1352,
                                                                       1373, 11264, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35498, 0, 3,
                                                                       32960, 9806, 33050, 1373,
                                                                       1394, 11327, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35624, 0, 3,
                                                                       33050, 9851, 33140, 1394,
                                                                       1415, 11390, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35750, 0, 3,
                                                                       33140, 9896, 33230, 1415,
                                                                       1436, 11453, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 35876, 0, 3,
                                                                       33230, 9941, 33320, 1436,
                                                                       1457, 11516, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36002, 0, 3,
                                                                       33320, 9986, 33410, 1457,
                                                                       1478, 11579, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36128, 0, 3,
                                                                       33410, 10031, 33500, 1478,
                                                                       1499, 11642, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36254, 0, 3,
                                                                       33500, 10076, 33590, 1499,
                                                                       1520, 11705, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36380, 0, 3,
                                                                       33590, 10121, 33680, 1520,
                                                                       1541, 11768, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36506, 0, 3,
                                                                       33680, 10166, 33770, 1541,
                                                                       1562, 11831, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 36632, 0, 3,
                                                                       33860, 10382, 33986, 1604,
                                                                       1632, 12062, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 36800, 0, 3,
                                                                       33986, 10445, 34112, 1632,
                                                                       1660, 12146, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 36968, 0, 3,
                                                                       34112, 10508, 34238, 1660,
                                                                       1688, 12230, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37136, 0, 3,
                                                                       34238, 10571, 34364, 1688,
                                                                       1716, 12314, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37304, 0, 3,
                                                                       34364, 10634, 34490, 1716,
                                                                       1744, 12398, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37472, 0, 3,
                                                                       34490, 10697, 34616, 1744,
                                                                       1772, 12482, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37640, 0, 3,
                                                                       34616, 10760, 34742, 1772,
                                                                       1800, 12566, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37808, 0, 3,
                                                                       34742, 10823, 34868, 1800,
                                                                       1828, 12650, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 37976, 0, 3,
                                                                       34868, 10886, 34994, 1828,
                                                                       1856, 12734, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38144, 0, 3,
                                                                       34994, 10949, 35120, 1856,
                                                                       1884, 12818, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38312, 0, 3,
                                                                       35246, 11201, 35372, 1940,
                                                                       1968, 13070, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38480, 0, 3,
                                                                       35372, 11264, 35498, 1968,
                                                                       1996, 13154, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38648, 0, 3,
                                                                       35498, 11327, 35624, 1996,
                                                                       2024, 13238, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38816, 0, 3,
                                                                       35624, 11390, 35750, 2024,
                                                                       2052, 13322, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 38984, 0, 3,
                                                                       35750, 11453, 35876, 2052,
                                                                       2080, 13406, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39152, 0, 3,
                                                                       35876, 11516, 36002, 2080,
                                                                       2108, 13490, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39320, 0, 3,
                                                                       36002, 11579, 36128, 2108,
                                                                       2136, 13574, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39488, 0, 3,
                                                                       36128, 11642, 36254, 2136,
                                                                       2164, 13658, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39656, 0, 3,
                                                                       36254, 11705, 36380, 2164,
                                                                       2192, 13742, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39824, 0, 3,
                                                                       36380, 11768, 36506, 2192,
                                                                       2220, 13826, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 39992, 0, 3,
                                                                       36632, 12062, 36800, 2276,
                                                                       2312, 14126, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 40208, 0, 3,
                                                                       36800, 12146, 36968, 2312,
                                                                       2348, 14234, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 40424, 0, 3,
                                                                       36968, 12230, 37136, 2348,
                                                                       2384, 14342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 40640, 0, 3,
                                                                       37136, 12314, 37304, 2384,
                                                                       2420, 14450, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 40856, 0, 3,
                                                                       37304, 12398, 37472, 2420,
                                                                       2456, 14558, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41072, 0, 3,
                                                                       37472, 12482, 37640, 2456,
                                                                       2492, 14666, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41288, 0, 3,
                                                                       37640, 12566, 37808, 2492,
                                                                       2528, 14774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41504, 0, 3,
                                                                       37808, 12650, 37976, 2528,
                                                                       2564, 14882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41720, 0, 3,
                                                                       37976, 12734, 38144, 2564,
                                                                       2600, 14990, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 41936, 0, 3,
                                                                       38312, 13070, 38480, 2672,
                                                                       2708, 15314, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42152, 0, 3,
                                                                       38480, 13154, 38648, 2708,
                                                                       2744, 15422, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42368, 0, 3,
                                                                       38648, 13238, 38816, 2744,
                                                                       2780, 15530, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42584, 0, 3,
                                                                       38816, 13322, 38984, 2780,
                                                                       2816, 15638, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 42800, 0, 3,
                                                                       38984, 13406, 39152, 2816,
                                                                       2852, 15746, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 43016, 0, 3,
                                                                       39152, 13490, 39320, 2852,
                                                                       2888, 15854, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 43232, 0, 3,
                                                                       39320, 13574, 39488, 2888,
                                                                       2924, 15962, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 43448, 0, 3,
                                                                       39488, 13658, 39656, 2924,
                                                                       2960, 16070, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 43664, 0, 3,
                                                                       39656, 13742, 39824, 2960,
                                                                       2996, 16178, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 43880, 0, 3,
                                                                       39992, 14126, 40208, 3068,
                                                                       3113, 16556, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 44150, 0, 3,
                                                                       40208, 14234, 40424, 3113,
                                                                       3158, 16691, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 44420, 0, 3,
                                                                       40424, 14342, 40640, 3158,
                                                                       3203, 16826, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 44690, 0, 3,
                                                                       40640, 14450, 40856, 3203,
                                                                       3248, 16961, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 44960, 0, 3,
                                                                       40856, 14558, 41072, 3248,
                                                                       3293, 17096, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 45230, 0, 3,
                                                                       41072, 14666, 41288, 3293,
                                                                       3338, 17231, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 45500, 0, 3,
                                                                       41288, 14774, 41504, 3338,
                                                                       3383, 17366, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 45770, 0, 3,
                                                                       41504, 14882, 41720, 3383,
                                                                       3428, 17501, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 46040, 0, 3,
                                                                       41936, 15314, 42152, 3518,
                                                                       3563, 17906, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 46310, 0, 3,
                                                                       42152, 15422, 42368, 3563,
                                                                       3608, 18041, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 46580, 0, 3,
                                                                       42368, 15530, 42584, 3608,
                                                                       3653, 18176, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 46850, 0, 3,
                                                                       42584, 15638, 42800, 3653,
                                                                       3698, 18311, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 47120, 0, 3,
                                                                       42800, 15746, 43016, 3698,
                                                                       3743, 18446, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 47390, 0, 3,
                                                                       43016, 15854, 43232, 3743,
                                                                       3788, 18581, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 47660, 0, 3,
                                                                       43232, 15962, 43448, 3788,
                                                                       3833, 18716, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 47930, 0, 3,
                                                                       43448, 16070, 43664, 3833,
                                                                       3878, 18851, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 48200, 0, 3,
                                                                       43880, 16556, 44150, 3968,
                                                                       4023, 19316, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 48530, 0, 3,
                                                                       44150, 16691, 44420, 4023,
                                                                       4078, 19481, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 48860, 0, 3,
                                                                       44420, 16826, 44690, 4078,
                                                                       4133, 19646, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 49190, 0, 3,
                                                                       44690, 16961, 44960, 4133,
                                                                       4188, 19811, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 49520, 0, 3,
                                                                       44960, 17096, 45230, 4188,
                                                                       4243, 19976, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 49850, 0, 3,
                                                                       45230, 17231, 45500, 4243,
                                                                       4298, 20141, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 50180, 0, 3,
                                                                       45500, 17366, 45770, 4298,
                                                                       4353, 20306, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 50510, 0, 3,
                                                                       46040, 17906, 46310, 4463,
                                                                       4518, 20801, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 50840, 0, 3,
                                                                       46310, 18041, 46580, 4518,
                                                                       4573, 20966, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 51170, 0, 3,
                                                                       46580, 18176, 46850, 4573,
                                                                       4628, 21131, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 51500, 0, 3,
                                                                       46850, 18311, 47120, 4628,
                                                                       4683, 21296, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 51830, 0, 3,
                                                                       47120, 18446, 47390, 4683,
                                                                       4738, 21461, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 52160, 0, 3,
                                                                       47390, 18581, 47660, 4738,
                                                                       4793, 21626, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 52490, 0, 3,
                                                                       47660, 18716, 47930, 4793,
                                                                       4848, 21791, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 52820, 0, 3,
                                                                       48200, 19316, 48530, 4958,
                                                                       5024, 22352, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 53216, 0, 3,
                                                                       48530, 19481, 48860, 5024,
                                                                       5090, 22550, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 53612, 0, 3,
                                                                       48860, 19646, 49190, 5090,
                                                                       5156, 22748, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 54008, 0, 3,
                                                                       49190, 19811, 49520, 5156,
                                                                       5222, 22946, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 54404, 0, 3,
                                                                       49520, 19976, 49850, 5222,
                                                                       5288, 23144, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 54800, 0, 3,
                                                                       49850, 20141, 50180, 5288,
                                                                       5354, 23342, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 55196, 0, 3,
                                                                       50510, 20801, 50840, 5486,
                                                                       5552, 23936, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 55592, 0, 3,
                                                                       50840, 20966, 51170, 5552,
                                                                       5618, 24134, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 55988, 0, 3,
                                                                       51170, 21131, 51500, 5618,
                                                                       5684, 24332, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 56384, 0, 3,
                                                                       51500, 21296, 51830, 5684,
                                                                       5750, 24530, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 56780, 0, 3,
                                                                       51830, 21461, 52160, 5750,
                                                                       5816, 24728, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 57176, 0, 3,
                                                                       52160, 21626, 52490, 5816,
                                                                       5882, 24926, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 57572, 0, 3,
                                                                       52820, 22352, 53216, 6014,
                                                                       6092, 25592, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 58040, 0, 3,
                                                                       53216, 22550, 53612, 6092,
                                                                       6170, 25826, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 58508, 0, 3,
                                                                       53612, 22748, 54008, 6170,
                                                                       6248, 26060, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 58976, 0, 3,
                                                                       54008, 22946, 54404, 6248,
                                                                       6326, 26294, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 59444, 0, 3,
                                                                       54404, 23144, 54800, 6326,
                                                                       6404, 26528, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 59912, 0, 3,
                                                                       55196, 23936, 55592, 6560,
                                                                       6638, 27230, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 60380, 0, 3,
                                                                       55592, 24134, 55988, 6638,
                                                                       6716, 27464, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 60848, 0, 3,
                                                                       55988, 24332, 56384, 6716,
                                                                       6794, 27698, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 61316, 0, 3,
                                                                       56384, 24530, 56780, 6794,
                                                                       6872, 27932, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 61784, 0, 3,
                                                                       56780, 24728, 57176, 6872,
                                                                       6950, 28166, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62252, 3, 7106,
                                                                       7109, 28400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62262, 3, 7109,
                                                                       7112, 28406, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62272, 3, 7112,
                                                                       7115, 28412, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62282, 3, 7115,
                                                                       7118, 28418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62292, 3, 7118,
                                                                       7121, 28424, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62302, 3, 7121,
                                                                       7124, 28430, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62312, 3, 7124,
                                                                       7127, 28436, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62322, 3, 7127,
                                                                       7130, 28442, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62332, 3, 7130,
                                                                       7133, 28448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62342, 3, 7133,
                                                                       7136, 28454, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62352, 3, 7136,
                                                                       7139, 28460, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62362, 3, 7139,
                                                                       7142, 28466, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62372, 3, 7142,
                                                                       7145, 28472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62382, 3, 7145,
                                                                       7148, 28478, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62392, 3, 7148,
                                                                       7151, 28484, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62402, 3, 7151,
                                                                       7154, 28490, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62412, 3, 7160,
                                                                       7163, 28496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62422, 3, 7163,
                                                                       7166, 28502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62432, 3, 7166,
                                                                       7169, 28508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62442, 3, 7169,
                                                                       7172, 28514, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62452, 3, 7172,
                                                                       7175, 28520, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62462, 3, 7175,
                                                                       7178, 28526, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62472, 3, 7178,
                                                                       7181, 28532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62482, 3, 7181,
                                                                       7184, 28538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62492, 3, 7184,
                                                                       7187, 28544, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62502, 3, 7187,
                                                                       7190, 28550, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62512, 3, 7190,
                                                                       7193, 28556, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62522, 3, 7193,
                                                                       7196, 28562, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62532, 3, 7196,
                                                                       7199, 28568, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62542, 3, 7199,
                                                                       7202, 28574, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62552, 3, 7202,
                                                                       7205, 28580, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62562, 3, 7205,
                                                                       7208, 28586, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62572, 0, 3,
                                                                       62252, 28400, 62262, 7214,
                                                                       7223, 28592, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62602, 0, 3,
                                                                       62262, 28406, 62272, 7223,
                                                                       7232, 28610, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62632, 0, 3,
                                                                       62272, 28412, 62282, 7232,
                                                                       7241, 28628, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62662, 0, 3,
                                                                       62282, 28418, 62292, 7241,
                                                                       7250, 28646, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62692, 0, 3,
                                                                       62292, 28424, 62302, 7250,
                                                                       7259, 28664, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62722, 0, 3,
                                                                       62302, 28430, 62312, 7259,
                                                                       7268, 28682, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62752, 0, 3,
                                                                       62312, 28436, 62322, 7268,
                                                                       7277, 28700, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62782, 0, 3,
                                                                       62322, 28442, 62332, 7277,
                                                                       7286, 28718, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62812, 0, 3,
                                                                       62332, 28448, 62342, 7286,
                                                                       7295, 28736, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62842, 0, 3,
                                                                       62342, 28454, 62352, 7295,
                                                                       7304, 28754, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62872, 0, 3,
                                                                       62352, 28460, 62362, 7304,
                                                                       7313, 28772, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62902, 0, 3,
                                                                       62362, 28466, 62372, 7313,
                                                                       7322, 28790, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62932, 0, 3,
                                                                       62372, 28472, 62382, 7322,
                                                                       7331, 28808, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62962, 0, 3,
                                                                       62382, 28478, 62392, 7331,
                                                                       7340, 28826, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 62992, 0, 3,
                                                                       62392, 28484, 62402, 7340,
                                                                       7349, 28844, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63022, 0, 3,
                                                                       62412, 28496, 62422, 7367,
                                                                       7376, 28862, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63052, 0, 3,
                                                                       62422, 28502, 62432, 7376,
                                                                       7385, 28880, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63082, 0, 3,
                                                                       62432, 28508, 62442, 7385,
                                                                       7394, 28898, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63112, 0, 3,
                                                                       62442, 28514, 62452, 7394,
                                                                       7403, 28916, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63142, 0, 3,
                                                                       62452, 28520, 62462, 7403,
                                                                       7412, 28934, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63172, 0, 3,
                                                                       62462, 28526, 62472, 7412,
                                                                       7421, 28952, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63202, 0, 3,
                                                                       62472, 28532, 62482, 7421,
                                                                       7430, 28970, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63232, 0, 3,
                                                                       62482, 28538, 62492, 7430,
                                                                       7439, 28988, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63262, 0, 3,
                                                                       62492, 28544, 62502, 7439,
                                                                       7448, 29006, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63292, 0, 3,
                                                                       62502, 28550, 62512, 7448,
                                                                       7457, 29024, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63322, 0, 3,
                                                                       62512, 28556, 62522, 7457,
                                                                       7466, 29042, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63352, 0, 3,
                                                                       62522, 28562, 62532, 7466,
                                                                       7475, 29060, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63382, 0, 3,
                                                                       62532, 28568, 62542, 7475,
                                                                       7484, 29078, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63412, 0, 3,
                                                                       62542, 28574, 62552, 7484,
                                                                       7493, 29096, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 63442, 0, 3,
                                                                       62552, 28580, 62562, 7493,
                                                                       7502, 29114, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63472, 0, 3,
                                                                       62572, 28592, 62602, 7520,
                                                                       7538, 29132, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63532, 0, 3,
                                                                       62602, 28610, 62632, 7538,
                                                                       7556, 29168, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63592, 0, 3,
                                                                       62632, 28628, 62662, 7556,
                                                                       7574, 29204, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63652, 0, 3,
                                                                       62662, 28646, 62692, 7574,
                                                                       7592, 29240, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63712, 0, 3,
                                                                       62692, 28664, 62722, 7592,
                                                                       7610, 29276, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63772, 0, 3,
                                                                       62722, 28682, 62752, 7610,
                                                                       7628, 29312, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63832, 0, 3,
                                                                       62752, 28700, 62782, 7628,
                                                                       7646, 29348, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63892, 0, 3,
                                                                       62782, 28718, 62812, 7646,
                                                                       7664, 29384, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 63952, 0, 3,
                                                                       62812, 28736, 62842, 7664,
                                                                       7682, 29420, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64012, 0, 3,
                                                                       62842, 28754, 62872, 7682,
                                                                       7700, 29456, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64072, 0, 3,
                                                                       62872, 28772, 62902, 7700,
                                                                       7718, 29492, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64132, 0, 3,
                                                                       62902, 28790, 62932, 7718,
                                                                       7736, 29528, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64192, 0, 3,
                                                                       62932, 28808, 62962, 7736,
                                                                       7754, 29564, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64252, 0, 3,
                                                                       62962, 28826, 62992, 7754,
                                                                       7772, 29600, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64312, 0, 3,
                                                                       63022, 28862, 63052, 7808,
                                                                       7826, 29636, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64372, 0, 3,
                                                                       63052, 28880, 63082, 7826,
                                                                       7844, 29672, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64432, 0, 3,
                                                                       63082, 28898, 63112, 7844,
                                                                       7862, 29708, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64492, 0, 3,
                                                                       63112, 28916, 63142, 7862,
                                                                       7880, 29744, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64552, 0, 3,
                                                                       63142, 28934, 63172, 7880,
                                                                       7898, 29780, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64612, 0, 3,
                                                                       63172, 28952, 63202, 7898,
                                                                       7916, 29816, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64672, 0, 3,
                                                                       63202, 28970, 63232, 7916,
                                                                       7934, 29852, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64732, 0, 3,
                                                                       63232, 28988, 63262, 7934,
                                                                       7952, 29888, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64792, 0, 3,
                                                                       63262, 29006, 63292, 7952,
                                                                       7970, 29924, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64852, 0, 3,
                                                                       63292, 29024, 63322, 7970,
                                                                       7988, 29960, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64912, 0, 3,
                                                                       63322, 29042, 63352, 7988,
                                                                       8006, 29996, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 64972, 0, 3,
                                                                       63352, 29060, 63382, 8006,
                                                                       8024, 30032, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 65032, 0, 3,
                                                                       63382, 29078, 63412, 8024,
                                                                       8042, 30068, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 65092, 0, 3,
                                                                       63412, 29096, 63442, 8042,
                                                                       8060, 30104, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65152, 0, 3,
                                                                       63472, 29132, 63532, 8096,
                                                                       8126, 30140, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65252, 0, 3,
                                                                       63532, 29168, 63592, 8126,
                                                                       8156, 30200, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65352, 0, 3,
                                                                       63592, 29204, 63652, 8156,
                                                                       8186, 30260, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65452, 0, 3,
                                                                       63652, 29240, 63712, 8186,
                                                                       8216, 30320, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65552, 0, 3,
                                                                       63712, 29276, 63772, 8216,
                                                                       8246, 30380, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65652, 0, 3,
                                                                       63772, 29312, 63832, 8246,
                                                                       8276, 30440, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65752, 0, 3,
                                                                       63832, 29348, 63892, 8276,
                                                                       8306, 30500, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65852, 0, 3,
                                                                       63892, 29384, 63952, 8306,
                                                                       8336, 30560, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 65952, 0, 3,
                                                                       63952, 29420, 64012, 8336,
                                                                       8366, 30620, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66052, 0, 3,
                                                                       64012, 29456, 64072, 8366,
                                                                       8396, 30680, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66152, 0, 3,
                                                                       64072, 29492, 64132, 8396,
                                                                       8426, 30740, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66252, 0, 3,
                                                                       64132, 29528, 64192, 8426,
                                                                       8456, 30800, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66352, 0, 3,
                                                                       64192, 29564, 64252, 8456,
                                                                       8486, 30860, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66452, 0, 3,
                                                                       64312, 29636, 64372, 8546,
                                                                       8576, 30920, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66552, 0, 3,
                                                                       64372, 29672, 64432, 8576,
                                                                       8606, 30980, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66652, 0, 3,
                                                                       64432, 29708, 64492, 8606,
                                                                       8636, 31040, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66752, 0, 3,
                                                                       64492, 29744, 64552, 8636,
                                                                       8666, 31100, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66852, 0, 3,
                                                                       64552, 29780, 64612, 8666,
                                                                       8696, 31160, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 66952, 0, 3,
                                                                       64612, 29816, 64672, 8696,
                                                                       8726, 31220, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 67052, 0, 3,
                                                                       64672, 29852, 64732, 8726,
                                                                       8756, 31280, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 67152, 0, 3,
                                                                       64732, 29888, 64792, 8756,
                                                                       8786, 31340, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 67252, 0, 3,
                                                                       64792, 29924, 64852, 8786,
                                                                       8816, 31400, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 67352, 0, 3,
                                                                       64852, 29960, 64912, 8816,
                                                                       8846, 31460, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 67452, 0, 3,
                                                                       64912, 29996, 64972, 8846,
                                                                       8876, 31520, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 67552, 0, 3,
                                                                       64972, 30032, 65032, 8876,
                                                                       8906, 31580, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 67652, 0, 3,
                                                                       65032, 30068, 65092, 8906,
                                                                       8936, 31640, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 67752, 0, 3,
                                                                       65152, 30140, 65252, 8996,
                                                                       9041, 31700, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 67902, 0, 3,
                                                                       65252, 30200, 65352, 9041,
                                                                       9086, 31790, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68052, 0, 3,
                                                                       65352, 30260, 65452, 9086,
                                                                       9131, 31880, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68202, 0, 3,
                                                                       65452, 30320, 65552, 9131,
                                                                       9176, 31970, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68352, 0, 3,
                                                                       65552, 30380, 65652, 9176,
                                                                       9221, 32060, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68502, 0, 3,
                                                                       65652, 30440, 65752, 9221,
                                                                       9266, 32150, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68652, 0, 3,
                                                                       65752, 30500, 65852, 9266,
                                                                       9311, 32240, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68802, 0, 3,
                                                                       65852, 30560, 65952, 9311,
                                                                       9356, 32330, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 68952, 0, 3,
                                                                       65952, 30620, 66052, 9356,
                                                                       9401, 32420, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69102, 0, 3,
                                                                       66052, 30680, 66152, 9401,
                                                                       9446, 32510, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69252, 0, 3,
                                                                       66152, 30740, 66252, 9446,
                                                                       9491, 32600, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69402, 0, 3,
                                                                       66252, 30800, 66352, 9491,
                                                                       9536, 32690, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69552, 0, 3,
                                                                       66452, 30920, 66552, 9626,
                                                                       9671, 32780, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69702, 0, 3,
                                                                       66552, 30980, 66652, 9671,
                                                                       9716, 32870, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 69852, 0, 3,
                                                                       66652, 31040, 66752, 9716,
                                                                       9761, 32960, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 70002, 0, 3,
                                                                       66752, 31100, 66852, 9761,
                                                                       9806, 33050, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 70152, 0, 3,
                                                                       66852, 31160, 66952, 9806,
                                                                       9851, 33140, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 70302, 0, 3,
                                                                       66952, 31220, 67052, 9851,
                                                                       9896, 33230, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 70452, 0, 3,
                                                                       67052, 31280, 67152, 9896,
                                                                       9941, 33320, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 70602, 0, 3,
                                                                       67152, 31340, 67252, 9941,
                                                                       9986, 33410, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 70752, 0, 3,
                                                                       67252, 31400, 67352, 9986,
                                                                       10031, 33500, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 70902, 0, 3,
                                                                       67352, 31460, 67452,
                                                                       10031, 10076, 33590,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 71052, 0, 3,
                                                                       67452, 31520, 67552,
                                                                       10076, 10121, 33680,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 71202, 0, 3,
                                                                       67552, 31580, 67652,
                                                                       10121, 10166, 33770,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71352, 0, 3,
                                                                       67752, 31700, 67902,
                                                                       10256, 10319, 33860,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71562, 0, 3,
                                                                       67902, 31790, 68052,
                                                                       10319, 10382, 33986,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71772, 0, 3,
                                                                       68052, 31880, 68202,
                                                                       10382, 10445, 34112,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 71982, 0, 3,
                                                                       68202, 31970, 68352,
                                                                       10445, 10508, 34238,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 72192, 0, 3,
                                                                       68352, 32060, 68502,
                                                                       10508, 10571, 34364,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 72402, 0, 3,
                                                                       68502, 32150, 68652,
                                                                       10571, 10634, 34490,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 72612, 0, 3,
                                                                       68652, 32240, 68802,
                                                                       10634, 10697, 34616,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 72822, 0, 3,
                                                                       68802, 32330, 68952,
                                                                       10697, 10760, 34742,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73032, 0, 3,
                                                                       68952, 32420, 69102,
                                                                       10760, 10823, 34868,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73242, 0, 3,
                                                                       69102, 32510, 69252,
                                                                       10823, 10886, 34994,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73452, 0, 3,
                                                                       69252, 32600, 69402,
                                                                       10886, 10949, 35120,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73662, 0, 3,
                                                                       69552, 32780, 69702,
                                                                       11075, 11138, 35246,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 73872, 0, 3,
                                                                       69702, 32870, 69852,
                                                                       11138, 11201, 35372,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 74082, 0, 3,
                                                                       69852, 32960, 70002,
                                                                       11201, 11264, 35498,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 74292, 0, 3,
                                                                       70002, 33050, 70152,
                                                                       11264, 11327, 35624,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 74502, 0, 3,
                                                                       70152, 33140, 70302,
                                                                       11327, 11390, 35750,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 74712, 0, 3,
                                                                       70302, 33230, 70452,
                                                                       11390, 11453, 35876,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 74922, 0, 3,
                                                                       70452, 33320, 70602,
                                                                       11453, 11516, 36002,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 75132, 0, 3,
                                                                       70602, 33410, 70752,
                                                                       11516, 11579, 36128,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 75342, 0, 3,
                                                                       70752, 33500, 70902,
                                                                       11579, 11642, 36254,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 75552, 0, 3,
                                                                       70902, 33590, 71052,
                                                                       11642, 11705, 36380,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 75762, 0, 3,
                                                                       71052, 33680, 71202,
                                                                       11705, 11768, 36506,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 75972, 0, 3,
                                                                       71352, 33860, 71562,
                                                                       11894, 11978, 36632,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 76252, 0, 3,
                                                                       71562, 33986, 71772,
                                                                       11978, 12062, 36800,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 76532, 0, 3,
                                                                       71772, 34112, 71982,
                                                                       12062, 12146, 36968,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 76812, 0, 3,
                                                                       71982, 34238, 72192,
                                                                       12146, 12230, 37136,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 77092, 0, 3,
                                                                       72192, 34364, 72402,
                                                                       12230, 12314, 37304,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 77372, 0, 3,
                                                                       72402, 34490, 72612,
                                                                       12314, 12398, 37472,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 77652, 0, 3,
                                                                       72612, 34616, 72822,
                                                                       12398, 12482, 37640,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 77932, 0, 3,
                                                                       72822, 34742, 73032,
                                                                       12482, 12566, 37808,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 78212, 0, 3,
                                                                       73032, 34868, 73242,
                                                                       12566, 12650, 37976,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 78492, 0, 3,
                                                                       73242, 34994, 73452,
                                                                       12650, 12734, 38144,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 78772, 0, 3,
                                                                       73662, 35246, 73872,
                                                                       12902, 12986, 38312,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 79052, 0, 3,
                                                                       73872, 35372, 74082,
                                                                       12986, 13070, 38480,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 79332, 0, 3,
                                                                       74082, 35498, 74292,
                                                                       13070, 13154, 38648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 79612, 0, 3,
                                                                       74292, 35624, 74502,
                                                                       13154, 13238, 38816,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 79892, 0, 3,
                                                                       74502, 35750, 74712,
                                                                       13238, 13322, 38984,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 80172, 0, 3,
                                                                       74712, 35876, 74922,
                                                                       13322, 13406, 39152,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 80452, 0, 3,
                                                                       74922, 36002, 75132,
                                                                       13406, 13490, 39320,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 80732, 0, 3,
                                                                       75132, 36128, 75342,
                                                                       13490, 13574, 39488,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 81012, 0, 3,
                                                                       75342, 36254, 75552,
                                                                       13574, 13658, 39656,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 81292, 0, 3,
                                                                       75552, 36380, 75762,
                                                                       13658, 13742, 39824,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 81572, 0, 3,
                                                                       75972, 36632, 76252,
                                                                       13910, 14018, 39992,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 81932, 0, 3,
                                                                       76252, 36800, 76532,
                                                                       14018, 14126, 40208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 82292, 0, 3,
                                                                       76532, 36968, 76812,
                                                                       14126, 14234, 40424,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 82652, 0, 3,
                                                                       76812, 37136, 77092,
                                                                       14234, 14342, 40640,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 83012, 0, 3,
                                                                       77092, 37304, 77372,
                                                                       14342, 14450, 40856,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 83372, 0, 3,
                                                                       77372, 37472, 77652,
                                                                       14450, 14558, 41072,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 83732, 0, 3,
                                                                       77652, 37640, 77932,
                                                                       14558, 14666, 41288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 84092, 0, 3,
                                                                       77932, 37808, 78212,
                                                                       14666, 14774, 41504,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 84452, 0, 3,
                                                                       78212, 37976, 78492,
                                                                       14774, 14882, 41720,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 84812, 0, 3,
                                                                       78772, 38312, 79052,
                                                                       15098, 15206, 41936,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 85172, 0, 3,
                                                                       79052, 38480, 79332,
                                                                       15206, 15314, 42152,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 85532, 0, 3,
                                                                       79332, 38648, 79612,
                                                                       15314, 15422, 42368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 85892, 0, 3,
                                                                       79612, 38816, 79892,
                                                                       15422, 15530, 42584,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 86252, 0, 3,
                                                                       79892, 38984, 80172,
                                                                       15530, 15638, 42800,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 86612, 0, 3,
                                                                       80172, 39152, 80452,
                                                                       15638, 15746, 43016,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 86972, 0, 3,
                                                                       80452, 39320, 80732,
                                                                       15746, 15854, 43232,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 87332, 0, 3,
                                                                       80732, 39488, 81012,
                                                                       15854, 15962, 43448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 87692, 0, 3,
                                                                       81012, 39656, 81292,
                                                                       15962, 16070, 43664,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 88052, 0, 3,
                                                                       81572, 39992, 81932,
                                                                       16286, 16421, 43880,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 88502, 0, 3,
                                                                       81932, 40208, 82292,
                                                                       16421, 16556, 44150,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 88952, 0, 3,
                                                                       82292, 40424, 82652,
                                                                       16556, 16691, 44420,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 89402, 0, 3,
                                                                       82652, 40640, 83012,
                                                                       16691, 16826, 44690,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 89852, 0, 3,
                                                                       83012, 40856, 83372,
                                                                       16826, 16961, 44960,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 90302, 0, 3,
                                                                       83372, 41072, 83732,
                                                                       16961, 17096, 45230,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 90752, 0, 3,
                                                                       83732, 41288, 84092,
                                                                       17096, 17231, 45500,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 91202, 0, 3,
                                                                       84092, 41504, 84452,
                                                                       17231, 17366, 45770,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 91652, 0, 3,
                                                                       84812, 41936, 85172,
                                                                       17636, 17771, 46040,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 92102, 0, 3,
                                                                       85172, 42152, 85532,
                                                                       17771, 17906, 46310,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 92552, 0, 3,
                                                                       85532, 42368, 85892,
                                                                       17906, 18041, 46580,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 93002, 0, 3,
                                                                       85892, 42584, 86252,
                                                                       18041, 18176, 46850,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 93452, 0, 3,
                                                                       86252, 42800, 86612,
                                                                       18176, 18311, 47120,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 93902, 0, 3,
                                                                       86612, 43016, 86972,
                                                                       18311, 18446, 47390,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 94352, 0, 3,
                                                                       86972, 43232, 87332,
                                                                       18446, 18581, 47660,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 94802, 0, 3,
                                                                       87332, 43448, 87692,
                                                                       18581, 18716, 47930,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 95252, 0, 3,
                                                                       88052, 43880, 88502,
                                                                       18986, 19151, 48200,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 95802, 0, 3,
                                                                       88502, 44150, 88952,
                                                                       19151, 19316, 48530,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 96352, 0, 3,
                                                                       88952, 44420, 89402,
                                                                       19316, 19481, 48860,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 96902, 0, 3,
                                                                       89402, 44690, 89852,
                                                                       19481, 19646, 49190,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 97452, 0, 3,
                                                                       89852, 44960, 90302,
                                                                       19646, 19811, 49520,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 98002, 0, 3,
                                                                       90302, 45230, 90752,
                                                                       19811, 19976, 49850,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 98552, 0, 3,
                                                                       90752, 45500, 91202,
                                                                       19976, 20141, 50180,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 99102, 0, 3,
                                                                       91652, 46040, 92102,
                                                                       20471, 20636, 50510,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 99652, 0, 3,
                                                                       92102, 46310, 92552,
                                                                       20636, 20801, 50840,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 100202, 0, 3,
                                                                       92552, 46580, 93002,
                                                                       20801, 20966, 51170,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 100752, 0, 3,
                                                                       93002, 46850, 93452,
                                                                       20966, 21131, 51500,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 101302, 0, 3,
                                                                       93452, 47120, 93902,
                                                                       21131, 21296, 51830,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 101852, 0, 3,
                                                                       93902, 47390, 94352,
                                                                       21296, 21461, 52160,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 102402, 0, 3,
                                                                       94352, 47660, 94802,
                                                                       21461, 21626, 52490,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 102952, 0, 3,
                                                                       95252, 48200, 95802,
                                                                       21956, 22154, 52820,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 103612, 0, 3,
                                                                       95802, 48530, 96352,
                                                                       22154, 22352, 53216,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 104272, 0, 3,
                                                                       96352, 48860, 96902,
                                                                       22352, 22550, 53612,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 104932, 0, 3,
                                                                       96902, 49190, 97452,
                                                                       22550, 22748, 54008,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 105592, 0, 3,
                                                                       97452, 49520, 98002,
                                                                       22748, 22946, 54404,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 106252, 0, 3,
                                                                       98002, 49850, 98552,
                                                                       22946, 23144, 54800,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 106912, 0, 3,
                                                                       99102, 50510, 99652,
                                                                       23540, 23738, 55196,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 107572, 0, 3,
                                                                       99652, 50840, 100202,
                                                                       23738, 23936, 55592,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 108232, 0, 3,
                                                                       100202, 51170, 100752,
                                                                       23936, 24134, 55988,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 108892, 0, 3,
                                                                       100752, 51500, 101302,
                                                                       24134, 24332, 56384,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 109552, 0, 3,
                                                                       101302, 51830, 101852,
                                                                       24332, 24530, 56780,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 110212, 0, 3,
                                                                       101852, 52160, 102402,
                                                                       24530, 24728, 57176,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 110872, 0, 3,
                                                                       102952, 52820, 103612,
                                                                       25124, 25358, 57572,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 111652, 0, 3,
                                                                       103612, 53216, 104272,
                                                                       25358, 25592, 58040,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 112432, 0, 3,
                                                                       104272, 53612, 104932,
                                                                       25592, 25826, 58508,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 113212, 0, 3,
                                                                       104932, 54008, 105592,
                                                                       25826, 26060, 58976,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 113992, 0, 3,
                                                                       105592, 54404, 106252,
                                                                       26060, 26294, 59444,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 114772, 0, 3,
                                                                       106912, 55196, 107572,
                                                                       26762, 26996, 59912,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 115552, 0, 3,
                                                                       107572, 55592, 108232,
                                                                       26996, 27230, 60380,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 116332, 0, 3,
                                                                       108232, 55988, 108892,
                                                                       27230, 27464, 60848,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 117112, 0, 3,
                                                                       108892, 56384, 109552,
                                                                       27464, 27698, 61316,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 117892, 0, 3,
                                                                       109552, 56780, 110212,
                                                                       27698, 27932, 61784,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118672, 3, 28400,
                                                                       28406, 62272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118687, 3, 28406,
                                                                       28412, 62282, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118702, 3, 28412,
                                                                       28418, 62292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118717, 3, 28418,
                                                                       28424, 62302, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118732, 3, 28424,
                                                                       28430, 62312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118747, 3, 28430,
                                                                       28436, 62322, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118762, 3, 28436,
                                                                       28442, 62332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118777, 3, 28442,
                                                                       28448, 62342, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118792, 3, 28448,
                                                                       28454, 62352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118807, 3, 28454,
                                                                       28460, 62362, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118822, 3, 28460,
                                                                       28466, 62372, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118837, 3, 28466,
                                                                       28472, 62382, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118852, 3, 28472,
                                                                       28478, 62392, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118867, 3, 28478,
                                                                       28484, 62402, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118882, 3, 28496,
                                                                       28502, 62432, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118897, 3, 28502,
                                                                       28508, 62442, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118912, 3, 28508,
                                                                       28514, 62452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118927, 3, 28514,
                                                                       28520, 62462, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118942, 3, 28520,
                                                                       28526, 62472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118957, 3, 28526,
                                                                       28532, 62482, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118972, 3, 28532,
                                                                       28538, 62492, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118987, 3, 28538,
                                                                       28544, 62502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119002, 3, 28544,
                                                                       28550, 62512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119017, 3, 28550,
                                                                       28556, 62522, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119032, 3, 28556,
                                                                       28562, 62532, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119047, 3, 28562,
                                                                       28568, 62542, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119062, 3, 28568,
                                                                       28574, 62552, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119077, 3, 28574,
                                                                       28580, 62562, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119092, 0, 3,
                                                                       118672, 62272, 118687,
                                                                       28592, 28610, 62632,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119137, 0, 3,
                                                                       118687, 62282, 118702,
                                                                       28610, 28628, 62662,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119182, 0, 3,
                                                                       118702, 62292, 118717,
                                                                       28628, 28646, 62692,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119227, 0, 3,
                                                                       118717, 62302, 118732,
                                                                       28646, 28664, 62722,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119272, 0, 3,
                                                                       118732, 62312, 118747,
                                                                       28664, 28682, 62752,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119317, 0, 3,
                                                                       118747, 62322, 118762,
                                                                       28682, 28700, 62782,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119362, 0, 3,
                                                                       118762, 62332, 118777,
                                                                       28700, 28718, 62812,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119407, 0, 3,
                                                                       118777, 62342, 118792,
                                                                       28718, 28736, 62842,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119452, 0, 3,
                                                                       118792, 62352, 118807,
                                                                       28736, 28754, 62872,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119497, 0, 3,
                                                                       118807, 62362, 118822,
                                                                       28754, 28772, 62902,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119542, 0, 3,
                                                                       118822, 62372, 118837,
                                                                       28772, 28790, 62932,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119587, 0, 3,
                                                                       118837, 62382, 118852,
                                                                       28790, 28808, 62962,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119632, 0, 3,
                                                                       118852, 62392, 118867,
                                                                       28808, 28826, 62992,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119677, 0, 3,
                                                                       118882, 62432, 118897,
                                                                       28862, 28880, 63082,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119722, 0, 3,
                                                                       118897, 62442, 118912,
                                                                       28880, 28898, 63112,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119767, 0, 3,
                                                                       118912, 62452, 118927,
                                                                       28898, 28916, 63142,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119812, 0, 3,
                                                                       118927, 62462, 118942,
                                                                       28916, 28934, 63172,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119857, 0, 3,
                                                                       118942, 62472, 118957,
                                                                       28934, 28952, 63202,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119902, 0, 3,
                                                                       118957, 62482, 118972,
                                                                       28952, 28970, 63232,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119947, 0, 3,
                                                                       118972, 62492, 118987,
                                                                       28970, 28988, 63262,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 119992, 0, 3,
                                                                       118987, 62502, 119002,
                                                                       28988, 29006, 63292,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 120037, 0, 3,
                                                                       119002, 62512, 119017,
                                                                       29006, 29024, 63322,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 120082, 0, 3,
                                                                       119017, 62522, 119032,
                                                                       29024, 29042, 63352,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 120127, 0, 3,
                                                                       119032, 62532, 119047,
                                                                       29042, 29060, 63382,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 120172, 0, 3,
                                                                       119047, 62542, 119062,
                                                                       29060, 29078, 63412,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 120217, 0, 3,
                                                                       119062, 62552, 119077,
                                                                       29078, 29096, 63442,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120262, 0, 3,
                                                                       119092, 62632, 119137,
                                                                       29132, 29168, 63592,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120352, 0, 3,
                                                                       119137, 62662, 119182,
                                                                       29168, 29204, 63652,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120442, 0, 3,
                                                                       119182, 62692, 119227,
                                                                       29204, 29240, 63712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120532, 0, 3,
                                                                       119227, 62722, 119272,
                                                                       29240, 29276, 63772,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120622, 0, 3,
                                                                       119272, 62752, 119317,
                                                                       29276, 29312, 63832,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120712, 0, 3,
                                                                       119317, 62782, 119362,
                                                                       29312, 29348, 63892,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120802, 0, 3,
                                                                       119362, 62812, 119407,
                                                                       29348, 29384, 63952,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120892, 0, 3,
                                                                       119407, 62842, 119452,
                                                                       29384, 29420, 64012,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 120982, 0, 3,
                                                                       119452, 62872, 119497,
                                                                       29420, 29456, 64072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121072, 0, 3,
                                                                       119497, 62902, 119542,
                                                                       29456, 29492, 64132,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121162, 0, 3,
                                                                       119542, 62932, 119587,
                                                                       29492, 29528, 64192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121252, 0, 3,
                                                                       119587, 62962, 119632,
                                                                       29528, 29564, 64252,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121342, 0, 3,
                                                                       119677, 63082, 119722,
                                                                       29636, 29672, 64432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121432, 0, 3,
                                                                       119722, 63112, 119767,
                                                                       29672, 29708, 64492,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121522, 0, 3,
                                                                       119767, 63142, 119812,
                                                                       29708, 29744, 64552,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121612, 0, 3,
                                                                       119812, 63172, 119857,
                                                                       29744, 29780, 64612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121702, 0, 3,
                                                                       119857, 63202, 119902,
                                                                       29780, 29816, 64672,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121792, 0, 3,
                                                                       119902, 63232, 119947,
                                                                       29816, 29852, 64732,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121882, 0, 3,
                                                                       119947, 63262, 119992,
                                                                       29852, 29888, 64792,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 121972, 0, 3,
                                                                       119992, 63292, 120037,
                                                                       29888, 29924, 64852,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 122062, 0, 3,
                                                                       120037, 63322, 120082,
                                                                       29924, 29960, 64912,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 122152, 0, 3,
                                                                       120082, 63352, 120127,
                                                                       29960, 29996, 64972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 122242, 0, 3,
                                                                       120127, 63382, 120172,
                                                                       29996, 30032, 65032,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 122332, 0, 3,
                                                                       120172, 63412, 120217,
                                                                       30032, 30068, 65092,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122422, 0, 3,
                                                                       120262, 63592, 120352,
                                                                       30140, 30200, 65352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122572, 0, 3,
                                                                       120352, 63652, 120442,
                                                                       30200, 30260, 65452,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122722, 0, 3,
                                                                       120442, 63712, 120532,
                                                                       30260, 30320, 65552,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 122872, 0, 3,
                                                                       120532, 63772, 120622,
                                                                       30320, 30380, 65652,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 123022, 0, 3,
                                                                       120622, 63832, 120712,
                                                                       30380, 30440, 65752,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 123172, 0, 3,
                                                                       120712, 63892, 120802,
                                                                       30440, 30500, 65852,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 123322, 0, 3,
                                                                       120802, 63952, 120892,
                                                                       30500, 30560, 65952,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 123472, 0, 3,
                                                                       120892, 64012, 120982,
                                                                       30560, 30620, 66052,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 123622, 0, 3,
                                                                       120982, 64072, 121072,
                                                                       30620, 30680, 66152,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 123772, 0, 3,
                                                                       121072, 64132, 121162,
                                                                       30680, 30740, 66252,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 123922, 0, 3,
                                                                       121162, 64192, 121252,
                                                                       30740, 30800, 66352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 124072, 0, 3,
                                                                       121342, 64432, 121432,
                                                                       30920, 30980, 66652,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 124222, 0, 3,
                                                                       121432, 64492, 121522,
                                                                       30980, 31040, 66752,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 124372, 0, 3,
                                                                       121522, 64552, 121612,
                                                                       31040, 31100, 66852,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 124522, 0, 3,
                                                                       121612, 64612, 121702,
                                                                       31100, 31160, 66952,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 124672, 0, 3,
                                                                       121702, 64672, 121792,
                                                                       31160, 31220, 67052,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 124822, 0, 3,
                                                                       121792, 64732, 121882,
                                                                       31220, 31280, 67152,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 124972, 0, 3,
                                                                       121882, 64792, 121972,
                                                                       31280, 31340, 67252,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 125122, 0, 3,
                                                                       121972, 64852, 122062,
                                                                       31340, 31400, 67352,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 125272, 0, 3,
                                                                       122062, 64912, 122152,
                                                                       31400, 31460, 67452,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 125422, 0, 3,
                                                                       122152, 64972, 122242,
                                                                       31460, 31520, 67552,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 125572, 0, 3,
                                                                       122242, 65032, 122332,
                                                                       31520, 31580, 67652,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 125722, 0, 3,
                                                                       122422, 65352, 122572,
                                                                       31700, 31790, 68052,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 125947, 0, 3,
                                                                       122572, 65452, 122722,
                                                                       31790, 31880, 68202,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 126172, 0, 3,
                                                                       122722, 65552, 122872,
                                                                       31880, 31970, 68352,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 126397, 0, 3,
                                                                       122872, 65652, 123022,
                                                                       31970, 32060, 68502,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 126622, 0, 3,
                                                                       123022, 65752, 123172,
                                                                       32060, 32150, 68652,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 126847, 0, 3,
                                                                       123172, 65852, 123322,
                                                                       32150, 32240, 68802,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 127072, 0, 3,
                                                                       123322, 65952, 123472,
                                                                       32240, 32330, 68952,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 127297, 0, 3,
                                                                       123472, 66052, 123622,
                                                                       32330, 32420, 69102,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 127522, 0, 3,
                                                                       123622, 66152, 123772,
                                                                       32420, 32510, 69252,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 127747, 0, 3,
                                                                       123772, 66252, 123922,
                                                                       32510, 32600, 69402,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 127972, 0, 3,
                                                                       124072, 66652, 124222,
                                                                       32780, 32870, 69852,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 128197, 0, 3,
                                                                       124222, 66752, 124372,
                                                                       32870, 32960, 70002,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 128422, 0, 3,
                                                                       124372, 66852, 124522,
                                                                       32960, 33050, 70152,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 128647, 0, 3,
                                                                       124522, 66952, 124672,
                                                                       33050, 33140, 70302,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 128872, 0, 3,
                                                                       124672, 67052, 124822,
                                                                       33140, 33230, 70452,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 129097, 0, 3,
                                                                       124822, 67152, 124972,
                                                                       33230, 33320, 70602,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 129322, 0, 3,
                                                                       124972, 67252, 125122,
                                                                       33320, 33410, 70752,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 129547, 0, 3,
                                                                       125122, 67352, 125272,
                                                                       33410, 33500, 70902,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 129772, 0, 3,
                                                                       125272, 67452, 125422,
                                                                       33500, 33590, 71052,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 129997, 0, 3,
                                                                       125422, 67552, 125572,
                                                                       33590, 33680, 71202,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 130222, 0, 3,
                                                                       125722, 68052, 125947,
                                                                       33860, 33986, 71772,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 130537, 0, 3,
                                                                       125947, 68202, 126172,
                                                                       33986, 34112, 71982,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 130852, 0, 3,
                                                                       126172, 68352, 126397,
                                                                       34112, 34238, 72192,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 131167, 0, 3,
                                                                       126397, 68502, 126622,
                                                                       34238, 34364, 72402,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 131482, 0, 3,
                                                                       126622, 68652, 126847,
                                                                       34364, 34490, 72612,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 131797, 0, 3,
                                                                       126847, 68802, 127072,
                                                                       34490, 34616, 72822,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 132112, 0, 3,
                                                                       127072, 68952, 127297,
                                                                       34616, 34742, 73032,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 132427, 0, 3,
                                                                       127297, 69102, 127522,
                                                                       34742, 34868, 73242,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 132742, 0, 3,
                                                                       127522, 69252, 127747,
                                                                       34868, 34994, 73452,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 133057, 0, 3,
                                                                       127972, 69852, 128197,
                                                                       35246, 35372, 74082,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 133372, 0, 3,
                                                                       128197, 70002, 128422,
                                                                       35372, 35498, 74292,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 133687, 0, 3,
                                                                       128422, 70152, 128647,
                                                                       35498, 35624, 74502,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 134002, 0, 3,
                                                                       128647, 70302, 128872,
                                                                       35624, 35750, 74712,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 134317, 0, 3,
                                                                       128872, 70452, 129097,
                                                                       35750, 35876, 74922,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 134632, 0, 3,
                                                                       129097, 70602, 129322,
                                                                       35876, 36002, 75132,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 134947, 0, 3,
                                                                       129322, 70752, 129547,
                                                                       36002, 36128, 75342,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 135262, 0, 3,
                                                                       129547, 70902, 129772,
                                                                       36128, 36254, 75552,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 135577, 0, 3,
                                                                       129772, 71052, 129997,
                                                                       36254, 36380, 75762,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 135892, 0, 3,
                                                                       130222, 71772, 130537,
                                                                       36632, 36800, 76532,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 136312, 0, 3,
                                                                       130537, 71982, 130852,
                                                                       36800, 36968, 76812,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 136732, 0, 3,
                                                                       130852, 72192, 131167,
                                                                       36968, 37136, 77092,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 137152, 0, 3,
                                                                       131167, 72402, 131482,
                                                                       37136, 37304, 77372,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 137572, 0, 3,
                                                                       131482, 72612, 131797,
                                                                       37304, 37472, 77652,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 137992, 0, 3,
                                                                       131797, 72822, 132112,
                                                                       37472, 37640, 77932,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 138412, 0, 3,
                                                                       132112, 73032, 132427,
                                                                       37640, 37808, 78212,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 138832, 0, 3,
                                                                       132427, 73242, 132742,
                                                                       37808, 37976, 78492,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 139252, 0, 3,
                                                                       133057, 74082, 133372,
                                                                       38312, 38480, 79332,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 139672, 0, 3,
                                                                       133372, 74292, 133687,
                                                                       38480, 38648, 79612,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 140092, 0, 3,
                                                                       133687, 74502, 134002,
                                                                       38648, 38816, 79892,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 140512, 0, 3,
                                                                       134002, 74712, 134317,
                                                                       38816, 38984, 80172,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 140932, 0, 3,
                                                                       134317, 74922, 134632,
                                                                       38984, 39152, 80452,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 141352, 0, 3,
                                                                       134632, 75132, 134947,
                                                                       39152, 39320, 80732,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 141772, 0, 3,
                                                                       134947, 75342, 135262,
                                                                       39320, 39488, 81012,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 142192, 0, 3,
                                                                       135262, 75552, 135577,
                                                                       39488, 39656, 81292,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 142612, 0, 3,
                                                                       135892, 76532, 136312,
                                                                       39992, 40208, 82292,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 143152, 0, 3,
                                                                       136312, 76812, 136732,
                                                                       40208, 40424, 82652,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 143692, 0, 3,
                                                                       136732, 77092, 137152,
                                                                       40424, 40640, 83012,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 144232, 0, 3,
                                                                       137152, 77372, 137572,
                                                                       40640, 40856, 83372,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 144772, 0, 3,
                                                                       137572, 77652, 137992,
                                                                       40856, 41072, 83732,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 145312, 0, 3,
                                                                       137992, 77932, 138412,
                                                                       41072, 41288, 84092,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 145852, 0, 3,
                                                                       138412, 78212, 138832,
                                                                       41288, 41504, 84452,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 146392, 0, 3,
                                                                       139252, 79332, 139672,
                                                                       41936, 42152, 85532,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 146932, 0, 3,
                                                                       139672, 79612, 140092,
                                                                       42152, 42368, 85892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 147472, 0, 3,
                                                                       140092, 79892, 140512,
                                                                       42368, 42584, 86252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 148012, 0, 3,
                                                                       140512, 80172, 140932,
                                                                       42584, 42800, 86612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 148552, 0, 3,
                                                                       140932, 80452, 141352,
                                                                       42800, 43016, 86972,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 149092, 0, 3,
                                                                       141352, 80732, 141772,
                                                                       43016, 43232, 87332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 149632, 0, 3,
                                                                       141772, 81012, 142192,
                                                                       43232, 43448, 87692,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 150172, 0, 3,
                                                                       142612, 82292, 143152,
                                                                       43880, 44150, 88952,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 150847, 0, 3,
                                                                       143152, 82652, 143692,
                                                                       44150, 44420, 89402,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 151522, 0, 3,
                                                                       143692, 83012, 144232,
                                                                       44420, 44690, 89852,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 152197, 0, 3,
                                                                       144232, 83372, 144772,
                                                                       44690, 44960, 90302,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 152872, 0, 3,
                                                                       144772, 83732, 145312,
                                                                       44960, 45230, 90752,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 153547, 0, 3,
                                                                       145312, 84092, 145852,
                                                                       45230, 45500, 91202,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 154222, 0, 3,
                                                                       146392, 85532, 146932,
                                                                       46040, 46310, 92552,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 154897, 0, 3,
                                                                       146932, 85892, 147472,
                                                                       46310, 46580, 93002,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 155572, 0, 3,
                                                                       147472, 86252, 148012,
                                                                       46580, 46850, 93452,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 156247, 0, 3,
                                                                       148012, 86612, 148552,
                                                                       46850, 47120, 93902,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 156922, 0, 3,
                                                                       148552, 86972, 149092,
                                                                       47120, 47390, 94352,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 157597, 0, 3,
                                                                       149092, 87332, 149632,
                                                                       47390, 47660, 94802,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 158272, 0, 3,
                                                                       150172, 88952, 150847,
                                                                       48200, 48530, 96352,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 159097, 0, 3,
                                                                       150847, 89402, 151522,
                                                                       48530, 48860, 96902,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 159922, 0, 3,
                                                                       151522, 89852, 152197,
                                                                       48860, 49190, 97452,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 160747, 0, 3,
                                                                       152197, 90302, 152872,
                                                                       49190, 49520, 98002,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 161572, 0, 3,
                                                                       152872, 90752, 153547,
                                                                       49520, 49850, 98552,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 162397, 0, 3,
                                                                       154222, 92552, 154897,
                                                                       50510, 50840, 100202,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 163222, 0, 3,
                                                                       154897, 93002, 155572,
                                                                       50840, 51170, 100752,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 164047, 0, 3,
                                                                       155572, 93452, 156247,
                                                                       51170, 51500, 101302,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 164872, 0, 3,
                                                                       156247, 93902, 156922,
                                                                       51500, 51830, 101852,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 165697, 0, 3,
                                                                       156922, 94352, 157597,
                                                                       51830, 52160, 102402,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 166522, 0, 3,
                                                                       158272, 96352, 159097,
                                                                       52820, 53216, 104272,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 167512, 0, 3,
                                                                       159097, 96902, 159922,
                                                                       53216, 53612, 104932,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 168502, 0, 3,
                                                                       159922, 97452, 160747,
                                                                       53612, 54008, 105592,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 169492, 0, 3,
                                                                       160747, 98002, 161572,
                                                                       54008, 54404, 106252,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 170482, 0, 3,
                                                                       162397, 100202, 163222,
                                                                       55196, 55592, 108232,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 171472, 0, 3,
                                                                       163222, 100752, 164047,
                                                                       55592, 55988, 108892,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 172462, 0, 3,
                                                                       164047, 101302, 164872,
                                                                       55988, 56384, 109552,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 173452, 0, 3,
                                                                       164872, 101852, 165697,
                                                                       56384, 56780, 110212,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 174442, 0, 3,
                                                                       166522, 104272, 167512,
                                                                       57572, 58040, 112432,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 175612, 0, 3,
                                                                       167512, 104932, 168502,
                                                                       58040, 58508, 113212,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 176782, 0, 3,
                                                                       168502, 105592, 169492,
                                                                       58508, 58976, 113992,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 177952, 0, 3,
                                                                       170482, 108232, 171472,
                                                                       59912, 60380, 116332,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 179122, 0, 3,
                                                                       171472, 108892, 172462,
                                                                       60380, 60848, 117112,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 180292, 0, 3,
                                                                       172462, 109552, 173452,
                                                                       60848, 61316, 117892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181462, 3, 62252,
                                                                       62262, 118672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181483, 3, 62262,
                                                                       62272, 118687, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181504, 3, 62272,
                                                                       62282, 118702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181525, 3, 62282,
                                                                       62292, 118717, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181546, 3, 62292,
                                                                       62302, 118732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181567, 3, 62302,
                                                                       62312, 118747, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181588, 3, 62312,
                                                                       62322, 118762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181609, 3, 62322,
                                                                       62332, 118777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181630, 3, 62332,
                                                                       62342, 118792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181651, 3, 62342,
                                                                       62352, 118807, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181672, 3, 62352,
                                                                       62362, 118822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181693, 3, 62362,
                                                                       62372, 118837, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181714, 3, 62372,
                                                                       62382, 118852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181735, 3, 62382,
                                                                       62392, 118867, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181756, 3, 62412,
                                                                       62422, 118882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181777, 3, 62422,
                                                                       62432, 118897, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181798, 3, 62432,
                                                                       62442, 118912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181819, 3, 62442,
                                                                       62452, 118927, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181840, 3, 62452,
                                                                       62462, 118942, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181861, 3, 62462,
                                                                       62472, 118957, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181882, 3, 62472,
                                                                       62482, 118972, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181903, 3, 62482,
                                                                       62492, 118987, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181924, 3, 62492,
                                                                       62502, 119002, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181945, 3, 62502,
                                                                       62512, 119017, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181966, 3, 62512,
                                                                       62522, 119032, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181987, 3, 62522,
                                                                       62532, 119047, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 182008, 3, 62532,
                                                                       62542, 119062, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 182029, 3, 62542,
                                                                       62552, 119077, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182050, 0, 3,
                                                                       181462, 118672, 181483,
                                                                       62572, 62602, 119092,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182113, 0, 3,
                                                                       181483, 118687, 181504,
                                                                       62602, 62632, 119137,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182176, 0, 3,
                                                                       181504, 118702, 181525,
                                                                       62632, 62662, 119182,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182239, 0, 3,
                                                                       181525, 118717, 181546,
                                                                       62662, 62692, 119227,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182302, 0, 3,
                                                                       181546, 118732, 181567,
                                                                       62692, 62722, 119272,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182365, 0, 3,
                                                                       181567, 118747, 181588,
                                                                       62722, 62752, 119317,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182428, 0, 3,
                                                                       181588, 118762, 181609,
                                                                       62752, 62782, 119362,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182491, 0, 3,
                                                                       181609, 118777, 181630,
                                                                       62782, 62812, 119407,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182554, 0, 3,
                                                                       181630, 118792, 181651,
                                                                       62812, 62842, 119452,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182617, 0, 3,
                                                                       181651, 118807, 181672,
                                                                       62842, 62872, 119497,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182680, 0, 3,
                                                                       181672, 118822, 181693,
                                                                       62872, 62902, 119542,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182743, 0, 3,
                                                                       181693, 118837, 181714,
                                                                       62902, 62932, 119587,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182806, 0, 3,
                                                                       181714, 118852, 181735,
                                                                       62932, 62962, 119632,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182869, 0, 3,
                                                                       181756, 118882, 181777,
                                                                       63022, 63052, 119677,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182932, 0, 3,
                                                                       181777, 118897, 181798,
                                                                       63052, 63082, 119722,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 182995, 0, 3,
                                                                       181798, 118912, 181819,
                                                                       63082, 63112, 119767,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183058, 0, 3,
                                                                       181819, 118927, 181840,
                                                                       63112, 63142, 119812,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183121, 0, 3,
                                                                       181840, 118942, 181861,
                                                                       63142, 63172, 119857,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183184, 0, 3,
                                                                       181861, 118957, 181882,
                                                                       63172, 63202, 119902,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183247, 0, 3,
                                                                       181882, 118972, 181903,
                                                                       63202, 63232, 119947,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183310, 0, 3,
                                                                       181903, 118987, 181924,
                                                                       63232, 63262, 119992,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183373, 0, 3,
                                                                       181924, 119002, 181945,
                                                                       63262, 63292, 120037,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183436, 0, 3,
                                                                       181945, 119017, 181966,
                                                                       63292, 63322, 120082,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183499, 0, 3,
                                                                       181966, 119032, 181987,
                                                                       63322, 63352, 120127,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183562, 0, 3,
                                                                       181987, 119047, 182008,
                                                                       63352, 63382, 120172,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 183625, 0, 3,
                                                                       182008, 119062, 182029,
                                                                       63382, 63412, 120217,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 183688, 0, 3,
                                                                       182050, 119092, 182113,
                                                                       63472, 63532, 120262,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 183814, 0, 3,
                                                                       182113, 119137, 182176,
                                                                       63532, 63592, 120352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 183940, 0, 3,
                                                                       182176, 119182, 182239,
                                                                       63592, 63652, 120442,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 184066, 0, 3,
                                                                       182239, 119227, 182302,
                                                                       63652, 63712, 120532,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 184192, 0, 3,
                                                                       182302, 119272, 182365,
                                                                       63712, 63772, 120622,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 184318, 0, 3,
                                                                       182365, 119317, 182428,
                                                                       63772, 63832, 120712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 184444, 0, 3,
                                                                       182428, 119362, 182491,
                                                                       63832, 63892, 120802,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 184570, 0, 3,
                                                                       182491, 119407, 182554,
                                                                       63892, 63952, 120892,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 184696, 0, 3,
                                                                       182554, 119452, 182617,
                                                                       63952, 64012, 120982,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 184822, 0, 3,
                                                                       182617, 119497, 182680,
                                                                       64012, 64072, 121072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 184948, 0, 3,
                                                                       182680, 119542, 182743,
                                                                       64072, 64132, 121162,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 185074, 0, 3,
                                                                       182743, 119587, 182806,
                                                                       64132, 64192, 121252,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 185200, 0, 3,
                                                                       182869, 119677, 182932,
                                                                       64312, 64372, 121342,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 185326, 0, 3,
                                                                       182932, 119722, 182995,
                                                                       64372, 64432, 121432,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 185452, 0, 3,
                                                                       182995, 119767, 183058,
                                                                       64432, 64492, 121522,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 185578, 0, 3,
                                                                       183058, 119812, 183121,
                                                                       64492, 64552, 121612,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 185704, 0, 3,
                                                                       183121, 119857, 183184,
                                                                       64552, 64612, 121702,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 185830, 0, 3,
                                                                       183184, 119902, 183247,
                                                                       64612, 64672, 121792,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 185956, 0, 3,
                                                                       183247, 119947, 183310,
                                                                       64672, 64732, 121882,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 186082, 0, 3,
                                                                       183310, 119992, 183373,
                                                                       64732, 64792, 121972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 186208, 0, 3,
                                                                       183373, 120037, 183436,
                                                                       64792, 64852, 122062,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 186334, 0, 3,
                                                                       183436, 120082, 183499,
                                                                       64852, 64912, 122152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 186460, 0, 3,
                                                                       183499, 120127, 183562,
                                                                       64912, 64972, 122242,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 186586, 0, 3,
                                                                       183562, 120172, 183625,
                                                                       64972, 65032, 122332,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 186712, 0, 3,
                                                                       183688, 120262, 183814,
                                                                       65152, 65252, 122422,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 186922, 0, 3,
                                                                       183814, 120352, 183940,
                                                                       65252, 65352, 122572,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 187132, 0, 3,
                                                                       183940, 120442, 184066,
                                                                       65352, 65452, 122722,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 187342, 0, 3,
                                                                       184066, 120532, 184192,
                                                                       65452, 65552, 122872,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 187552, 0, 3,
                                                                       184192, 120622, 184318,
                                                                       65552, 65652, 123022,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 187762, 0, 3,
                                                                       184318, 120712, 184444,
                                                                       65652, 65752, 123172,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 187972, 0, 3,
                                                                       184444, 120802, 184570,
                                                                       65752, 65852, 123322,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 188182, 0, 3,
                                                                       184570, 120892, 184696,
                                                                       65852, 65952, 123472,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 188392, 0, 3,
                                                                       184696, 120982, 184822,
                                                                       65952, 66052, 123622,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 188602, 0, 3,
                                                                       184822, 121072, 184948,
                                                                       66052, 66152, 123772,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 188812, 0, 3,
                                                                       184948, 121162, 185074,
                                                                       66152, 66252, 123922,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 189022, 0, 3,
                                                                       185200, 121342, 185326,
                                                                       66452, 66552, 124072,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 189232, 0, 3,
                                                                       185326, 121432, 185452,
                                                                       66552, 66652, 124222,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 189442, 0, 3,
                                                                       185452, 121522, 185578,
                                                                       66652, 66752, 124372,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 189652, 0, 3,
                                                                       185578, 121612, 185704,
                                                                       66752, 66852, 124522,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 189862, 0, 3,
                                                                       185704, 121702, 185830,
                                                                       66852, 66952, 124672,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 190072, 0, 3,
                                                                       185830, 121792, 185956,
                                                                       66952, 67052, 124822,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 190282, 0, 3,
                                                                       185956, 121882, 186082,
                                                                       67052, 67152, 124972,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 190492, 0, 3,
                                                                       186082, 121972, 186208,
                                                                       67152, 67252, 125122,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 190702, 0, 3,
                                                                       186208, 122062, 186334,
                                                                       67252, 67352, 125272,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 190912, 0, 3,
                                                                       186334, 122152, 186460,
                                                                       67352, 67452, 125422,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 191122, 0, 3,
                                                                       186460, 122242, 186586,
                                                                       67452, 67552, 125572,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 191332, 0, 3,
                                                                       186712, 122422, 186922,
                                                                       67752, 67902, 125722,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 191647, 0, 3,
                                                                       186922, 122572, 187132,
                                                                       67902, 68052, 125947,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 191962, 0, 3,
                                                                       187132, 122722, 187342,
                                                                       68052, 68202, 126172,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 192277, 0, 3,
                                                                       187342, 122872, 187552,
                                                                       68202, 68352, 126397,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 192592, 0, 3,
                                                                       187552, 123022, 187762,
                                                                       68352, 68502, 126622,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 192907, 0, 3,
                                                                       187762, 123172, 187972,
                                                                       68502, 68652, 126847,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 193222, 0, 3,
                                                                       187972, 123322, 188182,
                                                                       68652, 68802, 127072,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 193537, 0, 3,
                                                                       188182, 123472, 188392,
                                                                       68802, 68952, 127297,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 193852, 0, 3,
                                                                       188392, 123622, 188602,
                                                                       68952, 69102, 127522,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 194167, 0, 3,
                                                                       188602, 123772, 188812,
                                                                       69102, 69252, 127747,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 194482, 0, 3,
                                                                       189022, 124072, 189232,
                                                                       69552, 69702, 127972,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 194797, 0, 3,
                                                                       189232, 124222, 189442,
                                                                       69702, 69852, 128197,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 195112, 0, 3,
                                                                       189442, 124372, 189652,
                                                                       69852, 70002, 128422,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 195427, 0, 3,
                                                                       189652, 124522, 189862,
                                                                       70002, 70152, 128647,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 195742, 0, 3,
                                                                       189862, 124672, 190072,
                                                                       70152, 70302, 128872,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 196057, 0, 3,
                                                                       190072, 124822, 190282,
                                                                       70302, 70452, 129097,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 196372, 0, 3,
                                                                       190282, 124972, 190492,
                                                                       70452, 70602, 129322,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 196687, 0, 3,
                                                                       190492, 125122, 190702,
                                                                       70602, 70752, 129547,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 197002, 0, 3,
                                                                       190702, 125272, 190912,
                                                                       70752, 70902, 129772,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 197317, 0, 3,
                                                                       190912, 125422, 191122,
                                                                       70902, 71052, 129997,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 197632, 0, 3,
                                                                       191332, 125722, 191647,
                                                                       71352, 71562, 130222,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 198073, 0, 3,
                                                                       191647, 125947, 191962,
                                                                       71562, 71772, 130537,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 198514, 0, 3,
                                                                       191962, 126172, 192277,
                                                                       71772, 71982, 130852,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 198955, 0, 3,
                                                                       192277, 126397, 192592,
                                                                       71982, 72192, 131167,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 199396, 0, 3,
                                                                       192592, 126622, 192907,
                                                                       72192, 72402, 131482,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 199837, 0, 3,
                                                                       192907, 126847, 193222,
                                                                       72402, 72612, 131797,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 200278, 0, 3,
                                                                       193222, 127072, 193537,
                                                                       72612, 72822, 132112,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 200719, 0, 3,
                                                                       193537, 127297, 193852,
                                                                       72822, 73032, 132427,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 201160, 0, 3,
                                                                       193852, 127522, 194167,
                                                                       73032, 73242, 132742,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 201601, 0, 3,
                                                                       194482, 127972, 194797,
                                                                       73662, 73872, 133057,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 202042, 0, 3,
                                                                       194797, 128197, 195112,
                                                                       73872, 74082, 133372,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 202483, 0, 3,
                                                                       195112, 128422, 195427,
                                                                       74082, 74292, 133687,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 202924, 0, 3,
                                                                       195427, 128647, 195742,
                                                                       74292, 74502, 134002,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 203365, 0, 3,
                                                                       195742, 128872, 196057,
                                                                       74502, 74712, 134317,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 203806, 0, 3,
                                                                       196057, 129097, 196372,
                                                                       74712, 74922, 134632,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 204247, 0, 3,
                                                                       196372, 129322, 196687,
                                                                       74922, 75132, 134947,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 204688, 0, 3,
                                                                       196687, 129547, 197002,
                                                                       75132, 75342, 135262,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 205129, 0, 3,
                                                                       197002, 129772, 197317,
                                                                       75342, 75552, 135577,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 205570, 0, 3,
                                                                       197632, 130222, 198073,
                                                                       75972, 76252, 135892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 206158, 0, 3,
                                                                       198073, 130537, 198514,
                                                                       76252, 76532, 136312,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 206746, 0, 3,
                                                                       198514, 130852, 198955,
                                                                       76532, 76812, 136732,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 207334, 0, 3,
                                                                       198955, 131167, 199396,
                                                                       76812, 77092, 137152,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 207922, 0, 3,
                                                                       199396, 131482, 199837,
                                                                       77092, 77372, 137572,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 208510, 0, 3,
                                                                       199837, 131797, 200278,
                                                                       77372, 77652, 137992,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 209098, 0, 3,
                                                                       200278, 132112, 200719,
                                                                       77652, 77932, 138412,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 209686, 0, 3,
                                                                       200719, 132427, 201160,
                                                                       77932, 78212, 138832,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 210274, 0, 3,
                                                                       201601, 133057, 202042,
                                                                       78772, 79052, 139252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 210862, 0, 3,
                                                                       202042, 133372, 202483,
                                                                       79052, 79332, 139672,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 211450, 0, 3,
                                                                       202483, 133687, 202924,
                                                                       79332, 79612, 140092,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 212038, 0, 3,
                                                                       202924, 134002, 203365,
                                                                       79612, 79892, 140512,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 212626, 0, 3,
                                                                       203365, 134317, 203806,
                                                                       79892, 80172, 140932,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 213214, 0, 3,
                                                                       203806, 134632, 204247,
                                                                       80172, 80452, 141352,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 213802, 0, 3,
                                                                       204247, 134947, 204688,
                                                                       80452, 80732, 141772,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 214390, 0, 3,
                                                                       204688, 135262, 205129,
                                                                       80732, 81012, 142192,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 214978, 0, 3,
                                                                       205570, 135892, 206158,
                                                                       81572, 81932, 142612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 215734, 0, 3,
                                                                       206158, 136312, 206746,
                                                                       81932, 82292, 143152,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 216490, 0, 3,
                                                                       206746, 136732, 207334,
                                                                       82292, 82652, 143692,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 217246, 0, 3,
                                                                       207334, 137152, 207922,
                                                                       82652, 83012, 144232,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 218002, 0, 3,
                                                                       207922, 137572, 208510,
                                                                       83012, 83372, 144772,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 218758, 0, 3,
                                                                       208510, 137992, 209098,
                                                                       83372, 83732, 145312,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 219514, 0, 3,
                                                                       209098, 138412, 209686,
                                                                       83732, 84092, 145852,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 220270, 0, 3,
                                                                       210274, 139252, 210862,
                                                                       84812, 85172, 146392,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 221026, 0, 3,
                                                                       210862, 139672, 211450,
                                                                       85172, 85532, 146932,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 221782, 0, 3,
                                                                       211450, 140092, 212038,
                                                                       85532, 85892, 147472,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 222538, 0, 3,
                                                                       212038, 140512, 212626,
                                                                       85892, 86252, 148012,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 223294, 0, 3,
                                                                       212626, 140932, 213214,
                                                                       86252, 86612, 148552,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 224050, 0, 3,
                                                                       213214, 141352, 213802,
                                                                       86612, 86972, 149092,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 224806, 0, 3,
                                                                       213802, 141772, 214390,
                                                                       86972, 87332, 149632,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 225562, 0, 3,
                                                                       214978, 142612, 215734,
                                                                       88052, 88502, 150172,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 226507, 0, 3,
                                                                       215734, 143152, 216490,
                                                                       88502, 88952, 150847,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 227452, 0, 3,
                                                                       216490, 143692, 217246,
                                                                       88952, 89402, 151522,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 228397, 0, 3,
                                                                       217246, 144232, 218002,
                                                                       89402, 89852, 152197,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 229342, 0, 3,
                                                                       218002, 144772, 218758,
                                                                       89852, 90302, 152872,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 230287, 0, 3,
                                                                       218758, 145312, 219514,
                                                                       90302, 90752, 153547,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 231232, 0, 3,
                                                                       220270, 146392, 221026,
                                                                       91652, 92102, 154222,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 232177, 0, 3,
                                                                       221026, 146932, 221782,
                                                                       92102, 92552, 154897,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 233122, 0, 3,
                                                                       221782, 147472, 222538,
                                                                       92552, 93002, 155572,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 234067, 0, 3,
                                                                       222538, 148012, 223294,
                                                                       93002, 93452, 156247,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 235012, 0, 3,
                                                                       223294, 148552, 224050,
                                                                       93452, 93902, 156922,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 235957, 0, 3,
                                                                       224050, 149092, 224806,
                                                                       93902, 94352, 157597,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 236902, 0, 3,
                                                                       225562, 150172, 226507,
                                                                       95252, 95802, 158272,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 238057, 0, 3,
                                                                       226507, 150847, 227452,
                                                                       95802, 96352, 159097,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 239212, 0, 3,
                                                                       227452, 151522, 228397,
                                                                       96352, 96902, 159922,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 240367, 0, 3,
                                                                       228397, 152197, 229342,
                                                                       96902, 97452, 160747,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 241522, 0, 3,
                                                                       229342, 152872, 230287,
                                                                       97452, 98002, 161572,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 242677, 0, 3,
                                                                       231232, 154222, 232177,
                                                                       99102, 99652, 162397,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 243832, 0, 3,
                                                                       232177, 154897, 233122,
                                                                       99652, 100202, 163222,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 244987, 0, 3,
                                                                       233122, 155572, 234067,
                                                                       100202, 100752, 164047,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 246142, 0, 3,
                                                                       234067, 156247, 235012,
                                                                       100752, 101302, 164872,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 247297, 0, 3,
                                                                       235012, 156922, 235957,
                                                                       101302, 101852, 165697,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 248452, 0, 3,
                                                                       236902, 158272, 238057,
                                                                       102952, 103612, 166522,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 249838, 0, 3,
                                                                       238057, 159097, 239212,
                                                                       103612, 104272, 167512,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 251224, 0, 3,
                                                                       239212, 159922, 240367,
                                                                       104272, 104932, 168502,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 252610, 0, 3,
                                                                       240367, 160747, 241522,
                                                                       104932, 105592, 169492,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 253996, 0, 3,
                                                                       242677, 162397, 243832,
                                                                       106912, 107572, 170482,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 255382, 0, 3,
                                                                       243832, 163222, 244987,
                                                                       107572, 108232, 171472,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 256768, 0, 3,
                                                                       244987, 164047, 246142,
                                                                       108232, 108892, 172462,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 258154, 0, 3,
                                                                       246142, 164872, 247297,
                                                                       108892, 109552, 173452,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 259540, 0, 3,
                                                                       248452, 166522, 249838,
                                                                       110872, 111652, 174442,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 261178, 0, 3,
                                                                       249838, 167512, 251224,
                                                                       111652, 112432, 175612,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 262816, 0, 3,
                                                                       251224, 168502, 252610,
                                                                       112432, 113212, 176782,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 264454, 0, 3,
                                                                       253996, 170482, 255382,
                                                                       114772, 115552, 177952,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 266092, 0, 3,
                                                                       255382, 171472, 256768,
                                                                       115552, 116332, 179122,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 267730, 0, 3,
                                                                       256768, 172462, 258154,
                                                                       116332, 117112, 180292,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269368, 3, 118672,
                                                                       118687, 181504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269396, 3, 118687,
                                                                       118702, 181525, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269424, 3, 118702,
                                                                       118717, 181546, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269452, 3, 118717,
                                                                       118732, 181567, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269480, 3, 118732,
                                                                       118747, 181588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269508, 3, 118747,
                                                                       118762, 181609, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269536, 3, 118762,
                                                                       118777, 181630, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269564, 3, 118777,
                                                                       118792, 181651, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269592, 3, 118792,
                                                                       118807, 181672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269620, 3, 118807,
                                                                       118822, 181693, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269648, 3, 118822,
                                                                       118837, 181714, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269676, 3, 118837,
                                                                       118852, 181735, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269704, 3, 118882,
                                                                       118897, 181798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269732, 3, 118897,
                                                                       118912, 181819, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269760, 3, 118912,
                                                                       118927, 181840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269788, 3, 118927,
                                                                       118942, 181861, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269816, 3, 118942,
                                                                       118957, 181882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269844, 3, 118957,
                                                                       118972, 181903, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269872, 3, 118972,
                                                                       118987, 181924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269900, 3, 118987,
                                                                       119002, 181945, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269928, 3, 119002,
                                                                       119017, 181966, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269956, 3, 119017,
                                                                       119032, 181987, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269984, 3, 119032,
                                                                       119047, 182008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 270012, 3, 119047,
                                                                       119062, 182029, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270040, 0, 3,
                                                                       269368, 181504, 269396,
                                                                       119092, 119137, 182176,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270124, 0, 3,
                                                                       269396, 181525, 269424,
                                                                       119137, 119182, 182239,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270208, 0, 3,
                                                                       269424, 181546, 269452,
                                                                       119182, 119227, 182302,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270292, 0, 3,
                                                                       269452, 181567, 269480,
                                                                       119227, 119272, 182365,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270376, 0, 3,
                                                                       269480, 181588, 269508,
                                                                       119272, 119317, 182428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270460, 0, 3,
                                                                       269508, 181609, 269536,
                                                                       119317, 119362, 182491,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270544, 0, 3,
                                                                       269536, 181630, 269564,
                                                                       119362, 119407, 182554,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270628, 0, 3,
                                                                       269564, 181651, 269592,
                                                                       119407, 119452, 182617,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270712, 0, 3,
                                                                       269592, 181672, 269620,
                                                                       119452, 119497, 182680,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270796, 0, 3,
                                                                       269620, 181693, 269648,
                                                                       119497, 119542, 182743,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270880, 0, 3,
                                                                       269648, 181714, 269676,
                                                                       119542, 119587, 182806,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 270964, 0, 3,
                                                                       269704, 181798, 269732,
                                                                       119677, 119722, 182995,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271048, 0, 3,
                                                                       269732, 181819, 269760,
                                                                       119722, 119767, 183058,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271132, 0, 3,
                                                                       269760, 181840, 269788,
                                                                       119767, 119812, 183121,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271216, 0, 3,
                                                                       269788, 181861, 269816,
                                                                       119812, 119857, 183184,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271300, 0, 3,
                                                                       269816, 181882, 269844,
                                                                       119857, 119902, 183247,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271384, 0, 3,
                                                                       269844, 181903, 269872,
                                                                       119902, 119947, 183310,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271468, 0, 3,
                                                                       269872, 181924, 269900,
                                                                       119947, 119992, 183373,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271552, 0, 3,
                                                                       269900, 181945, 269928,
                                                                       119992, 120037, 183436,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271636, 0, 3,
                                                                       269928, 181966, 269956,
                                                                       120037, 120082, 183499,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271720, 0, 3,
                                                                       269956, 181987, 269984,
                                                                       120082, 120127, 183562,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 271804, 0, 3,
                                                                       269984, 182008, 270012,
                                                                       120127, 120172, 183625,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 271888, 0, 3,
                                                                       270040, 182176, 270124,
                                                                       120262, 120352, 183940,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 272056, 0, 3,
                                                                       270124, 182239, 270208,
                                                                       120352, 120442, 184066,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 272224, 0, 3,
                                                                       270208, 182302, 270292,
                                                                       120442, 120532, 184192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 272392, 0, 3,
                                                                       270292, 182365, 270376,
                                                                       120532, 120622, 184318,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 272560, 0, 3,
                                                                       270376, 182428, 270460,
                                                                       120622, 120712, 184444,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 272728, 0, 3,
                                                                       270460, 182491, 270544,
                                                                       120712, 120802, 184570,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 272896, 0, 3,
                                                                       270544, 182554, 270628,
                                                                       120802, 120892, 184696,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 273064, 0, 3,
                                                                       270628, 182617, 270712,
                                                                       120892, 120982, 184822,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 273232, 0, 3,
                                                                       270712, 182680, 270796,
                                                                       120982, 121072, 184948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 273400, 0, 3,
                                                                       270796, 182743, 270880,
                                                                       121072, 121162, 185074,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 273568, 0, 3,
                                                                       270964, 182995, 271048,
                                                                       121342, 121432, 185452,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 273736, 0, 3,
                                                                       271048, 183058, 271132,
                                                                       121432, 121522, 185578,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 273904, 0, 3,
                                                                       271132, 183121, 271216,
                                                                       121522, 121612, 185704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 274072, 0, 3,
                                                                       271216, 183184, 271300,
                                                                       121612, 121702, 185830,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 274240, 0, 3,
                                                                       271300, 183247, 271384,
                                                                       121702, 121792, 185956,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 274408, 0, 3,
                                                                       271384, 183310, 271468,
                                                                       121792, 121882, 186082,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 274576, 0, 3,
                                                                       271468, 183373, 271552,
                                                                       121882, 121972, 186208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 274744, 0, 3,
                                                                       271552, 183436, 271636,
                                                                       121972, 122062, 186334,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 274912, 0, 3,
                                                                       271636, 183499, 271720,
                                                                       122062, 122152, 186460,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 275080, 0, 3,
                                                                       271720, 183562, 271804,
                                                                       122152, 122242, 186586,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 275248, 0, 3,
                                                                       271888, 183940, 272056,
                                                                       122422, 122572, 187132,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 275528, 0, 3,
                                                                       272056, 184066, 272224,
                                                                       122572, 122722, 187342,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 275808, 0, 3,
                                                                       272224, 184192, 272392,
                                                                       122722, 122872, 187552,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 276088, 0, 3,
                                                                       272392, 184318, 272560,
                                                                       122872, 123022, 187762,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 276368, 0, 3,
                                                                       272560, 184444, 272728,
                                                                       123022, 123172, 187972,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 276648, 0, 3,
                                                                       272728, 184570, 272896,
                                                                       123172, 123322, 188182,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 276928, 0, 3,
                                                                       272896, 184696, 273064,
                                                                       123322, 123472, 188392,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 277208, 0, 3,
                                                                       273064, 184822, 273232,
                                                                       123472, 123622, 188602,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 277488, 0, 3,
                                                                       273232, 184948, 273400,
                                                                       123622, 123772, 188812,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 277768, 0, 3,
                                                                       273568, 185452, 273736,
                                                                       124072, 124222, 189442,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 278048, 0, 3,
                                                                       273736, 185578, 273904,
                                                                       124222, 124372, 189652,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 278328, 0, 3,
                                                                       273904, 185704, 274072,
                                                                       124372, 124522, 189862,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 278608, 0, 3,
                                                                       274072, 185830, 274240,
                                                                       124522, 124672, 190072,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 278888, 0, 3,
                                                                       274240, 185956, 274408,
                                                                       124672, 124822, 190282,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 279168, 0, 3,
                                                                       274408, 186082, 274576,
                                                                       124822, 124972, 190492,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 279448, 0, 3,
                                                                       274576, 186208, 274744,
                                                                       124972, 125122, 190702,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 279728, 0, 3,
                                                                       274744, 186334, 274912,
                                                                       125122, 125272, 190912,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 280008, 0, 3,
                                                                       274912, 186460, 275080,
                                                                       125272, 125422, 191122,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 280288, 0, 3,
                                                                       275248, 187132, 275528,
                                                                       125722, 125947, 191962,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 280708, 0, 3,
                                                                       275528, 187342, 275808,
                                                                       125947, 126172, 192277,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 281128, 0, 3,
                                                                       275808, 187552, 276088,
                                                                       126172, 126397, 192592,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 281548, 0, 3,
                                                                       276088, 187762, 276368,
                                                                       126397, 126622, 192907,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 281968, 0, 3,
                                                                       276368, 187972, 276648,
                                                                       126622, 126847, 193222,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 282388, 0, 3,
                                                                       276648, 188182, 276928,
                                                                       126847, 127072, 193537,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 282808, 0, 3,
                                                                       276928, 188392, 277208,
                                                                       127072, 127297, 193852,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 283228, 0, 3,
                                                                       277208, 188602, 277488,
                                                                       127297, 127522, 194167,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 283648, 0, 3,
                                                                       277768, 189442, 278048,
                                                                       127972, 128197, 195112,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 284068, 0, 3,
                                                                       278048, 189652, 278328,
                                                                       128197, 128422, 195427,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 284488, 0, 3,
                                                                       278328, 189862, 278608,
                                                                       128422, 128647, 195742,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 284908, 0, 3,
                                                                       278608, 190072, 278888,
                                                                       128647, 128872, 196057,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 285328, 0, 3,
                                                                       278888, 190282, 279168,
                                                                       128872, 129097, 196372,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 285748, 0, 3,
                                                                       279168, 190492, 279448,
                                                                       129097, 129322, 196687,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 286168, 0, 3,
                                                                       279448, 190702, 279728,
                                                                       129322, 129547, 197002,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 286588, 0, 3,
                                                                       279728, 190912, 280008,
                                                                       129547, 129772, 197317,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 287008, 0, 3,
                                                                       280288, 191962, 280708,
                                                                       130222, 130537, 198514,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 287596, 0, 3,
                                                                       280708, 192277, 281128,
                                                                       130537, 130852, 198955,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 288184, 0, 3,
                                                                       281128, 192592, 281548,
                                                                       130852, 131167, 199396,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 288772, 0, 3,
                                                                       281548, 192907, 281968,
                                                                       131167, 131482, 199837,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 289360, 0, 3,
                                                                       281968, 193222, 282388,
                                                                       131482, 131797, 200278,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 289948, 0, 3,
                                                                       282388, 193537, 282808,
                                                                       131797, 132112, 200719,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 290536, 0, 3,
                                                                       282808, 193852, 283228,
                                                                       132112, 132427, 201160,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 291124, 0, 3,
                                                                       283648, 195112, 284068,
                                                                       133057, 133372, 202483,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 291712, 0, 3,
                                                                       284068, 195427, 284488,
                                                                       133372, 133687, 202924,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 292300, 0, 3,
                                                                       284488, 195742, 284908,
                                                                       133687, 134002, 203365,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 292888, 0, 3,
                                                                       284908, 196057, 285328,
                                                                       134002, 134317, 203806,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 293476, 0, 3,
                                                                       285328, 196372, 285748,
                                                                       134317, 134632, 204247,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 294064, 0, 3,
                                                                       285748, 196687, 286168,
                                                                       134632, 134947, 204688,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 294652, 0, 3,
                                                                       286168, 197002, 286588,
                                                                       134947, 135262, 205129,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 295240, 0, 3,
                                                                       287008, 198514, 287596,
                                                                       135892, 136312, 206746,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 296024, 0, 3,
                                                                       287596, 198955, 288184,
                                                                       136312, 136732, 207334,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 296808, 0, 3,
                                                                       288184, 199396, 288772,
                                                                       136732, 137152, 207922,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 297592, 0, 3,
                                                                       288772, 199837, 289360,
                                                                       137152, 137572, 208510,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 298376, 0, 3,
                                                                       289360, 200278, 289948,
                                                                       137572, 137992, 209098,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 299160, 0, 3,
                                                                       289948, 200719, 290536,
                                                                       137992, 138412, 209686,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 299944, 0, 3,
                                                                       291124, 202483, 291712,
                                                                       139252, 139672, 211450,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 300728, 0, 3,
                                                                       291712, 202924, 292300,
                                                                       139672, 140092, 212038,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 301512, 0, 3,
                                                                       292300, 203365, 292888,
                                                                       140092, 140512, 212626,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 302296, 0, 3,
                                                                       292888, 203806, 293476,
                                                                       140512, 140932, 213214,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 303080, 0, 3,
                                                                       293476, 204247, 294064,
                                                                       140932, 141352, 213802,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 303864, 0, 3,
                                                                       294064, 204688, 294652,
                                                                       141352, 141772, 214390,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 304648, 0, 3,
                                                                       295240, 206746, 296024,
                                                                       142612, 143152, 216490,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 305656, 0, 3,
                                                                       296024, 207334, 296808,
                                                                       143152, 143692, 217246,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 306664, 0, 3,
                                                                       296808, 207922, 297592,
                                                                       143692, 144232, 218002,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 307672, 0, 3,
                                                                       297592, 208510, 298376,
                                                                       144232, 144772, 218758,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 308680, 0, 3,
                                                                       298376, 209098, 299160,
                                                                       144772, 145312, 219514,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 309688, 0, 3,
                                                                       299944, 211450, 300728,
                                                                       146392, 146932, 221782,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 310696, 0, 3,
                                                                       300728, 212038, 301512,
                                                                       146932, 147472, 222538,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 311704, 0, 3,
                                                                       301512, 212626, 302296,
                                                                       147472, 148012, 223294,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 312712, 0, 3,
                                                                       302296, 213214, 303080,
                                                                       148012, 148552, 224050,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 313720, 0, 3,
                                                                       303080, 213802, 303864,
                                                                       148552, 149092, 224806,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 314728, 0, 3,
                                                                       304648, 216490, 305656,
                                                                       150172, 150847, 227452,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 315988, 0, 3,
                                                                       305656, 217246, 306664,
                                                                       150847, 151522, 228397,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 317248, 0, 3,
                                                                       306664, 218002, 307672,
                                                                       151522, 152197, 229342,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 318508, 0, 3,
                                                                       307672, 218758, 308680,
                                                                       152197, 152872, 230287,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 319768, 0, 3,
                                                                       309688, 221782, 310696,
                                                                       154222, 154897, 233122,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 321028, 0, 3,
                                                                       310696, 222538, 311704,
                                                                       154897, 155572, 234067,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 322288, 0, 3,
                                                                       311704, 223294, 312712,
                                                                       155572, 156247, 235012,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 323548, 0, 3,
                                                                       312712, 224050, 313720,
                                                                       156247, 156922, 235957,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 324808, 0, 3,
                                                                       314728, 227452, 315988,
                                                                       158272, 159097, 239212,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 326348, 0, 3,
                                                                       315988, 228397, 317248,
                                                                       159097, 159922, 240367,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 327888, 0, 3,
                                                                       317248, 229342, 318508,
                                                                       159922, 160747, 241522,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 329428, 0, 3,
                                                                       319768, 233122, 321028,
                                                                       162397, 163222, 244987,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 330968, 0, 3,
                                                                       321028, 234067, 322288,
                                                                       163222, 164047, 246142,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 332508, 0, 3,
                                                                       322288, 235012, 323548,
                                                                       164047, 164872, 247297,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 334048, 0, 3,
                                                                       324808, 239212, 326348,
                                                                       166522, 167512, 251224,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 335896, 0, 3,
                                                                       326348, 240367, 327888,
                                                                       167512, 168502, 252610,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 337744, 0, 3,
                                                                       329428, 244987, 330968,
                                                                       170482, 171472, 256768,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 339592, 0, 3,
                                                                       330968, 246142, 332508,
                                                                       171472, 172462, 258154,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 341440, 0, 3,
                                                                       334048, 251224, 335896,
                                                                       174442, 175612, 262816,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 343624, 0, 3,
                                                                       337744, 256768, 339592,
                                                                       177952, 179122, 267730,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345808, 3, 181462,
                                                                       181483, 269368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345844, 3, 181483,
                                                                       181504, 269396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345880, 3, 181504,
                                                                       181525, 269424, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345916, 3, 181525,
                                                                       181546, 269452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345952, 3, 181546,
                                                                       181567, 269480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345988, 3, 181567,
                                                                       181588, 269508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346024, 3, 181588,
                                                                       181609, 269536, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346060, 3, 181609,
                                                                       181630, 269564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346096, 3, 181630,
                                                                       181651, 269592, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346132, 3, 181651,
                                                                       181672, 269620, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346168, 3, 181672,
                                                                       181693, 269648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346204, 3, 181693,
                                                                       181714, 269676, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346240, 3, 181756,
                                                                       181777, 269704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346276, 3, 181777,
                                                                       181798, 269732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346312, 3, 181798,
                                                                       181819, 269760, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346348, 3, 181819,
                                                                       181840, 269788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346384, 3, 181840,
                                                                       181861, 269816, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346420, 3, 181861,
                                                                       181882, 269844, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346456, 3, 181882,
                                                                       181903, 269872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346492, 3, 181903,
                                                                       181924, 269900, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346528, 3, 181924,
                                                                       181945, 269928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346564, 3, 181945,
                                                                       181966, 269956, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346600, 3, 181966,
                                                                       181987, 269984, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346636, 3, 181987,
                                                                       182008, 270012, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 346672, 0, 3,
                                                                       345808, 269368, 345844,
                                                                       182050, 182113, 270040,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 346780, 0, 3,
                                                                       345844, 269396, 345880,
                                                                       182113, 182176, 270124,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 346888, 0, 3,
                                                                       345880, 269424, 345916,
                                                                       182176, 182239, 270208,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 346996, 0, 3,
                                                                       345916, 269452, 345952,
                                                                       182239, 182302, 270292,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347104, 0, 3,
                                                                       345952, 269480, 345988,
                                                                       182302, 182365, 270376,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347212, 0, 3,
                                                                       345988, 269508, 346024,
                                                                       182365, 182428, 270460,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347320, 0, 3,
                                                                       346024, 269536, 346060,
                                                                       182428, 182491, 270544,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347428, 0, 3,
                                                                       346060, 269564, 346096,
                                                                       182491, 182554, 270628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347536, 0, 3,
                                                                       346096, 269592, 346132,
                                                                       182554, 182617, 270712,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347644, 0, 3,
                                                                       346132, 269620, 346168,
                                                                       182617, 182680, 270796,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347752, 0, 3,
                                                                       346168, 269648, 346204,
                                                                       182680, 182743, 270880,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347860, 0, 3,
                                                                       346240, 269704, 346276,
                                                                       182869, 182932, 270964,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 347968, 0, 3,
                                                                       346276, 269732, 346312,
                                                                       182932, 182995, 271048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348076, 0, 3,
                                                                       346312, 269760, 346348,
                                                                       182995, 183058, 271132,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348184, 0, 3,
                                                                       346348, 269788, 346384,
                                                                       183058, 183121, 271216,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348292, 0, 3,
                                                                       346384, 269816, 346420,
                                                                       183121, 183184, 271300,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348400, 0, 3,
                                                                       346420, 269844, 346456,
                                                                       183184, 183247, 271384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348508, 0, 3,
                                                                       346456, 269872, 346492,
                                                                       183247, 183310, 271468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348616, 0, 3,
                                                                       346492, 269900, 346528,
                                                                       183310, 183373, 271552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348724, 0, 3,
                                                                       346528, 269928, 346564,
                                                                       183373, 183436, 271636,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348832, 0, 3,
                                                                       346564, 269956, 346600,
                                                                       183436, 183499, 271720,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 348940, 0, 3,
                                                                       346600, 269984, 346636,
                                                                       183499, 183562, 271804,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 349048, 0, 3,
                                                                       346672, 270040, 346780,
                                                                       183688, 183814, 271888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 349264, 0, 3,
                                                                       346780, 270124, 346888,
                                                                       183814, 183940, 272056,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 349480, 0, 3,
                                                                       346888, 270208, 346996,
                                                                       183940, 184066, 272224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 349696, 0, 3,
                                                                       346996, 270292, 347104,
                                                                       184066, 184192, 272392,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 349912, 0, 3,
                                                                       347104, 270376, 347212,
                                                                       184192, 184318, 272560,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 350128, 0, 3,
                                                                       347212, 270460, 347320,
                                                                       184318, 184444, 272728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 350344, 0, 3,
                                                                       347320, 270544, 347428,
                                                                       184444, 184570, 272896,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 350560, 0, 3,
                                                                       347428, 270628, 347536,
                                                                       184570, 184696, 273064,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 350776, 0, 3,
                                                                       347536, 270712, 347644,
                                                                       184696, 184822, 273232,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 350992, 0, 3,
                                                                       347644, 270796, 347752,
                                                                       184822, 184948, 273400,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 351208, 0, 3,
                                                                       347860, 270964, 347968,
                                                                       185200, 185326, 273568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 351424, 0, 3,
                                                                       347968, 271048, 348076,
                                                                       185326, 185452, 273736,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 351640, 0, 3,
                                                                       348076, 271132, 348184,
                                                                       185452, 185578, 273904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 351856, 0, 3,
                                                                       348184, 271216, 348292,
                                                                       185578, 185704, 274072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 352072, 0, 3,
                                                                       348292, 271300, 348400,
                                                                       185704, 185830, 274240,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 352288, 0, 3,
                                                                       348400, 271384, 348508,
                                                                       185830, 185956, 274408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 352504, 0, 3,
                                                                       348508, 271468, 348616,
                                                                       185956, 186082, 274576,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 352720, 0, 3,
                                                                       348616, 271552, 348724,
                                                                       186082, 186208, 274744,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 352936, 0, 3,
                                                                       348724, 271636, 348832,
                                                                       186208, 186334, 274912,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 353152, 0, 3,
                                                                       348832, 271720, 348940,
                                                                       186334, 186460, 275080,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 353368, 0, 3,
                                                                       349048, 271888, 349264,
                                                                       186712, 186922, 275248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 353728, 0, 3,
                                                                       349264, 272056, 349480,
                                                                       186922, 187132, 275528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 354088, 0, 3,
                                                                       349480, 272224, 349696,
                                                                       187132, 187342, 275808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 354448, 0, 3,
                                                                       349696, 272392, 349912,
                                                                       187342, 187552, 276088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 354808, 0, 3,
                                                                       349912, 272560, 350128,
                                                                       187552, 187762, 276368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 355168, 0, 3,
                                                                       350128, 272728, 350344,
                                                                       187762, 187972, 276648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 355528, 0, 3,
                                                                       350344, 272896, 350560,
                                                                       187972, 188182, 276928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 355888, 0, 3,
                                                                       350560, 273064, 350776,
                                                                       188182, 188392, 277208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 356248, 0, 3,
                                                                       350776, 273232, 350992,
                                                                       188392, 188602, 277488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 356608, 0, 3,
                                                                       351208, 273568, 351424,
                                                                       189022, 189232, 277768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 356968, 0, 3,
                                                                       351424, 273736, 351640,
                                                                       189232, 189442, 278048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 357328, 0, 3,
                                                                       351640, 273904, 351856,
                                                                       189442, 189652, 278328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 357688, 0, 3,
                                                                       351856, 274072, 352072,
                                                                       189652, 189862, 278608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 358048, 0, 3,
                                                                       352072, 274240, 352288,
                                                                       189862, 190072, 278888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 358408, 0, 3,
                                                                       352288, 274408, 352504,
                                                                       190072, 190282, 279168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 358768, 0, 3,
                                                                       352504, 274576, 352720,
                                                                       190282, 190492, 279448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 359128, 0, 3,
                                                                       352720, 274744, 352936,
                                                                       190492, 190702, 279728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 359488, 0, 3,
                                                                       352936, 274912, 353152,
                                                                       190702, 190912, 280008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 359848, 0, 3,
                                                                       353368, 275248, 353728,
                                                                       191332, 191647, 280288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 360388, 0, 3,
                                                                       353728, 275528, 354088,
                                                                       191647, 191962, 280708,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 360928, 0, 3,
                                                                       354088, 275808, 354448,
                                                                       191962, 192277, 281128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 361468, 0, 3,
                                                                       354448, 276088, 354808,
                                                                       192277, 192592, 281548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 362008, 0, 3,
                                                                       354808, 276368, 355168,
                                                                       192592, 192907, 281968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 362548, 0, 3,
                                                                       355168, 276648, 355528,
                                                                       192907, 193222, 282388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 363088, 0, 3,
                                                                       355528, 276928, 355888,
                                                                       193222, 193537, 282808,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 363628, 0, 3,
                                                                       355888, 277208, 356248,
                                                                       193537, 193852, 283228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 364168, 0, 3,
                                                                       356608, 277768, 356968,
                                                                       194482, 194797, 283648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 364708, 0, 3,
                                                                       356968, 278048, 357328,
                                                                       194797, 195112, 284068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 365248, 0, 3,
                                                                       357328, 278328, 357688,
                                                                       195112, 195427, 284488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 365788, 0, 3,
                                                                       357688, 278608, 358048,
                                                                       195427, 195742, 284908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 366328, 0, 3,
                                                                       358048, 278888, 358408,
                                                                       195742, 196057, 285328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 366868, 0, 3,
                                                                       358408, 279168, 358768,
                                                                       196057, 196372, 285748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 367408, 0, 3,
                                                                       358768, 279448, 359128,
                                                                       196372, 196687, 286168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 367948, 0, 3,
                                                                       359128, 279728, 359488,
                                                                       196687, 197002, 286588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 368488, 0, 3,
                                                                       359848, 280288, 360388,
                                                                       197632, 198073, 287008,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 369244, 0, 3,
                                                                       360388, 280708, 360928,
                                                                       198073, 198514, 287596,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 370000, 0, 3,
                                                                       360928, 281128, 361468,
                                                                       198514, 198955, 288184,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 370756, 0, 3,
                                                                       361468, 281548, 362008,
                                                                       198955, 199396, 288772,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 371512, 0, 3,
                                                                       362008, 281968, 362548,
                                                                       199396, 199837, 289360,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 372268, 0, 3,
                                                                       362548, 282388, 363088,
                                                                       199837, 200278, 289948,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 373024, 0, 3,
                                                                       363088, 282808, 363628,
                                                                       200278, 200719, 290536,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 373780, 0, 3,
                                                                       364168, 283648, 364708,
                                                                       201601, 202042, 291124,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 374536, 0, 3,
                                                                       364708, 284068, 365248,
                                                                       202042, 202483, 291712,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 375292, 0, 3,
                                                                       365248, 284488, 365788,
                                                                       202483, 202924, 292300,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 376048, 0, 3,
                                                                       365788, 284908, 366328,
                                                                       202924, 203365, 292888,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 376804, 0, 3,
                                                                       366328, 285328, 366868,
                                                                       203365, 203806, 293476,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 377560, 0, 3,
                                                                       366868, 285748, 367408,
                                                                       203806, 204247, 294064,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 378316, 0, 3,
                                                                       367408, 286168, 367948,
                                                                       204247, 204688, 294652,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 379072, 0, 3,
                                                                       368488, 287008, 369244,
                                                                       205570, 206158, 295240,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 380080, 0, 3,
                                                                       369244, 287596, 370000,
                                                                       206158, 206746, 296024,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 381088, 0, 3,
                                                                       370000, 288184, 370756,
                                                                       206746, 207334, 296808,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 382096, 0, 3,
                                                                       370756, 288772, 371512,
                                                                       207334, 207922, 297592,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 383104, 0, 3,
                                                                       371512, 289360, 372268,
                                                                       207922, 208510, 298376,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 384112, 0, 3,
                                                                       372268, 289948, 373024,
                                                                       208510, 209098, 299160,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 385120, 0, 3,
                                                                       373780, 291124, 374536,
                                                                       210274, 210862, 299944,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 386128, 0, 3,
                                                                       374536, 291712, 375292,
                                                                       210862, 211450, 300728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 387136, 0, 3,
                                                                       375292, 292300, 376048,
                                                                       211450, 212038, 301512,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 388144, 0, 3,
                                                                       376048, 292888, 376804,
                                                                       212038, 212626, 302296,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 389152, 0, 3,
                                                                       376804, 293476, 377560,
                                                                       212626, 213214, 303080,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 390160, 0, 3,
                                                                       377560, 294064, 378316,
                                                                       213214, 213802, 303864,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 391168, 0, 3,
                                                                       379072, 295240, 380080,
                                                                       214978, 215734, 304648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 392464, 0, 3,
                                                                       380080, 296024, 381088,
                                                                       215734, 216490, 305656,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 393760, 0, 3,
                                                                       381088, 296808, 382096,
                                                                       216490, 217246, 306664,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 395056, 0, 3,
                                                                       382096, 297592, 383104,
                                                                       217246, 218002, 307672,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 396352, 0, 3,
                                                                       383104, 298376, 384112,
                                                                       218002, 218758, 308680,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 397648, 0, 3,
                                                                       385120, 299944, 386128,
                                                                       220270, 221026, 309688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 398944, 0, 3,
                                                                       386128, 300728, 387136,
                                                                       221026, 221782, 310696,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 400240, 0, 3,
                                                                       387136, 301512, 388144,
                                                                       221782, 222538, 311704,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 401536, 0, 3,
                                                                       388144, 302296, 389152,
                                                                       222538, 223294, 312712,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 402832, 0, 3,
                                                                       389152, 303080, 390160,
                                                                       223294, 224050, 313720,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 404128, 0, 3,
                                                                       391168, 304648, 392464,
                                                                       225562, 226507, 314728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 405748, 0, 3,
                                                                       392464, 305656, 393760,
                                                                       226507, 227452, 315988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 407368, 0, 3,
                                                                       393760, 306664, 395056,
                                                                       227452, 228397, 317248,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 408988, 0, 3,
                                                                       395056, 307672, 396352,
                                                                       228397, 229342, 318508,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 410608, 0, 3,
                                                                       397648, 309688, 398944,
                                                                       231232, 232177, 319768,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 412228, 0, 3,
                                                                       398944, 310696, 400240,
                                                                       232177, 233122, 321028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 413848, 0, 3,
                                                                       400240, 311704, 401536,
                                                                       233122, 234067, 322288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 415468, 0, 3,
                                                                       401536, 312712, 402832,
                                                                       234067, 235012, 323548,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 417088, 0, 3,
                                                                       404128, 314728, 405748,
                                                                       236902, 238057, 324808,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 419068, 0, 3,
                                                                       405748, 315988, 407368,
                                                                       238057, 239212, 326348,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 421048, 0, 3,
                                                                       407368, 317248, 408988,
                                                                       239212, 240367, 327888,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 423028, 0, 3,
                                                                       410608, 319768, 412228,
                                                                       242677, 243832, 329428,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 425008, 0, 3,
                                                                       412228, 321028, 413848,
                                                                       243832, 244987, 330968,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 426988, 0, 3,
                                                                       413848, 322288, 415468,
                                                                       244987, 246142, 332508,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 428968, 0, 3,
                                                                       417088, 324808, 419068,
                                                                       248452, 249838, 334048,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 431344, 0, 3,
                                                                       419068, 326348, 421048,
                                                                       249838, 251224, 335896,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 433720, 0, 3,
                                                                       423028, 329428, 425008,
                                                                       253996, 255382, 337744,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 436096, 0, 3,
                                                                       425008, 330968, 426988,
                                                                       255382, 256768, 339592,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 438472, 0, 3,
                                                                       428968, 334048, 431344,
                                                                       259540, 261178, 341440,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 441280, 0, 3,
                                                                       433720, 337744, 436096,
                                                                       264454, 266092, 343624,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 444088, 379072, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 445516, 385120, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 446944, 391168, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 448780, 397648, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 450616, 404128, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 452911, 410608, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 455206, 417088, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 458011, 423028, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 460816, 428968, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 464182, 433720, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 467548, 438472, 2808, ncols);

                    simdfunc::contract_primitives(buffer, 471526, 441280, 2808, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 445096, 444088, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 446524, 445516, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 448240, 446944, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 450076, 448780, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 452236, 450616, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 454531, 452911, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 457186, 455206, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 459991, 458011, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 463192, 460816, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 466558, 464182, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 470356, 467548, 78, 1, nmax);

        simdtrf::transform_k_inner(buffer, 474334, 471526, 78, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 475504, 445096, 448240, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 476764, 446524, 450076, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 478024, 448240, 452236, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 479644, 450076, 454531, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 481264, 452236, 457186, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 483289, 454531, 459991, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 485314, 457186, 463192, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 487789, 459991, 466558, 15,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 490264, 463192, 470356, 15,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 493234, 466558, 474334, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 496204, 475504, 478024, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 498724, 476764, 479644, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 501244, 478024, 481264, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 504484, 479644, 483289, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 507724, 481264, 485314, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 511774, 483289, 487789, 15,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 515824, 485314, 490264, 15,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 520774, 487789, 493234, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 525724, 496204, 501244, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 529924, 498724, 504484, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 534124, 501244, 507724, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 539524, 504484, 511774, 15,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 544924, 507724, 515824, 15,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 551674, 511774, 520774, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 558424, 525724, 534124, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 564724, 529924, 539524, 15,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 571024, 534124, 544924, 15,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 579124, 539524, 551674, 15,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 587224, 558424, 571024, 15,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 596044, 564724, 579124, 15,
                                             nmax);

        simdtrf::transform_h_inner(buffer, 604864, 596044, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 604864, 165, nmax);

        simdtrf::transform_h_inner(buffer, 604864, 587224, 28, 15, nmax);

        simdtrf::transform_i_outer(values + 2145 * nvalues + n * npairs, nvalues, buffer, 604864,
                                   165, nmax);
    }

    for (size_t m = 0; m < 4290; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
