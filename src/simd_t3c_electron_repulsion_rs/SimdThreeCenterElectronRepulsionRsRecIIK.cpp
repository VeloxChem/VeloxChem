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


#include "SimdThreeCenterElectronRepulsionRsRecIIK.hpp"

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
compute_rs_iik_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 859776, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 5070 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 859776, 582408, 39333, dimensions);

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
                                                            15, 16, 17, 18, 19}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 26, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17, 18, 19}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 109, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 115, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 121, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 127, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 139, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 41, 42,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 145, 0, 3, 42, 43,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 43, 44,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 151, 0, 3, 44, 45,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 3, 7, 8,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 3, 8, 9,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 3, 9, 10,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 10, 11,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 11, 12,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 184, 0, 3, 12, 13,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 190, 0, 3, 13, 14,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 196, 0, 3, 14, 15,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 15, 16,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 16, 17,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 214, 0, 3, 17, 18,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 220, 0, 3, 18, 19,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 226, 0, 3, 19, 20,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 20, 21,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 21, 22,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 244, 0, 3, 22, 23,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 250, 0, 3, 23, 24,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 256, 0, 3, 27, 28,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 28, 29,
                                                                       103, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 29, 30,
                                                                       106, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 274, 0, 3, 30, 31,
                                                                       109, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 280, 0, 3, 31, 32,
                                                                       112, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 286, 0, 3, 32, 33,
                                                                       115, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 292, 0, 3, 33, 34,
                                                                       118, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 34, 35,
                                                                       121, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 304, 0, 3, 35, 36,
                                                                       124, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 310, 0, 3, 36, 37,
                                                                       127, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 316, 0, 3, 37, 38,
                                                                       130, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 322, 0, 3, 38, 39,
                                                                       133, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 39, 40,
                                                                       136, 139, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 334, 0, 3, 40, 41,
                                                                       139, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 340, 0, 3, 41, 42,
                                                                       142, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 346, 0, 3, 42, 43,
                                                                       145, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 352, 0, 3, 43, 44,
                                                                       148, 151, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 46, 49,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 49, 52,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 52, 55,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 55, 58,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 58, 61,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 61, 64,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 64, 67,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 67, 70,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 70, 73,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 73, 76,
                                                                       208, 214, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 76, 79,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 79, 82,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 82, 85,
                                                                       226, 232, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 85, 88,
                                                                       232, 238, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 88, 91,
                                                                       238, 244, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 91, 94,
                                                                       244, 250, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 100,
                                                                       103, 256, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 103,
                                                                       106, 262, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 106,
                                                                       109, 268, 274, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 109,
                                                                       112, 274, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 112,
                                                                       115, 280, 286, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 115,
                                                                       118, 286, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 118,
                                                                       121, 292, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 121,
                                                                       124, 298, 304, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 598, 0, 3, 124,
                                                                       127, 304, 310, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 608, 0, 3, 127,
                                                                       130, 310, 316, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 618, 0, 3, 130,
                                                                       133, 316, 322, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 628, 0, 3, 133,
                                                                       136, 322, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 136,
                                                                       139, 328, 334, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 648, 0, 3, 139,
                                                                       142, 334, 340, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 658, 0, 3, 142,
                                                                       145, 340, 346, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 668, 0, 3, 145,
                                                                       148, 346, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 678, 0, 3, 154,
                                                                       160, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 160,
                                                                       166, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 708, 0, 3, 166,
                                                                       172, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 723, 0, 3, 172,
                                                                       178, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 738, 0, 3, 178,
                                                                       184, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 753, 0, 3, 184,
                                                                       190, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 768, 0, 3, 190,
                                                                       196, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 783, 0, 3, 196,
                                                                       202, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 798, 0, 3, 202,
                                                                       208, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 813, 0, 3, 208,
                                                                       214, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 828, 0, 3, 214,
                                                                       220, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 843, 0, 3, 220,
                                                                       226, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 858, 0, 3, 226,
                                                                       232, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 873, 0, 3, 232,
                                                                       238, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 888, 0, 3, 238,
                                                                       244, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 903, 0, 3, 256,
                                                                       262, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 262,
                                                                       268, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 933, 0, 3, 268,
                                                                       274, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 948, 0, 3, 274,
                                                                       280, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 963, 0, 3, 280,
                                                                       286, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 978, 0, 3, 286,
                                                                       292, 568, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 993, 0, 3, 292,
                                                                       298, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1008, 0, 3, 298,
                                                                       304, 588, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 304,
                                                                       310, 598, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1038, 0, 3, 310,
                                                                       316, 608, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1053, 0, 3, 316,
                                                                       322, 618, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1068, 0, 3, 322,
                                                                       328, 628, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1083, 0, 3, 328,
                                                                       334, 638, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1098, 0, 3, 334,
                                                                       340, 648, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1113, 0, 3, 340,
                                                                       346, 658, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 358,
                                                                       368, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1149, 0, 3, 368,
                                                                       378, 693, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 378,
                                                                       388, 708, 723, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1191, 0, 3, 388,
                                                                       398, 723, 738, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 398,
                                                                       408, 738, 753, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1233, 0, 3, 408,
                                                                       418, 753, 768, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 418,
                                                                       428, 768, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1275, 0, 3, 428,
                                                                       438, 783, 798, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 438,
                                                                       448, 798, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1317, 0, 3, 448,
                                                                       458, 813, 828, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1338, 0, 3, 458,
                                                                       468, 828, 843, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1359, 0, 3, 468,
                                                                       478, 843, 858, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 478,
                                                                       488, 858, 873, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1401, 0, 3, 488,
                                                                       498, 873, 888, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1422, 0, 3, 518,
                                                                       528, 903, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1443, 0, 3, 528,
                                                                       538, 918, 933, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 538,
                                                                       548, 933, 948, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1485, 0, 3, 548,
                                                                       558, 948, 963, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1506, 0, 3, 558,
                                                                       568, 963, 978, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1527, 0, 3, 568,
                                                                       578, 978, 993, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 578,
                                                                       588, 993, 1008, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1569, 0, 3, 588,
                                                                       598, 1008, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1590, 0, 3, 598,
                                                                       608, 1023, 1038, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1611, 0, 3, 608,
                                                                       618, 1038, 1053, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 618,
                                                                       628, 1053, 1068, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1653, 0, 3, 628,
                                                                       638, 1068, 1083, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1674, 0, 3, 638,
                                                                       648, 1083, 1098, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1695, 0, 3, 648,
                                                                       658, 1098, 1113, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 678,
                                                                       693, 1128, 1149, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 693,
                                                                       708, 1149, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 708,
                                                                       723, 1170, 1191, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 723,
                                                                       738, 1191, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 738,
                                                                       753, 1212, 1233, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 753,
                                                                       768, 1233, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 768,
                                                                       783, 1254, 1275, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 783,
                                                                       798, 1275, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 798,
                                                                       813, 1296, 1317, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 813,
                                                                       828, 1317, 1338, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 828,
                                                                       843, 1338, 1359, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 843,
                                                                       858, 1359, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 858,
                                                                       873, 1380, 1401, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 903,
                                                                       918, 1422, 1443, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 918,
                                                                       933, 1443, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2136, 0, 3, 933,
                                                                       948, 1464, 1485, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 948,
                                                                       963, 1485, 1506, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 963,
                                                                       978, 1506, 1527, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2220, 0, 3, 978,
                                                                       993, 1527, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 993,
                                                                       1008, 1548, 1569, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1008,
                                                                       1023, 1569, 1590, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2304, 0, 3, 1023,
                                                                       1038, 1590, 1611, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2332, 0, 3, 1038,
                                                                       1053, 1611, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2360, 0, 3, 1053,
                                                                       1068, 1632, 1653, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2388, 0, 3, 1068,
                                                                       1083, 1653, 1674, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 2416, 0, 3, 1083,
                                                                       1098, 1674, 1695, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 1128,
                                                                       1149, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2480, 0, 3, 1149,
                                                                       1170, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 1170,
                                                                       1191, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 1191,
                                                                       1212, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 1212,
                                                                       1233, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2624, 0, 3, 1233,
                                                                       1254, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2660, 0, 3, 1254,
                                                                       1275, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2696, 0, 3, 1275,
                                                                       1296, 1912, 1940, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2732, 0, 3, 1296,
                                                                       1317, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2768, 0, 3, 1317,
                                                                       1338, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2804, 0, 3, 1338,
                                                                       1359, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2840, 0, 3, 1359,
                                                                       1380, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2876, 0, 3, 1422,
                                                                       1443, 2080, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2912, 0, 3, 1443,
                                                                       1464, 2108, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 1464,
                                                                       1485, 2136, 2164, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2984, 0, 3, 1485,
                                                                       1506, 2164, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3020, 0, 3, 1506,
                                                                       1527, 2192, 2220, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3056, 0, 3, 1527,
                                                                       1548, 2220, 2248, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3092, 0, 3, 1548,
                                                                       1569, 2248, 2276, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 1569,
                                                                       1590, 2276, 2304, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3164, 0, 3, 1590,
                                                                       1611, 2304, 2332, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3200, 0, 3, 1611,
                                                                       1632, 2332, 2360, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3236, 0, 3, 1632,
                                                                       1653, 2360, 2388, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 3272, 0, 3, 1653,
                                                                       1674, 2388, 2416, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 1716,
                                                                       1744, 2444, 2480, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3353, 0, 3, 1744,
                                                                       1772, 2480, 2516, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3398, 0, 3, 1772,
                                                                       1800, 2516, 2552, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3443, 0, 3, 1800,
                                                                       1828, 2552, 2588, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3488, 0, 3, 1828,
                                                                       1856, 2588, 2624, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3533, 0, 3, 1856,
                                                                       1884, 2624, 2660, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3578, 0, 3, 1884,
                                                                       1912, 2660, 2696, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3623, 0, 3, 1912,
                                                                       1940, 2696, 2732, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3668, 0, 3, 1940,
                                                                       1968, 2732, 2768, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3713, 0, 3, 1968,
                                                                       1996, 2768, 2804, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3758, 0, 3, 1996,
                                                                       2024, 2804, 2840, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3803, 0, 3, 2080,
                                                                       2108, 2876, 2912, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3848, 0, 3, 2108,
                                                                       2136, 2912, 2948, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3893, 0, 3, 2136,
                                                                       2164, 2948, 2984, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3938, 0, 3, 2164,
                                                                       2192, 2984, 3020, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 3983, 0, 3, 2192,
                                                                       2220, 3020, 3056, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4028, 0, 3, 2220,
                                                                       2248, 3056, 3092, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4073, 0, 3, 2248,
                                                                       2276, 3092, 3128, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4118, 0, 3, 2276,
                                                                       2304, 3128, 3164, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4163, 0, 3, 2304,
                                                                       2332, 3164, 3200, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4208, 0, 3, 2332,
                                                                       2360, 3200, 3236, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 4253, 0, 3, 2360,
                                                                       2388, 3236, 3272, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 2444,
                                                                       2480, 3308, 3353, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4353, 0, 3, 2480,
                                                                       2516, 3353, 3398, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2516,
                                                                       2552, 3398, 3443, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4463, 0, 3, 2552,
                                                                       2588, 3443, 3488, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4518, 0, 3, 2588,
                                                                       2624, 3488, 3533, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4573, 0, 3, 2624,
                                                                       2660, 3533, 3578, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 2660,
                                                                       2696, 3578, 3623, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4683, 0, 3, 2696,
                                                                       2732, 3623, 3668, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4738, 0, 3, 2732,
                                                                       2768, 3668, 3713, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4793, 0, 3, 2768,
                                                                       2804, 3713, 3758, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 2876,
                                                                       2912, 3803, 3848, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4903, 0, 3, 2912,
                                                                       2948, 3848, 3893, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 2948,
                                                                       2984, 3893, 3938, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5013, 0, 3, 2984,
                                                                       3020, 3938, 3983, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5068, 0, 3, 3020,
                                                                       3056, 3983, 4028, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5123, 0, 3, 3056,
                                                                       3092, 4028, 4073, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5178, 0, 3, 3092,
                                                                       3128, 4073, 4118, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5233, 0, 3, 3128,
                                                                       3164, 4118, 4163, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 3164,
                                                                       3200, 4163, 4208, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 5343, 0, 3, 3200,
                                                                       3236, 4208, 4253, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5398, 0, 3, 3308,
                                                                       3353, 4298, 4353, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5464, 0, 3, 3353,
                                                                       3398, 4353, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5530, 0, 3, 3398,
                                                                       3443, 4408, 4463, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5596, 0, 3, 3443,
                                                                       3488, 4463, 4518, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5662, 0, 3, 3488,
                                                                       3533, 4518, 4573, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5728, 0, 3, 3533,
                                                                       3578, 4573, 4628, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5794, 0, 3, 3578,
                                                                       3623, 4628, 4683, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5860, 0, 3, 3623,
                                                                       3668, 4683, 4738, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5926, 0, 3, 3668,
                                                                       3713, 4738, 4793, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 5992, 0, 3, 3803,
                                                                       3848, 4848, 4903, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6058, 0, 3, 3848,
                                                                       3893, 4903, 4958, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6124, 0, 3, 3893,
                                                                       3938, 4958, 5013, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6190, 0, 3, 3938,
                                                                       3983, 5013, 5068, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6256, 0, 3, 3983,
                                                                       4028, 5068, 5123, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6322, 0, 3, 4028,
                                                                       4073, 5123, 5178, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6388, 0, 3, 4073,
                                                                       4118, 5178, 5233, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6454, 0, 3, 4118,
                                                                       4163, 5233, 5288, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 6520, 0, 3, 4163,
                                                                       4208, 5288, 5343, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6586, 0, 3, 4298,
                                                                       4353, 5398, 5464, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6664, 0, 3, 4353,
                                                                       4408, 5464, 5530, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6742, 0, 3, 4408,
                                                                       4463, 5530, 5596, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6820, 0, 3, 4463,
                                                                       4518, 5596, 5662, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6898, 0, 3, 4518,
                                                                       4573, 5662, 5728, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 6976, 0, 3, 4573,
                                                                       4628, 5728, 5794, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7054, 0, 3, 4628,
                                                                       4683, 5794, 5860, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7132, 0, 3, 4683,
                                                                       4738, 5860, 5926, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7210, 0, 3, 4848,
                                                                       4903, 5992, 6058, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7288, 0, 3, 4903,
                                                                       4958, 6058, 6124, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7366, 0, 3, 4958,
                                                                       5013, 6124, 6190, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7444, 0, 3, 5013,
                                                                       5068, 6190, 6256, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7522, 0, 3, 5068,
                                                                       5123, 6256, 6322, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7600, 0, 3, 5123,
                                                                       5178, 6322, 6388, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7678, 0, 3, 5178,
                                                                       5233, 6388, 6454, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 7756, 0, 3, 5233,
                                                                       5288, 6454, 6520, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 7834, 0, 3, 5398,
                                                                       5464, 6586, 6664, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 7925, 0, 3, 5464,
                                                                       5530, 6664, 6742, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8016, 0, 3, 5530,
                                                                       5596, 6742, 6820, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8107, 0, 3, 5596,
                                                                       5662, 6820, 6898, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8198, 0, 3, 5662,
                                                                       5728, 6898, 6976, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8289, 0, 3, 5728,
                                                                       5794, 6976, 7054, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8380, 0, 3, 5794,
                                                                       5860, 7054, 7132, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8471, 0, 3, 5992,
                                                                       6058, 7210, 7288, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8562, 0, 3, 6058,
                                                                       6124, 7288, 7366, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8653, 0, 3, 6124,
                                                                       6190, 7366, 7444, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8744, 0, 3, 6190,
                                                                       6256, 7444, 7522, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8835, 0, 3, 6256,
                                                                       6322, 7522, 7600, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 8926, 0, 3, 6322,
                                                                       6388, 7600, 7678, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 9017, 0, 3, 6388,
                                                                       6454, 7678, 7756, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9108, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9111, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9114, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9117, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9120, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9123, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9126, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9129, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9132, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9135, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9138, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9141, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9144, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9147, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9150, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9153, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9156, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9159, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9162, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9165, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9168, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9171, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9174, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9177, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9180, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9183, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9186, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9189, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9192, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9195, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9198, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9201, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9204, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9207, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9210, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9213, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9216, 3, 44,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9219, 3, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9222, 3, 7, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9231, 3, 8, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9240, 3, 9, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9249, 3, 10, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9258, 3, 11, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9267, 3, 12, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9276, 3, 13, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9285, 3, 14, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9294, 3, 15, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9303, 3, 16, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9312, 3, 17, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9321, 3, 18, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9330, 3, 19, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9339, 3, 20, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9348, 3, 21, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9357, 3, 22, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9366, 3, 23, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9375, 3, 24, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9384, 3, 27, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9393, 3, 28, 103,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9402, 3, 29, 106,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9411, 3, 30, 109,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9420, 3, 31, 112,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9429, 3, 32, 115,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9438, 3, 33, 118,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9447, 3, 34, 121,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9456, 3, 35, 124,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9465, 3, 36, 127,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9474, 3, 37, 130,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9483, 3, 38, 133,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9492, 3, 39, 136,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9501, 3, 40, 139,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9510, 3, 41, 142,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9519, 3, 42, 145,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9528, 3, 43, 148,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9537, 3, 44, 151,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9546, 3, 46, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9564, 3, 49, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9582, 3, 52, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9600, 3, 55, 172,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9618, 3, 58, 178,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9636, 3, 61, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9654, 3, 64, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9672, 3, 67, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9690, 3, 70, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9708, 3, 73, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9726, 3, 76, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9744, 3, 79, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9762, 3, 82, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9780, 3, 85, 232,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9798, 3, 88, 238,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9816, 3, 91, 244,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9834, 3, 94, 250,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9852, 3, 100, 256,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9870, 3, 103, 262,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9888, 3, 106, 268,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9906, 3, 109, 274,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9924, 3, 112, 280,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9942, 3, 115, 286,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9960, 3, 118, 292,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9978, 3, 121, 298,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9996, 3, 124, 304,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10014, 3, 127,
                                                                       310, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10032, 3, 130,
                                                                       316, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10050, 3, 133,
                                                                       322, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10068, 3, 136,
                                                                       328, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10086, 3, 139,
                                                                       334, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10104, 3, 142,
                                                                       340, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10122, 3, 145,
                                                                       346, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10140, 3, 148,
                                                                       352, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10158, 3, 154,
                                                                       358, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10188, 3, 160,
                                                                       368, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10218, 3, 166,
                                                                       378, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10248, 3, 172,
                                                                       388, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10278, 3, 178,
                                                                       398, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10308, 3, 184,
                                                                       408, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10338, 3, 190,
                                                                       418, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10368, 3, 196,
                                                                       428, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10398, 3, 202,
                                                                       438, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10428, 3, 208,
                                                                       448, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10458, 3, 214,
                                                                       458, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10488, 3, 220,
                                                                       468, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10518, 3, 226,
                                                                       478, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10548, 3, 232,
                                                                       488, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10578, 3, 238,
                                                                       498, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10608, 3, 244,
                                                                       508, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10638, 3, 256,
                                                                       518, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10668, 3, 262,
                                                                       528, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10698, 3, 268,
                                                                       538, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10728, 3, 274,
                                                                       548, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10758, 3, 280,
                                                                       558, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10788, 3, 286,
                                                                       568, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10818, 3, 292,
                                                                       578, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10848, 3, 298,
                                                                       588, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10878, 3, 304,
                                                                       598, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10908, 3, 310,
                                                                       608, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10938, 3, 316,
                                                                       618, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10968, 3, 322,
                                                                       628, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10998, 3, 328,
                                                                       638, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 11028, 3, 334,
                                                                       648, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 11058, 3, 340,
                                                                       658, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 11088, 3, 346,
                                                                       668, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11118, 3, 358,
                                                                       678, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11163, 3, 368,
                                                                       693, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11208, 3, 378,
                                                                       708, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11253, 3, 388,
                                                                       723, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11298, 3, 398,
                                                                       738, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11343, 3, 408,
                                                                       753, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11388, 3, 418,
                                                                       768, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11433, 3, 428,
                                                                       783, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11478, 3, 438,
                                                                       798, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11523, 3, 448,
                                                                       813, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11568, 3, 458,
                                                                       828, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11613, 3, 468,
                                                                       843, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11658, 3, 478,
                                                                       858, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11703, 3, 488,
                                                                       873, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11748, 3, 498,
                                                                       888, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11793, 3, 518,
                                                                       903, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11838, 3, 528,
                                                                       918, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11883, 3, 538,
                                                                       933, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11928, 3, 548,
                                                                       948, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11973, 3, 558,
                                                                       963, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12018, 3, 568,
                                                                       978, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12063, 3, 578,
                                                                       993, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12108, 3, 588,
                                                                       1008, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12153, 3, 598,
                                                                       1023, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12198, 3, 608,
                                                                       1038, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12243, 3, 618,
                                                                       1053, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12288, 3, 628,
                                                                       1068, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12333, 3, 638,
                                                                       1083, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12378, 3, 648,
                                                                       1098, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12423, 3, 658,
                                                                       1113, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12468, 3, 678,
                                                                       1128, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12531, 3, 693,
                                                                       1149, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12594, 3, 708,
                                                                       1170, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12657, 3, 723,
                                                                       1191, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12720, 3, 738,
                                                                       1212, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12783, 3, 753,
                                                                       1233, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12846, 3, 768,
                                                                       1254, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12909, 3, 783,
                                                                       1275, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12972, 3, 798,
                                                                       1296, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13035, 3, 813,
                                                                       1317, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13098, 3, 828,
                                                                       1338, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13161, 3, 843,
                                                                       1359, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13224, 3, 858,
                                                                       1380, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13287, 3, 873,
                                                                       1401, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13350, 3, 903,
                                                                       1422, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13413, 3, 918,
                                                                       1443, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13476, 3, 933,
                                                                       1464, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13539, 3, 948,
                                                                       1485, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13602, 3, 963,
                                                                       1506, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13665, 3, 978,
                                                                       1527, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13728, 3, 993,
                                                                       1548, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13791, 3, 1008,
                                                                       1569, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13854, 3, 1023,
                                                                       1590, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13917, 3, 1038,
                                                                       1611, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13980, 3, 1053,
                                                                       1632, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14043, 3, 1068,
                                                                       1653, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14106, 3, 1083,
                                                                       1674, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 14169, 3, 1098,
                                                                       1695, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14232, 3, 1128,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14316, 3, 1149,
                                                                       1744, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14400, 3, 1170,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14484, 3, 1191,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14568, 3, 1212,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14652, 3, 1233,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14736, 3, 1254,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14820, 3, 1275,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14904, 3, 1296,
                                                                       1940, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14988, 3, 1317,
                                                                       1968, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15072, 3, 1338,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15156, 3, 1359,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15240, 3, 1380,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15324, 3, 1422,
                                                                       2080, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15408, 3, 1443,
                                                                       2108, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15492, 3, 1464,
                                                                       2136, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15576, 3, 1485,
                                                                       2164, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15660, 3, 1506,
                                                                       2192, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15744, 3, 1527,
                                                                       2220, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15828, 3, 1548,
                                                                       2248, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15912, 3, 1569,
                                                                       2276, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15996, 3, 1590,
                                                                       2304, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16080, 3, 1611,
                                                                       2332, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16164, 3, 1632,
                                                                       2360, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16248, 3, 1653,
                                                                       2388, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 16332, 3, 1674,
                                                                       2416, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16416, 3, 1716,
                                                                       2444, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16524, 3, 1744,
                                                                       2480, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16632, 3, 1772,
                                                                       2516, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16740, 3, 1800,
                                                                       2552, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16848, 3, 1828,
                                                                       2588, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16956, 3, 1856,
                                                                       2624, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17064, 3, 1884,
                                                                       2660, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17172, 3, 1912,
                                                                       2696, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17280, 3, 1940,
                                                                       2732, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17388, 3, 1968,
                                                                       2768, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17496, 3, 1996,
                                                                       2804, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17604, 3, 2024,
                                                                       2840, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17712, 3, 2080,
                                                                       2876, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17820, 3, 2108,
                                                                       2912, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17928, 3, 2136,
                                                                       2948, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18036, 3, 2164,
                                                                       2984, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18144, 3, 2192,
                                                                       3020, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18252, 3, 2220,
                                                                       3056, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18360, 3, 2248,
                                                                       3092, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18468, 3, 2276,
                                                                       3128, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18576, 3, 2304,
                                                                       3164, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18684, 3, 2332,
                                                                       3200, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18792, 3, 2360,
                                                                       3236, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 18900, 3, 2388,
                                                                       3272, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19008, 3, 2444,
                                                                       3308, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19143, 3, 2480,
                                                                       3353, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19278, 3, 2516,
                                                                       3398, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19413, 3, 2552,
                                                                       3443, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19548, 3, 2588,
                                                                       3488, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19683, 3, 2624,
                                                                       3533, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19818, 3, 2660,
                                                                       3578, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19953, 3, 2696,
                                                                       3623, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20088, 3, 2732,
                                                                       3668, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20223, 3, 2768,
                                                                       3713, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20358, 3, 2804,
                                                                       3758, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20493, 3, 2876,
                                                                       3803, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20628, 3, 2912,
                                                                       3848, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20763, 3, 2948,
                                                                       3893, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 20898, 3, 2984,
                                                                       3938, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21033, 3, 3020,
                                                                       3983, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21168, 3, 3056,
                                                                       4028, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21303, 3, 3092,
                                                                       4073, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21438, 3, 3128,
                                                                       4118, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21573, 3, 3164,
                                                                       4163, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21708, 3, 3200,
                                                                       4208, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 21843, 3, 3236,
                                                                       4253, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21978, 3, 3308,
                                                                       4298, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22143, 3, 3353,
                                                                       4353, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22308, 3, 3398,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22473, 3, 3443,
                                                                       4463, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22638, 3, 3488,
                                                                       4518, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22803, 3, 3533,
                                                                       4573, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22968, 3, 3578,
                                                                       4628, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 23133, 3, 3623,
                                                                       4683, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 23298, 3, 3668,
                                                                       4738, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 23463, 3, 3713,
                                                                       4793, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 23628, 3, 3803,
                                                                       4848, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 23793, 3, 3848,
                                                                       4903, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 23958, 3, 3893,
                                                                       4958, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24123, 3, 3938,
                                                                       5013, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24288, 3, 3983,
                                                                       5068, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24453, 3, 4028,
                                                                       5123, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24618, 3, 4073,
                                                                       5178, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24783, 3, 4118,
                                                                       5233, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 24948, 3, 4163,
                                                                       5288, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 25113, 3, 4208,
                                                                       5343, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 25278, 3, 4298,
                                                                       5398, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 25476, 3, 4353,
                                                                       5464, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 25674, 3, 4408,
                                                                       5530, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 25872, 3, 4463,
                                                                       5596, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 26070, 3, 4518,
                                                                       5662, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 26268, 3, 4573,
                                                                       5728, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 26466, 3, 4628,
                                                                       5794, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 26664, 3, 4683,
                                                                       5860, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 26862, 3, 4738,
                                                                       5926, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27060, 3, 4848,
                                                                       5992, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27258, 3, 4903,
                                                                       6058, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27456, 3, 4958,
                                                                       6124, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27654, 3, 5013,
                                                                       6190, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 27852, 3, 5068,
                                                                       6256, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28050, 3, 5123,
                                                                       6322, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28248, 3, 5178,
                                                                       6388, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28446, 3, 5233,
                                                                       6454, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 28644, 3, 5288,
                                                                       6520, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 28842, 3, 5398,
                                                                       6586, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 29076, 3, 5464,
                                                                       6664, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 29310, 3, 5530,
                                                                       6742, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 29544, 3, 5596,
                                                                       6820, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 29778, 3, 5662,
                                                                       6898, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 30012, 3, 5728,
                                                                       6976, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 30246, 3, 5794,
                                                                       7054, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 30480, 3, 5860,
                                                                       7132, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 30714, 3, 5992,
                                                                       7210, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 30948, 3, 6058,
                                                                       7288, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 31182, 3, 6124,
                                                                       7366, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 31416, 3, 6190,
                                                                       7444, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 31650, 3, 6256,
                                                                       7522, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 31884, 3, 6322,
                                                                       7600, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 32118, 3, 6388,
                                                                       7678, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 32352, 3, 6454,
                                                                       7756, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 32586, 3, 6586,
                                                                       7834, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 32859, 3, 6664,
                                                                       7925, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 33132, 3, 6742,
                                                                       8016, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 33405, 3, 6820,
                                                                       8107, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 33678, 3, 6898,
                                                                       8198, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 33951, 3, 6976,
                                                                       8289, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 34224, 3, 7054,
                                                                       8380, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 34497, 3, 7210,
                                                                       8471, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 34770, 3, 7288,
                                                                       8562, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 35043, 3, 7366,
                                                                       8653, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 35316, 3, 7444,
                                                                       8744, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 35589, 3, 7522,
                                                                       8835, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 35862, 3, 7600,
                                                                       8926, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 36135, 3, 7678,
                                                                       9017, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36408, 3, 7, 8,
                                                                       9114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36414, 3, 8, 9,
                                                                       9117, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36420, 3, 9, 10,
                                                                       9120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36426, 3, 10, 11,
                                                                       9123, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36432, 3, 11, 12,
                                                                       9126, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36438, 3, 12, 13,
                                                                       9129, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36444, 3, 13, 14,
                                                                       9132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36450, 3, 14, 15,
                                                                       9135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36456, 3, 15, 16,
                                                                       9138, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36462, 3, 16, 17,
                                                                       9141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36468, 3, 17, 18,
                                                                       9144, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36474, 3, 18, 19,
                                                                       9147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36480, 3, 19, 20,
                                                                       9150, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36486, 3, 20, 21,
                                                                       9153, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36492, 3, 21, 22,
                                                                       9156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36498, 3, 22, 23,
                                                                       9159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36504, 3, 23, 24,
                                                                       9162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36510, 3, 27, 28,
                                                                       9171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36516, 3, 28, 29,
                                                                       9174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36522, 3, 29, 30,
                                                                       9177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36528, 3, 30, 31,
                                                                       9180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36534, 3, 31, 32,
                                                                       9183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36540, 3, 32, 33,
                                                                       9186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36546, 3, 33, 34,
                                                                       9189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36552, 3, 34, 35,
                                                                       9192, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36558, 3, 35, 36,
                                                                       9195, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36564, 3, 36, 37,
                                                                       9198, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36570, 3, 37, 38,
                                                                       9201, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36576, 3, 38, 39,
                                                                       9204, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36582, 3, 39, 40,
                                                                       9207, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36588, 3, 40, 41,
                                                                       9210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36594, 3, 41, 42,
                                                                       9213, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36600, 3, 42, 43,
                                                                       9216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 36606, 3, 43, 44,
                                                                       9219, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36612, 0, 3,
                                                                       36408, 9114, 36414, 9240,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36630, 0, 3,
                                                                       36414, 9117, 36420, 9249,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36648, 0, 3,
                                                                       36420, 9120, 36426, 9258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36666, 0, 3,
                                                                       36426, 9123, 36432, 9267,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36684, 0, 3,
                                                                       36432, 9126, 36438, 9276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36702, 0, 3,
                                                                       36438, 9129, 36444, 9285,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36720, 0, 3,
                                                                       36444, 9132, 36450, 9294,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36738, 0, 3,
                                                                       36450, 9135, 36456, 9303,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36756, 0, 3,
                                                                       36456, 9138, 36462, 9312,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36774, 0, 3,
                                                                       36462, 9141, 36468, 9321,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36792, 0, 3,
                                                                       36468, 9144, 36474, 9330,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36810, 0, 3,
                                                                       36474, 9147, 36480, 9339,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36828, 0, 3,
                                                                       36480, 9150, 36486, 9348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36846, 0, 3,
                                                                       36486, 9153, 36492, 9357,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36864, 0, 3,
                                                                       36492, 9156, 36498, 9366,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36882, 0, 3,
                                                                       36498, 9159, 36504, 9375,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36900, 0, 3,
                                                                       36510, 9171, 36516, 9402,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36918, 0, 3,
                                                                       36516, 9174, 36522, 9411,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36936, 0, 3,
                                                                       36522, 9177, 36528, 9420,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36954, 0, 3,
                                                                       36528, 9180, 36534, 9429,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36972, 0, 3,
                                                                       36534, 9183, 36540, 9438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 36990, 0, 3,
                                                                       36540, 9186, 36546, 9447,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37008, 0, 3,
                                                                       36546, 9189, 36552, 9456,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37026, 0, 3,
                                                                       36552, 9192, 36558, 9465,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37044, 0, 3,
                                                                       36558, 9195, 36564, 9474,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37062, 0, 3,
                                                                       36564, 9198, 36570, 9483,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37080, 0, 3,
                                                                       36570, 9201, 36576, 9492,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37098, 0, 3,
                                                                       36576, 9204, 36582, 9501,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37116, 0, 3,
                                                                       36582, 9207, 36588, 9510,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37134, 0, 3,
                                                                       36588, 9210, 36594, 9519,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37152, 0, 3,
                                                                       36594, 9213, 36600, 9528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 37170, 0, 3,
                                                                       36600, 9216, 36606, 9537,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37188, 0, 3,
                                                                       36612, 9240, 36630, 154,
                                                                       160, 9582, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37224, 0, 3,
                                                                       36630, 9249, 36648, 160,
                                                                       166, 9600, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37260, 0, 3,
                                                                       36648, 9258, 36666, 166,
                                                                       172, 9618, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37296, 0, 3,
                                                                       36666, 9267, 36684, 172,
                                                                       178, 9636, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37332, 0, 3,
                                                                       36684, 9276, 36702, 178,
                                                                       184, 9654, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37368, 0, 3,
                                                                       36702, 9285, 36720, 184,
                                                                       190, 9672, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37404, 0, 3,
                                                                       36720, 9294, 36738, 190,
                                                                       196, 9690, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37440, 0, 3,
                                                                       36738, 9303, 36756, 196,
                                                                       202, 9708, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37476, 0, 3,
                                                                       36756, 9312, 36774, 202,
                                                                       208, 9726, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37512, 0, 3,
                                                                       36774, 9321, 36792, 208,
                                                                       214, 9744, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37548, 0, 3,
                                                                       36792, 9330, 36810, 214,
                                                                       220, 9762, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37584, 0, 3,
                                                                       36810, 9339, 36828, 220,
                                                                       226, 9780, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37620, 0, 3,
                                                                       36828, 9348, 36846, 226,
                                                                       232, 9798, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37656, 0, 3,
                                                                       36846, 9357, 36864, 232,
                                                                       238, 9816, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37692, 0, 3,
                                                                       36864, 9366, 36882, 238,
                                                                       244, 9834, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37728, 0, 3,
                                                                       36900, 9402, 36918, 256,
                                                                       262, 9888, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37764, 0, 3,
                                                                       36918, 9411, 36936, 262,
                                                                       268, 9906, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37800, 0, 3,
                                                                       36936, 9420, 36954, 268,
                                                                       274, 9924, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37836, 0, 3,
                                                                       36954, 9429, 36972, 274,
                                                                       280, 9942, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37872, 0, 3,
                                                                       36972, 9438, 36990, 280,
                                                                       286, 9960, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37908, 0, 3,
                                                                       36990, 9447, 37008, 286,
                                                                       292, 9978, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37944, 0, 3,
                                                                       37008, 9456, 37026, 292,
                                                                       298, 9996, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 37980, 0, 3,
                                                                       37026, 9465, 37044, 298,
                                                                       304, 10014, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 38016, 0, 3,
                                                                       37044, 9474, 37062, 304,
                                                                       310, 10032, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 38052, 0, 3,
                                                                       37062, 9483, 37080, 310,
                                                                       316, 10050, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 38088, 0, 3,
                                                                       37080, 9492, 37098, 316,
                                                                       322, 10068, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 38124, 0, 3,
                                                                       37098, 9501, 37116, 322,
                                                                       328, 10086, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 38160, 0, 3,
                                                                       37116, 9510, 37134, 328,
                                                                       334, 10104, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 38196, 0, 3,
                                                                       37134, 9519, 37152, 334,
                                                                       340, 10122, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 38232, 0, 3,
                                                                       37152, 9528, 37170, 340,
                                                                       346, 10140, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38268, 0, 3,
                                                                       37188, 9582, 37224, 358,
                                                                       368, 10218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38328, 0, 3,
                                                                       37224, 9600, 37260, 368,
                                                                       378, 10248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38388, 0, 3,
                                                                       37260, 9618, 37296, 378,
                                                                       388, 10278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38448, 0, 3,
                                                                       37296, 9636, 37332, 388,
                                                                       398, 10308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38508, 0, 3,
                                                                       37332, 9654, 37368, 398,
                                                                       408, 10338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38568, 0, 3,
                                                                       37368, 9672, 37404, 408,
                                                                       418, 10368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38628, 0, 3,
                                                                       37404, 9690, 37440, 418,
                                                                       428, 10398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38688, 0, 3,
                                                                       37440, 9708, 37476, 428,
                                                                       438, 10428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38748, 0, 3,
                                                                       37476, 9726, 37512, 438,
                                                                       448, 10458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38808, 0, 3,
                                                                       37512, 9744, 37548, 448,
                                                                       458, 10488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38868, 0, 3,
                                                                       37548, 9762, 37584, 458,
                                                                       468, 10518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38928, 0, 3,
                                                                       37584, 9780, 37620, 468,
                                                                       478, 10548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 38988, 0, 3,
                                                                       37620, 9798, 37656, 478,
                                                                       488, 10578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39048, 0, 3,
                                                                       37656, 9816, 37692, 488,
                                                                       498, 10608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39108, 0, 3,
                                                                       37728, 9888, 37764, 518,
                                                                       528, 10698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39168, 0, 3,
                                                                       37764, 9906, 37800, 528,
                                                                       538, 10728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39228, 0, 3,
                                                                       37800, 9924, 37836, 538,
                                                                       548, 10758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39288, 0, 3,
                                                                       37836, 9942, 37872, 548,
                                                                       558, 10788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39348, 0, 3,
                                                                       37872, 9960, 37908, 558,
                                                                       568, 10818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39408, 0, 3,
                                                                       37908, 9978, 37944, 568,
                                                                       578, 10848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39468, 0, 3,
                                                                       37944, 9996, 37980, 578,
                                                                       588, 10878, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39528, 0, 3,
                                                                       37980, 10014, 38016, 588,
                                                                       598, 10908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39588, 0, 3,
                                                                       38016, 10032, 38052, 598,
                                                                       608, 10938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39648, 0, 3,
                                                                       38052, 10050, 38088, 608,
                                                                       618, 10968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39708, 0, 3,
                                                                       38088, 10068, 38124, 618,
                                                                       628, 10998, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39768, 0, 3,
                                                                       38124, 10086, 38160, 628,
                                                                       638, 11028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39828, 0, 3,
                                                                       38160, 10104, 38196, 638,
                                                                       648, 11058, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 39888, 0, 3,
                                                                       38196, 10122, 38232, 648,
                                                                       658, 11088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 39948, 0, 3,
                                                                       38268, 10218, 38328, 678,
                                                                       693, 11208, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40038, 0, 3,
                                                                       38328, 10248, 38388, 693,
                                                                       708, 11253, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40128, 0, 3,
                                                                       38388, 10278, 38448, 708,
                                                                       723, 11298, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40218, 0, 3,
                                                                       38448, 10308, 38508, 723,
                                                                       738, 11343, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40308, 0, 3,
                                                                       38508, 10338, 38568, 738,
                                                                       753, 11388, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40398, 0, 3,
                                                                       38568, 10368, 38628, 753,
                                                                       768, 11433, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40488, 0, 3,
                                                                       38628, 10398, 38688, 768,
                                                                       783, 11478, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40578, 0, 3,
                                                                       38688, 10428, 38748, 783,
                                                                       798, 11523, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40668, 0, 3,
                                                                       38748, 10458, 38808, 798,
                                                                       813, 11568, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40758, 0, 3,
                                                                       38808, 10488, 38868, 813,
                                                                       828, 11613, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40848, 0, 3,
                                                                       38868, 10518, 38928, 828,
                                                                       843, 11658, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 40938, 0, 3,
                                                                       38928, 10548, 38988, 843,
                                                                       858, 11703, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41028, 0, 3,
                                                                       38988, 10578, 39048, 858,
                                                                       873, 11748, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41118, 0, 3,
                                                                       39108, 10698, 39168, 903,
                                                                       918, 11883, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41208, 0, 3,
                                                                       39168, 10728, 39228, 918,
                                                                       933, 11928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41298, 0, 3,
                                                                       39228, 10758, 39288, 933,
                                                                       948, 11973, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41388, 0, 3,
                                                                       39288, 10788, 39348, 948,
                                                                       963, 12018, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41478, 0, 3,
                                                                       39348, 10818, 39408, 963,
                                                                       978, 12063, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41568, 0, 3,
                                                                       39408, 10848, 39468, 978,
                                                                       993, 12108, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41658, 0, 3,
                                                                       39468, 10878, 39528, 993,
                                                                       1008, 12153, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41748, 0, 3,
                                                                       39528, 10908, 39588, 1008,
                                                                       1023, 12198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41838, 0, 3,
                                                                       39588, 10938, 39648, 1023,
                                                                       1038, 12243, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 41928, 0, 3,
                                                                       39648, 10968, 39708, 1038,
                                                                       1053, 12288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42018, 0, 3,
                                                                       39708, 10998, 39768, 1053,
                                                                       1068, 12333, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42108, 0, 3,
                                                                       39768, 11028, 39828, 1068,
                                                                       1083, 12378, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 42198, 0, 3,
                                                                       39828, 11058, 39888, 1083,
                                                                       1098, 12423, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 42288, 0, 3,
                                                                       39948, 11208, 40038, 1128,
                                                                       1149, 12594, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 42414, 0, 3,
                                                                       40038, 11253, 40128, 1149,
                                                                       1170, 12657, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 42540, 0, 3,
                                                                       40128, 11298, 40218, 1170,
                                                                       1191, 12720, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 42666, 0, 3,
                                                                       40218, 11343, 40308, 1191,
                                                                       1212, 12783, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 42792, 0, 3,
                                                                       40308, 11388, 40398, 1212,
                                                                       1233, 12846, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 42918, 0, 3,
                                                                       40398, 11433, 40488, 1233,
                                                                       1254, 12909, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 43044, 0, 3,
                                                                       40488, 11478, 40578, 1254,
                                                                       1275, 12972, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 43170, 0, 3,
                                                                       40578, 11523, 40668, 1275,
                                                                       1296, 13035, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 43296, 0, 3,
                                                                       40668, 11568, 40758, 1296,
                                                                       1317, 13098, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 43422, 0, 3,
                                                                       40758, 11613, 40848, 1317,
                                                                       1338, 13161, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 43548, 0, 3,
                                                                       40848, 11658, 40938, 1338,
                                                                       1359, 13224, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 43674, 0, 3,
                                                                       40938, 11703, 41028, 1359,
                                                                       1380, 13287, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 43800, 0, 3,
                                                                       41118, 11883, 41208, 1422,
                                                                       1443, 13476, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 43926, 0, 3,
                                                                       41208, 11928, 41298, 1443,
                                                                       1464, 13539, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44052, 0, 3,
                                                                       41298, 11973, 41388, 1464,
                                                                       1485, 13602, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44178, 0, 3,
                                                                       41388, 12018, 41478, 1485,
                                                                       1506, 13665, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44304, 0, 3,
                                                                       41478, 12063, 41568, 1506,
                                                                       1527, 13728, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44430, 0, 3,
                                                                       41568, 12108, 41658, 1527,
                                                                       1548, 13791, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44556, 0, 3,
                                                                       41658, 12153, 41748, 1548,
                                                                       1569, 13854, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44682, 0, 3,
                                                                       41748, 12198, 41838, 1569,
                                                                       1590, 13917, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44808, 0, 3,
                                                                       41838, 12243, 41928, 1590,
                                                                       1611, 13980, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 44934, 0, 3,
                                                                       41928, 12288, 42018, 1611,
                                                                       1632, 14043, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45060, 0, 3,
                                                                       42018, 12333, 42108, 1632,
                                                                       1653, 14106, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 45186, 0, 3,
                                                                       42108, 12378, 42198, 1653,
                                                                       1674, 14169, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 45312, 0, 3,
                                                                       42288, 12594, 42414, 1716,
                                                                       1744, 14400, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 45480, 0, 3,
                                                                       42414, 12657, 42540, 1744,
                                                                       1772, 14484, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 45648, 0, 3,
                                                                       42540, 12720, 42666, 1772,
                                                                       1800, 14568, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 45816, 0, 3,
                                                                       42666, 12783, 42792, 1800,
                                                                       1828, 14652, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 45984, 0, 3,
                                                                       42792, 12846, 42918, 1828,
                                                                       1856, 14736, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 46152, 0, 3,
                                                                       42918, 12909, 43044, 1856,
                                                                       1884, 14820, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 46320, 0, 3,
                                                                       43044, 12972, 43170, 1884,
                                                                       1912, 14904, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 46488, 0, 3,
                                                                       43170, 13035, 43296, 1912,
                                                                       1940, 14988, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 46656, 0, 3,
                                                                       43296, 13098, 43422, 1940,
                                                                       1968, 15072, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 46824, 0, 3,
                                                                       43422, 13161, 43548, 1968,
                                                                       1996, 15156, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 46992, 0, 3,
                                                                       43548, 13224, 43674, 1996,
                                                                       2024, 15240, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 47160, 0, 3,
                                                                       43800, 13476, 43926, 2080,
                                                                       2108, 15492, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 47328, 0, 3,
                                                                       43926, 13539, 44052, 2108,
                                                                       2136, 15576, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 47496, 0, 3,
                                                                       44052, 13602, 44178, 2136,
                                                                       2164, 15660, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 47664, 0, 3,
                                                                       44178, 13665, 44304, 2164,
                                                                       2192, 15744, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 47832, 0, 3,
                                                                       44304, 13728, 44430, 2192,
                                                                       2220, 15828, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48000, 0, 3,
                                                                       44430, 13791, 44556, 2220,
                                                                       2248, 15912, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48168, 0, 3,
                                                                       44556, 13854, 44682, 2248,
                                                                       2276, 15996, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48336, 0, 3,
                                                                       44682, 13917, 44808, 2276,
                                                                       2304, 16080, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48504, 0, 3,
                                                                       44808, 13980, 44934, 2304,
                                                                       2332, 16164, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48672, 0, 3,
                                                                       44934, 14043, 45060, 2332,
                                                                       2360, 16248, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 48840, 0, 3,
                                                                       45060, 14106, 45186, 2360,
                                                                       2388, 16332, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 49008, 0, 3,
                                                                       45312, 14400, 45480, 2444,
                                                                       2480, 16632, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 49224, 0, 3,
                                                                       45480, 14484, 45648, 2480,
                                                                       2516, 16740, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 49440, 0, 3,
                                                                       45648, 14568, 45816, 2516,
                                                                       2552, 16848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 49656, 0, 3,
                                                                       45816, 14652, 45984, 2552,
                                                                       2588, 16956, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 49872, 0, 3,
                                                                       45984, 14736, 46152, 2588,
                                                                       2624, 17064, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 50088, 0, 3,
                                                                       46152, 14820, 46320, 2624,
                                                                       2660, 17172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 50304, 0, 3,
                                                                       46320, 14904, 46488, 2660,
                                                                       2696, 17280, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 50520, 0, 3,
                                                                       46488, 14988, 46656, 2696,
                                                                       2732, 17388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 50736, 0, 3,
                                                                       46656, 15072, 46824, 2732,
                                                                       2768, 17496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 50952, 0, 3,
                                                                       46824, 15156, 46992, 2768,
                                                                       2804, 17604, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 51168, 0, 3,
                                                                       47160, 15492, 47328, 2876,
                                                                       2912, 17928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 51384, 0, 3,
                                                                       47328, 15576, 47496, 2912,
                                                                       2948, 18036, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 51600, 0, 3,
                                                                       47496, 15660, 47664, 2948,
                                                                       2984, 18144, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 51816, 0, 3,
                                                                       47664, 15744, 47832, 2984,
                                                                       3020, 18252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 52032, 0, 3,
                                                                       47832, 15828, 48000, 3020,
                                                                       3056, 18360, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 52248, 0, 3,
                                                                       48000, 15912, 48168, 3056,
                                                                       3092, 18468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 52464, 0, 3,
                                                                       48168, 15996, 48336, 3092,
                                                                       3128, 18576, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 52680, 0, 3,
                                                                       48336, 16080, 48504, 3128,
                                                                       3164, 18684, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 52896, 0, 3,
                                                                       48504, 16164, 48672, 3164,
                                                                       3200, 18792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 53112, 0, 3,
                                                                       48672, 16248, 48840, 3200,
                                                                       3236, 18900, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 53328, 0, 3,
                                                                       49008, 16632, 49224, 3308,
                                                                       3353, 19278, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 53598, 0, 3,
                                                                       49224, 16740, 49440, 3353,
                                                                       3398, 19413, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 53868, 0, 3,
                                                                       49440, 16848, 49656, 3398,
                                                                       3443, 19548, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 54138, 0, 3,
                                                                       49656, 16956, 49872, 3443,
                                                                       3488, 19683, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 54408, 0, 3,
                                                                       49872, 17064, 50088, 3488,
                                                                       3533, 19818, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 54678, 0, 3,
                                                                       50088, 17172, 50304, 3533,
                                                                       3578, 19953, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 54948, 0, 3,
                                                                       50304, 17280, 50520, 3578,
                                                                       3623, 20088, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 55218, 0, 3,
                                                                       50520, 17388, 50736, 3623,
                                                                       3668, 20223, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 55488, 0, 3,
                                                                       50736, 17496, 50952, 3668,
                                                                       3713, 20358, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 55758, 0, 3,
                                                                       51168, 17928, 51384, 3803,
                                                                       3848, 20763, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 56028, 0, 3,
                                                                       51384, 18036, 51600, 3848,
                                                                       3893, 20898, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 56298, 0, 3,
                                                                       51600, 18144, 51816, 3893,
                                                                       3938, 21033, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 56568, 0, 3,
                                                                       51816, 18252, 52032, 3938,
                                                                       3983, 21168, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 56838, 0, 3,
                                                                       52032, 18360, 52248, 3983,
                                                                       4028, 21303, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 57108, 0, 3,
                                                                       52248, 18468, 52464, 4028,
                                                                       4073, 21438, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 57378, 0, 3,
                                                                       52464, 18576, 52680, 4073,
                                                                       4118, 21573, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 57648, 0, 3,
                                                                       52680, 18684, 52896, 4118,
                                                                       4163, 21708, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 57918, 0, 3,
                                                                       52896, 18792, 53112, 4163,
                                                                       4208, 21843, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 58188, 0, 3,
                                                                       53328, 19278, 53598, 4298,
                                                                       4353, 22308, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 58518, 0, 3,
                                                                       53598, 19413, 53868, 4353,
                                                                       4408, 22473, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 58848, 0, 3,
                                                                       53868, 19548, 54138, 4408,
                                                                       4463, 22638, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 59178, 0, 3,
                                                                       54138, 19683, 54408, 4463,
                                                                       4518, 22803, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 59508, 0, 3,
                                                                       54408, 19818, 54678, 4518,
                                                                       4573, 22968, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 59838, 0, 3,
                                                                       54678, 19953, 54948, 4573,
                                                                       4628, 23133, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 60168, 0, 3,
                                                                       54948, 20088, 55218, 4628,
                                                                       4683, 23298, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 60498, 0, 3,
                                                                       55218, 20223, 55488, 4683,
                                                                       4738, 23463, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 60828, 0, 3,
                                                                       55758, 20763, 56028, 4848,
                                                                       4903, 23958, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 61158, 0, 3,
                                                                       56028, 20898, 56298, 4903,
                                                                       4958, 24123, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 61488, 0, 3,
                                                                       56298, 21033, 56568, 4958,
                                                                       5013, 24288, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 61818, 0, 3,
                                                                       56568, 21168, 56838, 5013,
                                                                       5068, 24453, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 62148, 0, 3,
                                                                       56838, 21303, 57108, 5068,
                                                                       5123, 24618, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 62478, 0, 3,
                                                                       57108, 21438, 57378, 5123,
                                                                       5178, 24783, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 62808, 0, 3,
                                                                       57378, 21573, 57648, 5178,
                                                                       5233, 24948, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 63138, 0, 3,
                                                                       57648, 21708, 57918, 5233,
                                                                       5288, 25113, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 63468, 0, 3,
                                                                       58188, 22308, 58518, 5398,
                                                                       5464, 25674, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 63864, 0, 3,
                                                                       58518, 22473, 58848, 5464,
                                                                       5530, 25872, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 64260, 0, 3,
                                                                       58848, 22638, 59178, 5530,
                                                                       5596, 26070, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 64656, 0, 3,
                                                                       59178, 22803, 59508, 5596,
                                                                       5662, 26268, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 65052, 0, 3,
                                                                       59508, 22968, 59838, 5662,
                                                                       5728, 26466, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 65448, 0, 3,
                                                                       59838, 23133, 60168, 5728,
                                                                       5794, 26664, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 65844, 0, 3,
                                                                       60168, 23298, 60498, 5794,
                                                                       5860, 26862, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 66240, 0, 3,
                                                                       60828, 23958, 61158, 5992,
                                                                       6058, 27456, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 66636, 0, 3,
                                                                       61158, 24123, 61488, 6058,
                                                                       6124, 27654, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 67032, 0, 3,
                                                                       61488, 24288, 61818, 6124,
                                                                       6190, 27852, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 67428, 0, 3,
                                                                       61818, 24453, 62148, 6190,
                                                                       6256, 28050, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 67824, 0, 3,
                                                                       62148, 24618, 62478, 6256,
                                                                       6322, 28248, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 68220, 0, 3,
                                                                       62478, 24783, 62808, 6322,
                                                                       6388, 28446, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 68616, 0, 3,
                                                                       62808, 24948, 63138, 6388,
                                                                       6454, 28644, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 69012, 0, 3,
                                                                       63468, 25674, 63864, 6586,
                                                                       6664, 29310, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 69480, 0, 3,
                                                                       63864, 25872, 64260, 6664,
                                                                       6742, 29544, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 69948, 0, 3,
                                                                       64260, 26070, 64656, 6742,
                                                                       6820, 29778, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 70416, 0, 3,
                                                                       64656, 26268, 65052, 6820,
                                                                       6898, 30012, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 70884, 0, 3,
                                                                       65052, 26466, 65448, 6898,
                                                                       6976, 30246, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 71352, 0, 3,
                                                                       65448, 26664, 65844, 6976,
                                                                       7054, 30480, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 71820, 0, 3,
                                                                       66240, 27456, 66636, 7210,
                                                                       7288, 31182, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 72288, 0, 3,
                                                                       66636, 27654, 67032, 7288,
                                                                       7366, 31416, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 72756, 0, 3,
                                                                       67032, 27852, 67428, 7366,
                                                                       7444, 31650, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 73224, 0, 3,
                                                                       67428, 28050, 67824, 7444,
                                                                       7522, 31884, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 73692, 0, 3,
                                                                       67824, 28248, 68220, 7522,
                                                                       7600, 32118, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 74160, 0, 3,
                                                                       68220, 28446, 68616, 7600,
                                                                       7678, 32352, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 74628, 0, 3,
                                                                       69012, 29310, 69480, 7834,
                                                                       7925, 33132, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 75174, 0, 3,
                                                                       69480, 29544, 69948, 7925,
                                                                       8016, 33405, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 75720, 0, 3,
                                                                       69948, 29778, 70416, 8016,
                                                                       8107, 33678, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 76266, 0, 3,
                                                                       70416, 30012, 70884, 8107,
                                                                       8198, 33951, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 76812, 0, 3,
                                                                       70884, 30246, 71352, 8198,
                                                                       8289, 34224, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 77358, 0, 3,
                                                                       71820, 31182, 72288, 8471,
                                                                       8562, 35043, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 77904, 0, 3,
                                                                       72288, 31416, 72756, 8562,
                                                                       8653, 35316, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 78450, 0, 3,
                                                                       72756, 31650, 73224, 8653,
                                                                       8744, 35589, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 78996, 0, 3,
                                                                       73224, 31884, 73692, 8744,
                                                                       8835, 35862, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 79542, 0, 3,
                                                                       73692, 32118, 74160, 8835,
                                                                       8926, 36135, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80088, 3, 9108,
                                                                       9111, 36408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80098, 3, 9111,
                                                                       9114, 36414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80108, 3, 9114,
                                                                       9117, 36420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80118, 3, 9117,
                                                                       9120, 36426, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80128, 3, 9120,
                                                                       9123, 36432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80138, 3, 9123,
                                                                       9126, 36438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80148, 3, 9126,
                                                                       9129, 36444, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80158, 3, 9129,
                                                                       9132, 36450, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80168, 3, 9132,
                                                                       9135, 36456, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80178, 3, 9135,
                                                                       9138, 36462, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80188, 3, 9138,
                                                                       9141, 36468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80198, 3, 9141,
                                                                       9144, 36474, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80208, 3, 9144,
                                                                       9147, 36480, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80218, 3, 9147,
                                                                       9150, 36486, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80228, 3, 9150,
                                                                       9153, 36492, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80238, 3, 9153,
                                                                       9156, 36498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80248, 3, 9156,
                                                                       9159, 36504, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80258, 3, 9165,
                                                                       9168, 36510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80268, 3, 9168,
                                                                       9171, 36516, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80278, 3, 9171,
                                                                       9174, 36522, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80288, 3, 9174,
                                                                       9177, 36528, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80298, 3, 9177,
                                                                       9180, 36534, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80308, 3, 9180,
                                                                       9183, 36540, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80318, 3, 9183,
                                                                       9186, 36546, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80328, 3, 9186,
                                                                       9189, 36552, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80338, 3, 9189,
                                                                       9192, 36558, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80348, 3, 9192,
                                                                       9195, 36564, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80358, 3, 9195,
                                                                       9198, 36570, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80368, 3, 9198,
                                                                       9201, 36576, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80378, 3, 9201,
                                                                       9204, 36582, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80388, 3, 9204,
                                                                       9207, 36588, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80398, 3, 9207,
                                                                       9210, 36594, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80408, 3, 9210,
                                                                       9213, 36600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 80418, 3, 9213,
                                                                       9216, 36606, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80428, 0, 3,
                                                                       80088, 36408, 80098, 9222,
                                                                       9231, 36612, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80458, 0, 3,
                                                                       80098, 36414, 80108, 9231,
                                                                       9240, 36630, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80488, 0, 3,
                                                                       80108, 36420, 80118, 9240,
                                                                       9249, 36648, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80518, 0, 3,
                                                                       80118, 36426, 80128, 9249,
                                                                       9258, 36666, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80548, 0, 3,
                                                                       80128, 36432, 80138, 9258,
                                                                       9267, 36684, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80578, 0, 3,
                                                                       80138, 36438, 80148, 9267,
                                                                       9276, 36702, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80608, 0, 3,
                                                                       80148, 36444, 80158, 9276,
                                                                       9285, 36720, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80638, 0, 3,
                                                                       80158, 36450, 80168, 9285,
                                                                       9294, 36738, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80668, 0, 3,
                                                                       80168, 36456, 80178, 9294,
                                                                       9303, 36756, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80698, 0, 3,
                                                                       80178, 36462, 80188, 9303,
                                                                       9312, 36774, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80728, 0, 3,
                                                                       80188, 36468, 80198, 9312,
                                                                       9321, 36792, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80758, 0, 3,
                                                                       80198, 36474, 80208, 9321,
                                                                       9330, 36810, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80788, 0, 3,
                                                                       80208, 36480, 80218, 9330,
                                                                       9339, 36828, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80818, 0, 3,
                                                                       80218, 36486, 80228, 9339,
                                                                       9348, 36846, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80848, 0, 3,
                                                                       80228, 36492, 80238, 9348,
                                                                       9357, 36864, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80878, 0, 3,
                                                                       80238, 36498, 80248, 9357,
                                                                       9366, 36882, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80908, 0, 3,
                                                                       80258, 36510, 80268, 9384,
                                                                       9393, 36900, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80938, 0, 3,
                                                                       80268, 36516, 80278, 9393,
                                                                       9402, 36918, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80968, 0, 3,
                                                                       80278, 36522, 80288, 9402,
                                                                       9411, 36936, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 80998, 0, 3,
                                                                       80288, 36528, 80298, 9411,
                                                                       9420, 36954, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81028, 0, 3,
                                                                       80298, 36534, 80308, 9420,
                                                                       9429, 36972, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81058, 0, 3,
                                                                       80308, 36540, 80318, 9429,
                                                                       9438, 36990, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81088, 0, 3,
                                                                       80318, 36546, 80328, 9438,
                                                                       9447, 37008, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81118, 0, 3,
                                                                       80328, 36552, 80338, 9447,
                                                                       9456, 37026, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81148, 0, 3,
                                                                       80338, 36558, 80348, 9456,
                                                                       9465, 37044, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81178, 0, 3,
                                                                       80348, 36564, 80358, 9465,
                                                                       9474, 37062, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81208, 0, 3,
                                                                       80358, 36570, 80368, 9474,
                                                                       9483, 37080, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81238, 0, 3,
                                                                       80368, 36576, 80378, 9483,
                                                                       9492, 37098, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81268, 0, 3,
                                                                       80378, 36582, 80388, 9492,
                                                                       9501, 37116, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81298, 0, 3,
                                                                       80388, 36588, 80398, 9501,
                                                                       9510, 37134, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81328, 0, 3,
                                                                       80398, 36594, 80408, 9510,
                                                                       9519, 37152, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 81358, 0, 3,
                                                                       80408, 36600, 80418, 9519,
                                                                       9528, 37170, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81388, 0, 3,
                                                                       80428, 36612, 80458, 9546,
                                                                       9564, 37188, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81448, 0, 3,
                                                                       80458, 36630, 80488, 9564,
                                                                       9582, 37224, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81508, 0, 3,
                                                                       80488, 36648, 80518, 9582,
                                                                       9600, 37260, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81568, 0, 3,
                                                                       80518, 36666, 80548, 9600,
                                                                       9618, 37296, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81628, 0, 3,
                                                                       80548, 36684, 80578, 9618,
                                                                       9636, 37332, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81688, 0, 3,
                                                                       80578, 36702, 80608, 9636,
                                                                       9654, 37368, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81748, 0, 3,
                                                                       80608, 36720, 80638, 9654,
                                                                       9672, 37404, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81808, 0, 3,
                                                                       80638, 36738, 80668, 9672,
                                                                       9690, 37440, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81868, 0, 3,
                                                                       80668, 36756, 80698, 9690,
                                                                       9708, 37476, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81928, 0, 3,
                                                                       80698, 36774, 80728, 9708,
                                                                       9726, 37512, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 81988, 0, 3,
                                                                       80728, 36792, 80758, 9726,
                                                                       9744, 37548, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82048, 0, 3,
                                                                       80758, 36810, 80788, 9744,
                                                                       9762, 37584, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82108, 0, 3,
                                                                       80788, 36828, 80818, 9762,
                                                                       9780, 37620, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82168, 0, 3,
                                                                       80818, 36846, 80848, 9780,
                                                                       9798, 37656, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82228, 0, 3,
                                                                       80848, 36864, 80878, 9798,
                                                                       9816, 37692, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82288, 0, 3,
                                                                       80908, 36900, 80938, 9852,
                                                                       9870, 37728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82348, 0, 3,
                                                                       80938, 36918, 80968, 9870,
                                                                       9888, 37764, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82408, 0, 3,
                                                                       80968, 36936, 80998, 9888,
                                                                       9906, 37800, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82468, 0, 3,
                                                                       80998, 36954, 81028, 9906,
                                                                       9924, 37836, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82528, 0, 3,
                                                                       81028, 36972, 81058, 9924,
                                                                       9942, 37872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82588, 0, 3,
                                                                       81058, 36990, 81088, 9942,
                                                                       9960, 37908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82648, 0, 3,
                                                                       81088, 37008, 81118, 9960,
                                                                       9978, 37944, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82708, 0, 3,
                                                                       81118, 37026, 81148, 9978,
                                                                       9996, 37980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82768, 0, 3,
                                                                       81148, 37044, 81178, 9996,
                                                                       10014, 38016, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82828, 0, 3,
                                                                       81178, 37062, 81208,
                                                                       10014, 10032, 38052,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82888, 0, 3,
                                                                       81208, 37080, 81238,
                                                                       10032, 10050, 38088,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 82948, 0, 3,
                                                                       81238, 37098, 81268,
                                                                       10050, 10068, 38124,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 83008, 0, 3,
                                                                       81268, 37116, 81298,
                                                                       10068, 10086, 38160,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 83068, 0, 3,
                                                                       81298, 37134, 81328,
                                                                       10086, 10104, 38196,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 83128, 0, 3,
                                                                       81328, 37152, 81358,
                                                                       10104, 10122, 38232,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83188, 0, 3,
                                                                       81388, 37188, 81448,
                                                                       10158, 10188, 38268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83288, 0, 3,
                                                                       81448, 37224, 81508,
                                                                       10188, 10218, 38328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83388, 0, 3,
                                                                       81508, 37260, 81568,
                                                                       10218, 10248, 38388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83488, 0, 3,
                                                                       81568, 37296, 81628,
                                                                       10248, 10278, 38448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83588, 0, 3,
                                                                       81628, 37332, 81688,
                                                                       10278, 10308, 38508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83688, 0, 3,
                                                                       81688, 37368, 81748,
                                                                       10308, 10338, 38568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83788, 0, 3,
                                                                       81748, 37404, 81808,
                                                                       10338, 10368, 38628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83888, 0, 3,
                                                                       81808, 37440, 81868,
                                                                       10368, 10398, 38688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 83988, 0, 3,
                                                                       81868, 37476, 81928,
                                                                       10398, 10428, 38748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84088, 0, 3,
                                                                       81928, 37512, 81988,
                                                                       10428, 10458, 38808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84188, 0, 3,
                                                                       81988, 37548, 82048,
                                                                       10458, 10488, 38868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84288, 0, 3,
                                                                       82048, 37584, 82108,
                                                                       10488, 10518, 38928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84388, 0, 3,
                                                                       82108, 37620, 82168,
                                                                       10518, 10548, 38988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84488, 0, 3,
                                                                       82168, 37656, 82228,
                                                                       10548, 10578, 39048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84588, 0, 3,
                                                                       82288, 37728, 82348,
                                                                       10638, 10668, 39108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84688, 0, 3,
                                                                       82348, 37764, 82408,
                                                                       10668, 10698, 39168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84788, 0, 3,
                                                                       82408, 37800, 82468,
                                                                       10698, 10728, 39228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84888, 0, 3,
                                                                       82468, 37836, 82528,
                                                                       10728, 10758, 39288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 84988, 0, 3,
                                                                       82528, 37872, 82588,
                                                                       10758, 10788, 39348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85088, 0, 3,
                                                                       82588, 37908, 82648,
                                                                       10788, 10818, 39408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85188, 0, 3,
                                                                       82648, 37944, 82708,
                                                                       10818, 10848, 39468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85288, 0, 3,
                                                                       82708, 37980, 82768,
                                                                       10848, 10878, 39528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85388, 0, 3,
                                                                       82768, 38016, 82828,
                                                                       10878, 10908, 39588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85488, 0, 3,
                                                                       82828, 38052, 82888,
                                                                       10908, 10938, 39648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85588, 0, 3,
                                                                       82888, 38088, 82948,
                                                                       10938, 10968, 39708,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85688, 0, 3,
                                                                       82948, 38124, 83008,
                                                                       10968, 10998, 39768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85788, 0, 3,
                                                                       83008, 38160, 83068,
                                                                       10998, 11028, 39828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 85888, 0, 3,
                                                                       83068, 38196, 83128,
                                                                       11028, 11058, 39888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 85988, 0, 3,
                                                                       83188, 38268, 83288,
                                                                       11118, 11163, 39948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 86138, 0, 3,
                                                                       83288, 38328, 83388,
                                                                       11163, 11208, 40038,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 86288, 0, 3,
                                                                       83388, 38388, 83488,
                                                                       11208, 11253, 40128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 86438, 0, 3,
                                                                       83488, 38448, 83588,
                                                                       11253, 11298, 40218,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 86588, 0, 3,
                                                                       83588, 38508, 83688,
                                                                       11298, 11343, 40308,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 86738, 0, 3,
                                                                       83688, 38568, 83788,
                                                                       11343, 11388, 40398,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 86888, 0, 3,
                                                                       83788, 38628, 83888,
                                                                       11388, 11433, 40488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 87038, 0, 3,
                                                                       83888, 38688, 83988,
                                                                       11433, 11478, 40578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 87188, 0, 3,
                                                                       83988, 38748, 84088,
                                                                       11478, 11523, 40668,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 87338, 0, 3,
                                                                       84088, 38808, 84188,
                                                                       11523, 11568, 40758,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 87488, 0, 3,
                                                                       84188, 38868, 84288,
                                                                       11568, 11613, 40848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 87638, 0, 3,
                                                                       84288, 38928, 84388,
                                                                       11613, 11658, 40938,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 87788, 0, 3,
                                                                       84388, 38988, 84488,
                                                                       11658, 11703, 41028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 87938, 0, 3,
                                                                       84588, 39108, 84688,
                                                                       11793, 11838, 41118,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 88088, 0, 3,
                                                                       84688, 39168, 84788,
                                                                       11838, 11883, 41208,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 88238, 0, 3,
                                                                       84788, 39228, 84888,
                                                                       11883, 11928, 41298,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 88388, 0, 3,
                                                                       84888, 39288, 84988,
                                                                       11928, 11973, 41388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 88538, 0, 3,
                                                                       84988, 39348, 85088,
                                                                       11973, 12018, 41478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 88688, 0, 3,
                                                                       85088, 39408, 85188,
                                                                       12018, 12063, 41568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 88838, 0, 3,
                                                                       85188, 39468, 85288,
                                                                       12063, 12108, 41658,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 88988, 0, 3,
                                                                       85288, 39528, 85388,
                                                                       12108, 12153, 41748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 89138, 0, 3,
                                                                       85388, 39588, 85488,
                                                                       12153, 12198, 41838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 89288, 0, 3,
                                                                       85488, 39648, 85588,
                                                                       12198, 12243, 41928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 89438, 0, 3,
                                                                       85588, 39708, 85688,
                                                                       12243, 12288, 42018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 89588, 0, 3,
                                                                       85688, 39768, 85788,
                                                                       12288, 12333, 42108,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 89738, 0, 3,
                                                                       85788, 39828, 85888,
                                                                       12333, 12378, 42198,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 89888, 0, 3,
                                                                       85988, 39948, 86138,
                                                                       12468, 12531, 42288,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 90098, 0, 3,
                                                                       86138, 40038, 86288,
                                                                       12531, 12594, 42414,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 90308, 0, 3,
                                                                       86288, 40128, 86438,
                                                                       12594, 12657, 42540,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 90518, 0, 3,
                                                                       86438, 40218, 86588,
                                                                       12657, 12720, 42666,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 90728, 0, 3,
                                                                       86588, 40308, 86738,
                                                                       12720, 12783, 42792,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 90938, 0, 3,
                                                                       86738, 40398, 86888,
                                                                       12783, 12846, 42918,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 91148, 0, 3,
                                                                       86888, 40488, 87038,
                                                                       12846, 12909, 43044,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 91358, 0, 3,
                                                                       87038, 40578, 87188,
                                                                       12909, 12972, 43170,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 91568, 0, 3,
                                                                       87188, 40668, 87338,
                                                                       12972, 13035, 43296,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 91778, 0, 3,
                                                                       87338, 40758, 87488,
                                                                       13035, 13098, 43422,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 91988, 0, 3,
                                                                       87488, 40848, 87638,
                                                                       13098, 13161, 43548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 92198, 0, 3,
                                                                       87638, 40938, 87788,
                                                                       13161, 13224, 43674,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 92408, 0, 3,
                                                                       87938, 41118, 88088,
                                                                       13350, 13413, 43800,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 92618, 0, 3,
                                                                       88088, 41208, 88238,
                                                                       13413, 13476, 43926,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 92828, 0, 3,
                                                                       88238, 41298, 88388,
                                                                       13476, 13539, 44052,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 93038, 0, 3,
                                                                       88388, 41388, 88538,
                                                                       13539, 13602, 44178,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 93248, 0, 3,
                                                                       88538, 41478, 88688,
                                                                       13602, 13665, 44304,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 93458, 0, 3,
                                                                       88688, 41568, 88838,
                                                                       13665, 13728, 44430,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 93668, 0, 3,
                                                                       88838, 41658, 88988,
                                                                       13728, 13791, 44556,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 93878, 0, 3,
                                                                       88988, 41748, 89138,
                                                                       13791, 13854, 44682,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 94088, 0, 3,
                                                                       89138, 41838, 89288,
                                                                       13854, 13917, 44808,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 94298, 0, 3,
                                                                       89288, 41928, 89438,
                                                                       13917, 13980, 44934,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 94508, 0, 3,
                                                                       89438, 42018, 89588,
                                                                       13980, 14043, 45060,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 94718, 0, 3,
                                                                       89588, 42108, 89738,
                                                                       14043, 14106, 45186,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 94928, 0, 3,
                                                                       89888, 42288, 90098,
                                                                       14232, 14316, 45312,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 95208, 0, 3,
                                                                       90098, 42414, 90308,
                                                                       14316, 14400, 45480,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 95488, 0, 3,
                                                                       90308, 42540, 90518,
                                                                       14400, 14484, 45648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 95768, 0, 3,
                                                                       90518, 42666, 90728,
                                                                       14484, 14568, 45816,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 96048, 0, 3,
                                                                       90728, 42792, 90938,
                                                                       14568, 14652, 45984,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 96328, 0, 3,
                                                                       90938, 42918, 91148,
                                                                       14652, 14736, 46152,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 96608, 0, 3,
                                                                       91148, 43044, 91358,
                                                                       14736, 14820, 46320,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 96888, 0, 3,
                                                                       91358, 43170, 91568,
                                                                       14820, 14904, 46488,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 97168, 0, 3,
                                                                       91568, 43296, 91778,
                                                                       14904, 14988, 46656,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 97448, 0, 3,
                                                                       91778, 43422, 91988,
                                                                       14988, 15072, 46824,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 97728, 0, 3,
                                                                       91988, 43548, 92198,
                                                                       15072, 15156, 46992,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 98008, 0, 3,
                                                                       92408, 43800, 92618,
                                                                       15324, 15408, 47160,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 98288, 0, 3,
                                                                       92618, 43926, 92828,
                                                                       15408, 15492, 47328,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 98568, 0, 3,
                                                                       92828, 44052, 93038,
                                                                       15492, 15576, 47496,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 98848, 0, 3,
                                                                       93038, 44178, 93248,
                                                                       15576, 15660, 47664,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 99128, 0, 3,
                                                                       93248, 44304, 93458,
                                                                       15660, 15744, 47832,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 99408, 0, 3,
                                                                       93458, 44430, 93668,
                                                                       15744, 15828, 48000,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 99688, 0, 3,
                                                                       93668, 44556, 93878,
                                                                       15828, 15912, 48168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 99968, 0, 3,
                                                                       93878, 44682, 94088,
                                                                       15912, 15996, 48336,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 100248, 0, 3,
                                                                       94088, 44808, 94298,
                                                                       15996, 16080, 48504,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 100528, 0, 3,
                                                                       94298, 44934, 94508,
                                                                       16080, 16164, 48672,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 100808, 0, 3,
                                                                       94508, 45060, 94718,
                                                                       16164, 16248, 48840,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 101088, 0, 3,
                                                                       94928, 45312, 95208,
                                                                       16416, 16524, 49008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 101448, 0, 3,
                                                                       95208, 45480, 95488,
                                                                       16524, 16632, 49224,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 101808, 0, 3,
                                                                       95488, 45648, 95768,
                                                                       16632, 16740, 49440,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 102168, 0, 3,
                                                                       95768, 45816, 96048,
                                                                       16740, 16848, 49656,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 102528, 0, 3,
                                                                       96048, 45984, 96328,
                                                                       16848, 16956, 49872,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 102888, 0, 3,
                                                                       96328, 46152, 96608,
                                                                       16956, 17064, 50088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 103248, 0, 3,
                                                                       96608, 46320, 96888,
                                                                       17064, 17172, 50304,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 103608, 0, 3,
                                                                       96888, 46488, 97168,
                                                                       17172, 17280, 50520,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 103968, 0, 3,
                                                                       97168, 46656, 97448,
                                                                       17280, 17388, 50736,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 104328, 0, 3,
                                                                       97448, 46824, 97728,
                                                                       17388, 17496, 50952,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 104688, 0, 3,
                                                                       98008, 47160, 98288,
                                                                       17712, 17820, 51168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 105048, 0, 3,
                                                                       98288, 47328, 98568,
                                                                       17820, 17928, 51384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 105408, 0, 3,
                                                                       98568, 47496, 98848,
                                                                       17928, 18036, 51600,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 105768, 0, 3,
                                                                       98848, 47664, 99128,
                                                                       18036, 18144, 51816,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 106128, 0, 3,
                                                                       99128, 47832, 99408,
                                                                       18144, 18252, 52032,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 106488, 0, 3,
                                                                       99408, 48000, 99688,
                                                                       18252, 18360, 52248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 106848, 0, 3,
                                                                       99688, 48168, 99968,
                                                                       18360, 18468, 52464,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 107208, 0, 3,
                                                                       99968, 48336, 100248,
                                                                       18468, 18576, 52680,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 107568, 0, 3,
                                                                       100248, 48504, 100528,
                                                                       18576, 18684, 52896,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 107928, 0, 3,
                                                                       100528, 48672, 100808,
                                                                       18684, 18792, 53112,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 108288, 0, 3,
                                                                       101088, 49008, 101448,
                                                                       19008, 19143, 53328,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 108738, 0, 3,
                                                                       101448, 49224, 101808,
                                                                       19143, 19278, 53598,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 109188, 0, 3,
                                                                       101808, 49440, 102168,
                                                                       19278, 19413, 53868,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 109638, 0, 3,
                                                                       102168, 49656, 102528,
                                                                       19413, 19548, 54138,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 110088, 0, 3,
                                                                       102528, 49872, 102888,
                                                                       19548, 19683, 54408,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 110538, 0, 3,
                                                                       102888, 50088, 103248,
                                                                       19683, 19818, 54678,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 110988, 0, 3,
                                                                       103248, 50304, 103608,
                                                                       19818, 19953, 54948,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 111438, 0, 3,
                                                                       103608, 50520, 103968,
                                                                       19953, 20088, 55218,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 111888, 0, 3,
                                                                       103968, 50736, 104328,
                                                                       20088, 20223, 55488,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 112338, 0, 3,
                                                                       104688, 51168, 105048,
                                                                       20493, 20628, 55758,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 112788, 0, 3,
                                                                       105048, 51384, 105408,
                                                                       20628, 20763, 56028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 113238, 0, 3,
                                                                       105408, 51600, 105768,
                                                                       20763, 20898, 56298,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 113688, 0, 3,
                                                                       105768, 51816, 106128,
                                                                       20898, 21033, 56568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 114138, 0, 3,
                                                                       106128, 52032, 106488,
                                                                       21033, 21168, 56838,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 114588, 0, 3,
                                                                       106488, 52248, 106848,
                                                                       21168, 21303, 57108,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 115038, 0, 3,
                                                                       106848, 52464, 107208,
                                                                       21303, 21438, 57378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 115488, 0, 3,
                                                                       107208, 52680, 107568,
                                                                       21438, 21573, 57648,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 115938, 0, 3,
                                                                       107568, 52896, 107928,
                                                                       21573, 21708, 57918,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 116388, 0, 3,
                                                                       108288, 53328, 108738,
                                                                       21978, 22143, 58188,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 116938, 0, 3,
                                                                       108738, 53598, 109188,
                                                                       22143, 22308, 58518,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 117488, 0, 3,
                                                                       109188, 53868, 109638,
                                                                       22308, 22473, 58848,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 118038, 0, 3,
                                                                       109638, 54138, 110088,
                                                                       22473, 22638, 59178,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 118588, 0, 3,
                                                                       110088, 54408, 110538,
                                                                       22638, 22803, 59508,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 119138, 0, 3,
                                                                       110538, 54678, 110988,
                                                                       22803, 22968, 59838,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 119688, 0, 3,
                                                                       110988, 54948, 111438,
                                                                       22968, 23133, 60168,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 120238, 0, 3,
                                                                       111438, 55218, 111888,
                                                                       23133, 23298, 60498,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 120788, 0, 3,
                                                                       112338, 55758, 112788,
                                                                       23628, 23793, 60828,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 121338, 0, 3,
                                                                       112788, 56028, 113238,
                                                                       23793, 23958, 61158,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 121888, 0, 3,
                                                                       113238, 56298, 113688,
                                                                       23958, 24123, 61488,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 122438, 0, 3,
                                                                       113688, 56568, 114138,
                                                                       24123, 24288, 61818,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 122988, 0, 3,
                                                                       114138, 56838, 114588,
                                                                       24288, 24453, 62148,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 123538, 0, 3,
                                                                       114588, 57108, 115038,
                                                                       24453, 24618, 62478,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 124088, 0, 3,
                                                                       115038, 57378, 115488,
                                                                       24618, 24783, 62808,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 124638, 0, 3,
                                                                       115488, 57648, 115938,
                                                                       24783, 24948, 63138,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 125188, 0, 3,
                                                                       116388, 58188, 116938,
                                                                       25278, 25476, 63468,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 125848, 0, 3,
                                                                       116938, 58518, 117488,
                                                                       25476, 25674, 63864,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 126508, 0, 3,
                                                                       117488, 58848, 118038,
                                                                       25674, 25872, 64260,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 127168, 0, 3,
                                                                       118038, 59178, 118588,
                                                                       25872, 26070, 64656,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 127828, 0, 3,
                                                                       118588, 59508, 119138,
                                                                       26070, 26268, 65052,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 128488, 0, 3,
                                                                       119138, 59838, 119688,
                                                                       26268, 26466, 65448,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 129148, 0, 3,
                                                                       119688, 60168, 120238,
                                                                       26466, 26664, 65844,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 129808, 0, 3,
                                                                       120788, 60828, 121338,
                                                                       27060, 27258, 66240,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 130468, 0, 3,
                                                                       121338, 61158, 121888,
                                                                       27258, 27456, 66636,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 131128, 0, 3,
                                                                       121888, 61488, 122438,
                                                                       27456, 27654, 67032,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 131788, 0, 3,
                                                                       122438, 61818, 122988,
                                                                       27654, 27852, 67428,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 132448, 0, 3,
                                                                       122988, 62148, 123538,
                                                                       27852, 28050, 67824,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 133108, 0, 3,
                                                                       123538, 62478, 124088,
                                                                       28050, 28248, 68220,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 133768, 0, 3,
                                                                       124088, 62808, 124638,
                                                                       28248, 28446, 68616,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 134428, 0, 3,
                                                                       125188, 63468, 125848,
                                                                       28842, 29076, 69012,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 135208, 0, 3,
                                                                       125848, 63864, 126508,
                                                                       29076, 29310, 69480,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 135988, 0, 3,
                                                                       126508, 64260, 127168,
                                                                       29310, 29544, 69948,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 136768, 0, 3,
                                                                       127168, 64656, 127828,
                                                                       29544, 29778, 70416,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 137548, 0, 3,
                                                                       127828, 65052, 128488,
                                                                       29778, 30012, 70884,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 138328, 0, 3,
                                                                       128488, 65448, 129148,
                                                                       30012, 30246, 71352,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 139108, 0, 3,
                                                                       129808, 66240, 130468,
                                                                       30714, 30948, 71820,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 139888, 0, 3,
                                                                       130468, 66636, 131128,
                                                                       30948, 31182, 72288,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 140668, 0, 3,
                                                                       131128, 67032, 131788,
                                                                       31182, 31416, 72756,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 141448, 0, 3,
                                                                       131788, 67428, 132448,
                                                                       31416, 31650, 73224,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 142228, 0, 3,
                                                                       132448, 67824, 133108,
                                                                       31650, 31884, 73692,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 143008, 0, 3,
                                                                       133108, 68220, 133768,
                                                                       31884, 32118, 74160,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 143788, 0, 3,
                                                                       134428, 69012, 135208,
                                                                       32586, 32859, 74628,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 144698, 0, 3,
                                                                       135208, 69480, 135988,
                                                                       32859, 33132, 75174,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 145608, 0, 3,
                                                                       135988, 69948, 136768,
                                                                       33132, 33405, 75720,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 146518, 0, 3,
                                                                       136768, 70416, 137548,
                                                                       33405, 33678, 76266,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 147428, 0, 3,
                                                                       137548, 70884, 138328,
                                                                       33678, 33951, 76812,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 148338, 0, 3,
                                                                       139108, 71820, 139888,
                                                                       34497, 34770, 77358,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 149248, 0, 3,
                                                                       139888, 72288, 140668,
                                                                       34770, 35043, 77904,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 150158, 0, 3,
                                                                       140668, 72756, 141448,
                                                                       35043, 35316, 78450,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 151068, 0, 3,
                                                                       141448, 73224, 142228,
                                                                       35316, 35589, 78996,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 151978, 0, 3,
                                                                       142228, 73692, 143008,
                                                                       35589, 35862, 79542,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 152888, 3, 36408,
                                                                       36414, 80108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 152903, 3, 36414,
                                                                       36420, 80118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 152918, 3, 36420,
                                                                       36426, 80128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 152933, 3, 36426,
                                                                       36432, 80138, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 152948, 3, 36432,
                                                                       36438, 80148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 152963, 3, 36438,
                                                                       36444, 80158, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 152978, 3, 36444,
                                                                       36450, 80168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 152993, 3, 36450,
                                                                       36456, 80178, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153008, 3, 36456,
                                                                       36462, 80188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153023, 3, 36462,
                                                                       36468, 80198, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153038, 3, 36468,
                                                                       36474, 80208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153053, 3, 36474,
                                                                       36480, 80218, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153068, 3, 36480,
                                                                       36486, 80228, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153083, 3, 36486,
                                                                       36492, 80238, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153098, 3, 36492,
                                                                       36498, 80248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153113, 3, 36510,
                                                                       36516, 80278, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153128, 3, 36516,
                                                                       36522, 80288, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153143, 3, 36522,
                                                                       36528, 80298, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153158, 3, 36528,
                                                                       36534, 80308, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153173, 3, 36534,
                                                                       36540, 80318, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153188, 3, 36540,
                                                                       36546, 80328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153203, 3, 36546,
                                                                       36552, 80338, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153218, 3, 36552,
                                                                       36558, 80348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153233, 3, 36558,
                                                                       36564, 80358, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153248, 3, 36564,
                                                                       36570, 80368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153263, 3, 36570,
                                                                       36576, 80378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153278, 3, 36576,
                                                                       36582, 80388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153293, 3, 36582,
                                                                       36588, 80398, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153308, 3, 36588,
                                                                       36594, 80408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 153323, 3, 36594,
                                                                       36600, 80418, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153338, 0, 3,
                                                                       152888, 80108, 152903,
                                                                       36612, 36630, 80488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153383, 0, 3,
                                                                       152903, 80118, 152918,
                                                                       36630, 36648, 80518,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153428, 0, 3,
                                                                       152918, 80128, 152933,
                                                                       36648, 36666, 80548,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153473, 0, 3,
                                                                       152933, 80138, 152948,
                                                                       36666, 36684, 80578,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153518, 0, 3,
                                                                       152948, 80148, 152963,
                                                                       36684, 36702, 80608,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153563, 0, 3,
                                                                       152963, 80158, 152978,
                                                                       36702, 36720, 80638,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153608, 0, 3,
                                                                       152978, 80168, 152993,
                                                                       36720, 36738, 80668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153653, 0, 3,
                                                                       152993, 80178, 153008,
                                                                       36738, 36756, 80698,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153698, 0, 3,
                                                                       153008, 80188, 153023,
                                                                       36756, 36774, 80728,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153743, 0, 3,
                                                                       153023, 80198, 153038,
                                                                       36774, 36792, 80758,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153788, 0, 3,
                                                                       153038, 80208, 153053,
                                                                       36792, 36810, 80788,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153833, 0, 3,
                                                                       153053, 80218, 153068,
                                                                       36810, 36828, 80818,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153878, 0, 3,
                                                                       153068, 80228, 153083,
                                                                       36828, 36846, 80848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153923, 0, 3,
                                                                       153083, 80238, 153098,
                                                                       36846, 36864, 80878,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 153968, 0, 3,
                                                                       153113, 80278, 153128,
                                                                       36900, 36918, 80968,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154013, 0, 3,
                                                                       153128, 80288, 153143,
                                                                       36918, 36936, 80998,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154058, 0, 3,
                                                                       153143, 80298, 153158,
                                                                       36936, 36954, 81028,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154103, 0, 3,
                                                                       153158, 80308, 153173,
                                                                       36954, 36972, 81058,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154148, 0, 3,
                                                                       153173, 80318, 153188,
                                                                       36972, 36990, 81088,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154193, 0, 3,
                                                                       153188, 80328, 153203,
                                                                       36990, 37008, 81118,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154238, 0, 3,
                                                                       153203, 80338, 153218,
                                                                       37008, 37026, 81148,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154283, 0, 3,
                                                                       153218, 80348, 153233,
                                                                       37026, 37044, 81178,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154328, 0, 3,
                                                                       153233, 80358, 153248,
                                                                       37044, 37062, 81208,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154373, 0, 3,
                                                                       153248, 80368, 153263,
                                                                       37062, 37080, 81238,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154418, 0, 3,
                                                                       153263, 80378, 153278,
                                                                       37080, 37098, 81268,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154463, 0, 3,
                                                                       153278, 80388, 153293,
                                                                       37098, 37116, 81298,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154508, 0, 3,
                                                                       153293, 80398, 153308,
                                                                       37116, 37134, 81328,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 154553, 0, 3,
                                                                       153308, 80408, 153323,
                                                                       37134, 37152, 81358,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 154598, 0, 3,
                                                                       153338, 80488, 153383,
                                                                       37188, 37224, 81508,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 154688, 0, 3,
                                                                       153383, 80518, 153428,
                                                                       37224, 37260, 81568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 154778, 0, 3,
                                                                       153428, 80548, 153473,
                                                                       37260, 37296, 81628,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 154868, 0, 3,
                                                                       153473, 80578, 153518,
                                                                       37296, 37332, 81688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 154958, 0, 3,
                                                                       153518, 80608, 153563,
                                                                       37332, 37368, 81748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155048, 0, 3,
                                                                       153563, 80638, 153608,
                                                                       37368, 37404, 81808,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155138, 0, 3,
                                                                       153608, 80668, 153653,
                                                                       37404, 37440, 81868,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155228, 0, 3,
                                                                       153653, 80698, 153698,
                                                                       37440, 37476, 81928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155318, 0, 3,
                                                                       153698, 80728, 153743,
                                                                       37476, 37512, 81988,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155408, 0, 3,
                                                                       153743, 80758, 153788,
                                                                       37512, 37548, 82048,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155498, 0, 3,
                                                                       153788, 80788, 153833,
                                                                       37548, 37584, 82108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155588, 0, 3,
                                                                       153833, 80818, 153878,
                                                                       37584, 37620, 82168,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155678, 0, 3,
                                                                       153878, 80848, 153923,
                                                                       37620, 37656, 82228,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155768, 0, 3,
                                                                       153968, 80968, 154013,
                                                                       37728, 37764, 82408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155858, 0, 3,
                                                                       154013, 80998, 154058,
                                                                       37764, 37800, 82468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 155948, 0, 3,
                                                                       154058, 81028, 154103,
                                                                       37800, 37836, 82528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156038, 0, 3,
                                                                       154103, 81058, 154148,
                                                                       37836, 37872, 82588,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156128, 0, 3,
                                                                       154148, 81088, 154193,
                                                                       37872, 37908, 82648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156218, 0, 3,
                                                                       154193, 81118, 154238,
                                                                       37908, 37944, 82708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156308, 0, 3,
                                                                       154238, 81148, 154283,
                                                                       37944, 37980, 82768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156398, 0, 3,
                                                                       154283, 81178, 154328,
                                                                       37980, 38016, 82828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156488, 0, 3,
                                                                       154328, 81208, 154373,
                                                                       38016, 38052, 82888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156578, 0, 3,
                                                                       154373, 81238, 154418,
                                                                       38052, 38088, 82948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156668, 0, 3,
                                                                       154418, 81268, 154463,
                                                                       38088, 38124, 83008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156758, 0, 3,
                                                                       154463, 81298, 154508,
                                                                       38124, 38160, 83068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 156848, 0, 3,
                                                                       154508, 81328, 154553,
                                                                       38160, 38196, 83128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 156938, 0, 3,
                                                                       154598, 81508, 154688,
                                                                       38268, 38328, 83388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 157088, 0, 3,
                                                                       154688, 81568, 154778,
                                                                       38328, 38388, 83488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 157238, 0, 3,
                                                                       154778, 81628, 154868,
                                                                       38388, 38448, 83588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 157388, 0, 3,
                                                                       154868, 81688, 154958,
                                                                       38448, 38508, 83688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 157538, 0, 3,
                                                                       154958, 81748, 155048,
                                                                       38508, 38568, 83788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 157688, 0, 3,
                                                                       155048, 81808, 155138,
                                                                       38568, 38628, 83888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 157838, 0, 3,
                                                                       155138, 81868, 155228,
                                                                       38628, 38688, 83988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 157988, 0, 3,
                                                                       155228, 81928, 155318,
                                                                       38688, 38748, 84088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 158138, 0, 3,
                                                                       155318, 81988, 155408,
                                                                       38748, 38808, 84188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 158288, 0, 3,
                                                                       155408, 82048, 155498,
                                                                       38808, 38868, 84288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 158438, 0, 3,
                                                                       155498, 82108, 155588,
                                                                       38868, 38928, 84388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 158588, 0, 3,
                                                                       155588, 82168, 155678,
                                                                       38928, 38988, 84488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 158738, 0, 3,
                                                                       155768, 82408, 155858,
                                                                       39108, 39168, 84788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 158888, 0, 3,
                                                                       155858, 82468, 155948,
                                                                       39168, 39228, 84888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 159038, 0, 3,
                                                                       155948, 82528, 156038,
                                                                       39228, 39288, 84988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 159188, 0, 3,
                                                                       156038, 82588, 156128,
                                                                       39288, 39348, 85088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 159338, 0, 3,
                                                                       156128, 82648, 156218,
                                                                       39348, 39408, 85188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 159488, 0, 3,
                                                                       156218, 82708, 156308,
                                                                       39408, 39468, 85288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 159638, 0, 3,
                                                                       156308, 82768, 156398,
                                                                       39468, 39528, 85388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 159788, 0, 3,
                                                                       156398, 82828, 156488,
                                                                       39528, 39588, 85488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 159938, 0, 3,
                                                                       156488, 82888, 156578,
                                                                       39588, 39648, 85588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 160088, 0, 3,
                                                                       156578, 82948, 156668,
                                                                       39648, 39708, 85688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 160238, 0, 3,
                                                                       156668, 83008, 156758,
                                                                       39708, 39768, 85788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 160388, 0, 3,
                                                                       156758, 83068, 156848,
                                                                       39768, 39828, 85888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 160538, 0, 3,
                                                                       156938, 83388, 157088,
                                                                       39948, 40038, 86288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 160763, 0, 3,
                                                                       157088, 83488, 157238,
                                                                       40038, 40128, 86438,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 160988, 0, 3,
                                                                       157238, 83588, 157388,
                                                                       40128, 40218, 86588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 161213, 0, 3,
                                                                       157388, 83688, 157538,
                                                                       40218, 40308, 86738,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 161438, 0, 3,
                                                                       157538, 83788, 157688,
                                                                       40308, 40398, 86888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 161663, 0, 3,
                                                                       157688, 83888, 157838,
                                                                       40398, 40488, 87038,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 161888, 0, 3,
                                                                       157838, 83988, 157988,
                                                                       40488, 40578, 87188,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 162113, 0, 3,
                                                                       157988, 84088, 158138,
                                                                       40578, 40668, 87338,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 162338, 0, 3,
                                                                       158138, 84188, 158288,
                                                                       40668, 40758, 87488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 162563, 0, 3,
                                                                       158288, 84288, 158438,
                                                                       40758, 40848, 87638,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 162788, 0, 3,
                                                                       158438, 84388, 158588,
                                                                       40848, 40938, 87788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 163013, 0, 3,
                                                                       158738, 84788, 158888,
                                                                       41118, 41208, 88238,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 163238, 0, 3,
                                                                       158888, 84888, 159038,
                                                                       41208, 41298, 88388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 163463, 0, 3,
                                                                       159038, 84988, 159188,
                                                                       41298, 41388, 88538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 163688, 0, 3,
                                                                       159188, 85088, 159338,
                                                                       41388, 41478, 88688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 163913, 0, 3,
                                                                       159338, 85188, 159488,
                                                                       41478, 41568, 88838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 164138, 0, 3,
                                                                       159488, 85288, 159638,
                                                                       41568, 41658, 88988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 164363, 0, 3,
                                                                       159638, 85388, 159788,
                                                                       41658, 41748, 89138,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 164588, 0, 3,
                                                                       159788, 85488, 159938,
                                                                       41748, 41838, 89288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 164813, 0, 3,
                                                                       159938, 85588, 160088,
                                                                       41838, 41928, 89438,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 165038, 0, 3,
                                                                       160088, 85688, 160238,
                                                                       41928, 42018, 89588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 165263, 0, 3,
                                                                       160238, 85788, 160388,
                                                                       42018, 42108, 89738,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 165488, 0, 3,
                                                                       160538, 86288, 160763,
                                                                       42288, 42414, 90308,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 165803, 0, 3,
                                                                       160763, 86438, 160988,
                                                                       42414, 42540, 90518,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 166118, 0, 3,
                                                                       160988, 86588, 161213,
                                                                       42540, 42666, 90728,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 166433, 0, 3,
                                                                       161213, 86738, 161438,
                                                                       42666, 42792, 90938,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 166748, 0, 3,
                                                                       161438, 86888, 161663,
                                                                       42792, 42918, 91148,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 167063, 0, 3,
                                                                       161663, 87038, 161888,
                                                                       42918, 43044, 91358,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 167378, 0, 3,
                                                                       161888, 87188, 162113,
                                                                       43044, 43170, 91568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 167693, 0, 3,
                                                                       162113, 87338, 162338,
                                                                       43170, 43296, 91778,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 168008, 0, 3,
                                                                       162338, 87488, 162563,
                                                                       43296, 43422, 91988,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 168323, 0, 3,
                                                                       162563, 87638, 162788,
                                                                       43422, 43548, 92198,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 168638, 0, 3,
                                                                       163013, 88238, 163238,
                                                                       43800, 43926, 92828,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 168953, 0, 3,
                                                                       163238, 88388, 163463,
                                                                       43926, 44052, 93038,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 169268, 0, 3,
                                                                       163463, 88538, 163688,
                                                                       44052, 44178, 93248,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 169583, 0, 3,
                                                                       163688, 88688, 163913,
                                                                       44178, 44304, 93458,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 169898, 0, 3,
                                                                       163913, 88838, 164138,
                                                                       44304, 44430, 93668,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 170213, 0, 3,
                                                                       164138, 88988, 164363,
                                                                       44430, 44556, 93878,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 170528, 0, 3,
                                                                       164363, 89138, 164588,
                                                                       44556, 44682, 94088,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 170843, 0, 3,
                                                                       164588, 89288, 164813,
                                                                       44682, 44808, 94298,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 171158, 0, 3,
                                                                       164813, 89438, 165038,
                                                                       44808, 44934, 94508,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 171473, 0, 3,
                                                                       165038, 89588, 165263,
                                                                       44934, 45060, 94718,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 171788, 0, 3,
                                                                       165488, 90308, 165803,
                                                                       45312, 45480, 95488,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 172208, 0, 3,
                                                                       165803, 90518, 166118,
                                                                       45480, 45648, 95768,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 172628, 0, 3,
                                                                       166118, 90728, 166433,
                                                                       45648, 45816, 96048,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 173048, 0, 3,
                                                                       166433, 90938, 166748,
                                                                       45816, 45984, 96328,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 173468, 0, 3,
                                                                       166748, 91148, 167063,
                                                                       45984, 46152, 96608,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 173888, 0, 3,
                                                                       167063, 91358, 167378,
                                                                       46152, 46320, 96888,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 174308, 0, 3,
                                                                       167378, 91568, 167693,
                                                                       46320, 46488, 97168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 174728, 0, 3,
                                                                       167693, 91778, 168008,
                                                                       46488, 46656, 97448,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 175148, 0, 3,
                                                                       168008, 91988, 168323,
                                                                       46656, 46824, 97728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 175568, 0, 3,
                                                                       168638, 92828, 168953,
                                                                       47160, 47328, 98568,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 175988, 0, 3,
                                                                       168953, 93038, 169268,
                                                                       47328, 47496, 98848,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 176408, 0, 3,
                                                                       169268, 93248, 169583,
                                                                       47496, 47664, 99128,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 176828, 0, 3,
                                                                       169583, 93458, 169898,
                                                                       47664, 47832, 99408,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 177248, 0, 3,
                                                                       169898, 93668, 170213,
                                                                       47832, 48000, 99688,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 177668, 0, 3,
                                                                       170213, 93878, 170528,
                                                                       48000, 48168, 99968,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 178088, 0, 3,
                                                                       170528, 94088, 170843,
                                                                       48168, 48336, 100248,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 178508, 0, 3,
                                                                       170843, 94298, 171158,
                                                                       48336, 48504, 100528,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 178928, 0, 3,
                                                                       171158, 94508, 171473,
                                                                       48504, 48672, 100808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 179348, 0, 3,
                                                                       171788, 95488, 172208,
                                                                       49008, 49224, 101808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 179888, 0, 3,
                                                                       172208, 95768, 172628,
                                                                       49224, 49440, 102168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 180428, 0, 3,
                                                                       172628, 96048, 173048,
                                                                       49440, 49656, 102528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 180968, 0, 3,
                                                                       173048, 96328, 173468,
                                                                       49656, 49872, 102888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 181508, 0, 3,
                                                                       173468, 96608, 173888,
                                                                       49872, 50088, 103248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 182048, 0, 3,
                                                                       173888, 96888, 174308,
                                                                       50088, 50304, 103608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 182588, 0, 3,
                                                                       174308, 97168, 174728,
                                                                       50304, 50520, 103968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 183128, 0, 3,
                                                                       174728, 97448, 175148,
                                                                       50520, 50736, 104328,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 183668, 0, 3,
                                                                       175568, 98568, 175988,
                                                                       51168, 51384, 105408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 184208, 0, 3,
                                                                       175988, 98848, 176408,
                                                                       51384, 51600, 105768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 184748, 0, 3,
                                                                       176408, 99128, 176828,
                                                                       51600, 51816, 106128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 185288, 0, 3,
                                                                       176828, 99408, 177248,
                                                                       51816, 52032, 106488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 185828, 0, 3,
                                                                       177248, 99688, 177668,
                                                                       52032, 52248, 106848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 186368, 0, 3,
                                                                       177668, 99968, 178088,
                                                                       52248, 52464, 107208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 186908, 0, 3,
                                                                       178088, 100248, 178508,
                                                                       52464, 52680, 107568,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 187448, 0, 3,
                                                                       178508, 100528, 178928,
                                                                       52680, 52896, 107928,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 187988, 0, 3,
                                                                       179348, 101808, 179888,
                                                                       53328, 53598, 109188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 188663, 0, 3,
                                                                       179888, 102168, 180428,
                                                                       53598, 53868, 109638,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 189338, 0, 3,
                                                                       180428, 102528, 180968,
                                                                       53868, 54138, 110088,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 190013, 0, 3,
                                                                       180968, 102888, 181508,
                                                                       54138, 54408, 110538,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 190688, 0, 3,
                                                                       181508, 103248, 182048,
                                                                       54408, 54678, 110988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 191363, 0, 3,
                                                                       182048, 103608, 182588,
                                                                       54678, 54948, 111438,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 192038, 0, 3,
                                                                       182588, 103968, 183128,
                                                                       54948, 55218, 111888,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 192713, 0, 3,
                                                                       183668, 105408, 184208,
                                                                       55758, 56028, 113238,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 193388, 0, 3,
                                                                       184208, 105768, 184748,
                                                                       56028, 56298, 113688,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 194063, 0, 3,
                                                                       184748, 106128, 185288,
                                                                       56298, 56568, 114138,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 194738, 0, 3,
                                                                       185288, 106488, 185828,
                                                                       56568, 56838, 114588,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 195413, 0, 3,
                                                                       185828, 106848, 186368,
                                                                       56838, 57108, 115038,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 196088, 0, 3,
                                                                       186368, 107208, 186908,
                                                                       57108, 57378, 115488,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 196763, 0, 3,
                                                                       186908, 107568, 187448,
                                                                       57378, 57648, 115938,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 197438, 0, 3,
                                                                       187988, 109188, 188663,
                                                                       58188, 58518, 117488,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 198263, 0, 3,
                                                                       188663, 109638, 189338,
                                                                       58518, 58848, 118038,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 199088, 0, 3,
                                                                       189338, 110088, 190013,
                                                                       58848, 59178, 118588,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 199913, 0, 3,
                                                                       190013, 110538, 190688,
                                                                       59178, 59508, 119138,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 200738, 0, 3,
                                                                       190688, 110988, 191363,
                                                                       59508, 59838, 119688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 201563, 0, 3,
                                                                       191363, 111438, 192038,
                                                                       59838, 60168, 120238,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 202388, 0, 3,
                                                                       192713, 113238, 193388,
                                                                       60828, 61158, 121888,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 203213, 0, 3,
                                                                       193388, 113688, 194063,
                                                                       61158, 61488, 122438,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 204038, 0, 3,
                                                                       194063, 114138, 194738,
                                                                       61488, 61818, 122988,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 204863, 0, 3,
                                                                       194738, 114588, 195413,
                                                                       61818, 62148, 123538,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 205688, 0, 3,
                                                                       195413, 115038, 196088,
                                                                       62148, 62478, 124088,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 206513, 0, 3,
                                                                       196088, 115488, 196763,
                                                                       62478, 62808, 124638,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 207338, 0, 3,
                                                                       197438, 117488, 198263,
                                                                       63468, 63864, 126508,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 208328, 0, 3,
                                                                       198263, 118038, 199088,
                                                                       63864, 64260, 127168,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 209318, 0, 3,
                                                                       199088, 118588, 199913,
                                                                       64260, 64656, 127828,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 210308, 0, 3,
                                                                       199913, 119138, 200738,
                                                                       64656, 65052, 128488,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 211298, 0, 3,
                                                                       200738, 119688, 201563,
                                                                       65052, 65448, 129148,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 212288, 0, 3,
                                                                       202388, 121888, 203213,
                                                                       66240, 66636, 131128,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 213278, 0, 3,
                                                                       203213, 122438, 204038,
                                                                       66636, 67032, 131788,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 214268, 0, 3,
                                                                       204038, 122988, 204863,
                                                                       67032, 67428, 132448,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 215258, 0, 3,
                                                                       204863, 123538, 205688,
                                                                       67428, 67824, 133108,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 216248, 0, 3,
                                                                       205688, 124088, 206513,
                                                                       67824, 68220, 133768,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 217238, 0, 3,
                                                                       207338, 126508, 208328,
                                                                       69012, 69480, 135988,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 218408, 0, 3,
                                                                       208328, 127168, 209318,
                                                                       69480, 69948, 136768,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 219578, 0, 3,
                                                                       209318, 127828, 210308,
                                                                       69948, 70416, 137548,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 220748, 0, 3,
                                                                       210308, 128488, 211298,
                                                                       70416, 70884, 138328,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 221918, 0, 3,
                                                                       212288, 131128, 213278,
                                                                       71820, 72288, 140668,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 223088, 0, 3,
                                                                       213278, 131788, 214268,
                                                                       72288, 72756, 141448,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 224258, 0, 3,
                                                                       214268, 132448, 215258,
                                                                       72756, 73224, 142228,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 225428, 0, 3,
                                                                       215258, 133108, 216248,
                                                                       73224, 73692, 143008,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 226598, 0, 3,
                                                                       217238, 135988, 218408,
                                                                       74628, 75174, 145608,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 227963, 0, 3,
                                                                       218408, 136768, 219578,
                                                                       75174, 75720, 146518,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 229328, 0, 3,
                                                                       219578, 137548, 220748,
                                                                       75720, 76266, 147428,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 230693, 0, 3,
                                                                       221918, 140668, 223088,
                                                                       77358, 77904, 150158,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 232058, 0, 3,
                                                                       223088, 141448, 224258,
                                                                       77904, 78450, 151068,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 233423, 0, 3,
                                                                       224258, 142228, 225428,
                                                                       78450, 78996, 151978,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234788, 3, 80088,
                                                                       80098, 152888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234809, 3, 80098,
                                                                       80108, 152903, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234830, 3, 80108,
                                                                       80118, 152918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234851, 3, 80118,
                                                                       80128, 152933, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234872, 3, 80128,
                                                                       80138, 152948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234893, 3, 80138,
                                                                       80148, 152963, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234914, 3, 80148,
                                                                       80158, 152978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234935, 3, 80158,
                                                                       80168, 152993, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234956, 3, 80168,
                                                                       80178, 153008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234977, 3, 80178,
                                                                       80188, 153023, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 234998, 3, 80188,
                                                                       80198, 153038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235019, 3, 80198,
                                                                       80208, 153053, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235040, 3, 80208,
                                                                       80218, 153068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235061, 3, 80218,
                                                                       80228, 153083, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235082, 3, 80228,
                                                                       80238, 153098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235103, 3, 80258,
                                                                       80268, 153113, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235124, 3, 80268,
                                                                       80278, 153128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235145, 3, 80278,
                                                                       80288, 153143, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235166, 3, 80288,
                                                                       80298, 153158, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235187, 3, 80298,
                                                                       80308, 153173, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235208, 3, 80308,
                                                                       80318, 153188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235229, 3, 80318,
                                                                       80328, 153203, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235250, 3, 80328,
                                                                       80338, 153218, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235271, 3, 80338,
                                                                       80348, 153233, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235292, 3, 80348,
                                                                       80358, 153248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235313, 3, 80358,
                                                                       80368, 153263, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235334, 3, 80368,
                                                                       80378, 153278, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235355, 3, 80378,
                                                                       80388, 153293, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235376, 3, 80388,
                                                                       80398, 153308, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 235397, 3, 80398,
                                                                       80408, 153323, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235418, 0, 3,
                                                                       234788, 152888, 234809,
                                                                       80428, 80458, 153338,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235481, 0, 3,
                                                                       234809, 152903, 234830,
                                                                       80458, 80488, 153383,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235544, 0, 3,
                                                                       234830, 152918, 234851,
                                                                       80488, 80518, 153428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235607, 0, 3,
                                                                       234851, 152933, 234872,
                                                                       80518, 80548, 153473,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235670, 0, 3,
                                                                       234872, 152948, 234893,
                                                                       80548, 80578, 153518,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235733, 0, 3,
                                                                       234893, 152963, 234914,
                                                                       80578, 80608, 153563,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235796, 0, 3,
                                                                       234914, 152978, 234935,
                                                                       80608, 80638, 153608,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235859, 0, 3,
                                                                       234935, 152993, 234956,
                                                                       80638, 80668, 153653,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235922, 0, 3,
                                                                       234956, 153008, 234977,
                                                                       80668, 80698, 153698,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 235985, 0, 3,
                                                                       234977, 153023, 234998,
                                                                       80698, 80728, 153743,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236048, 0, 3,
                                                                       234998, 153038, 235019,
                                                                       80728, 80758, 153788,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236111, 0, 3,
                                                                       235019, 153053, 235040,
                                                                       80758, 80788, 153833,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236174, 0, 3,
                                                                       235040, 153068, 235061,
                                                                       80788, 80818, 153878,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236237, 0, 3,
                                                                       235061, 153083, 235082,
                                                                       80818, 80848, 153923,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236300, 0, 3,
                                                                       235103, 153113, 235124,
                                                                       80908, 80938, 153968,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236363, 0, 3,
                                                                       235124, 153128, 235145,
                                                                       80938, 80968, 154013,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236426, 0, 3,
                                                                       235145, 153143, 235166,
                                                                       80968, 80998, 154058,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236489, 0, 3,
                                                                       235166, 153158, 235187,
                                                                       80998, 81028, 154103,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236552, 0, 3,
                                                                       235187, 153173, 235208,
                                                                       81028, 81058, 154148,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236615, 0, 3,
                                                                       235208, 153188, 235229,
                                                                       81058, 81088, 154193,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236678, 0, 3,
                                                                       235229, 153203, 235250,
                                                                       81088, 81118, 154238,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236741, 0, 3,
                                                                       235250, 153218, 235271,
                                                                       81118, 81148, 154283,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236804, 0, 3,
                                                                       235271, 153233, 235292,
                                                                       81148, 81178, 154328,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236867, 0, 3,
                                                                       235292, 153248, 235313,
                                                                       81178, 81208, 154373,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236930, 0, 3,
                                                                       235313, 153263, 235334,
                                                                       81208, 81238, 154418,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 236993, 0, 3,
                                                                       235334, 153278, 235355,
                                                                       81238, 81268, 154463,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 237056, 0, 3,
                                                                       235355, 153293, 235376,
                                                                       81268, 81298, 154508,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 237119, 0, 3,
                                                                       235376, 153308, 235397,
                                                                       81298, 81328, 154553,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 237182, 0, 3,
                                                                       235418, 153338, 235481,
                                                                       81388, 81448, 154598,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 237308, 0, 3,
                                                                       235481, 153383, 235544,
                                                                       81448, 81508, 154688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 237434, 0, 3,
                                                                       235544, 153428, 235607,
                                                                       81508, 81568, 154778,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 237560, 0, 3,
                                                                       235607, 153473, 235670,
                                                                       81568, 81628, 154868,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 237686, 0, 3,
                                                                       235670, 153518, 235733,
                                                                       81628, 81688, 154958,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 237812, 0, 3,
                                                                       235733, 153563, 235796,
                                                                       81688, 81748, 155048,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 237938, 0, 3,
                                                                       235796, 153608, 235859,
                                                                       81748, 81808, 155138,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 238064, 0, 3,
                                                                       235859, 153653, 235922,
                                                                       81808, 81868, 155228,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 238190, 0, 3,
                                                                       235922, 153698, 235985,
                                                                       81868, 81928, 155318,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 238316, 0, 3,
                                                                       235985, 153743, 236048,
                                                                       81928, 81988, 155408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 238442, 0, 3,
                                                                       236048, 153788, 236111,
                                                                       81988, 82048, 155498,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 238568, 0, 3,
                                                                       236111, 153833, 236174,
                                                                       82048, 82108, 155588,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 238694, 0, 3,
                                                                       236174, 153878, 236237,
                                                                       82108, 82168, 155678,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 238820, 0, 3,
                                                                       236300, 153968, 236363,
                                                                       82288, 82348, 155768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 238946, 0, 3,
                                                                       236363, 154013, 236426,
                                                                       82348, 82408, 155858,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 239072, 0, 3,
                                                                       236426, 154058, 236489,
                                                                       82408, 82468, 155948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 239198, 0, 3,
                                                                       236489, 154103, 236552,
                                                                       82468, 82528, 156038,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 239324, 0, 3,
                                                                       236552, 154148, 236615,
                                                                       82528, 82588, 156128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 239450, 0, 3,
                                                                       236615, 154193, 236678,
                                                                       82588, 82648, 156218,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 239576, 0, 3,
                                                                       236678, 154238, 236741,
                                                                       82648, 82708, 156308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 239702, 0, 3,
                                                                       236741, 154283, 236804,
                                                                       82708, 82768, 156398,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 239828, 0, 3,
                                                                       236804, 154328, 236867,
                                                                       82768, 82828, 156488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 239954, 0, 3,
                                                                       236867, 154373, 236930,
                                                                       82828, 82888, 156578,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 240080, 0, 3,
                                                                       236930, 154418, 236993,
                                                                       82888, 82948, 156668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 240206, 0, 3,
                                                                       236993, 154463, 237056,
                                                                       82948, 83008, 156758,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 240332, 0, 3,
                                                                       237056, 154508, 237119,
                                                                       83008, 83068, 156848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 240458, 0, 3,
                                                                       237182, 154598, 237308,
                                                                       83188, 83288, 156938,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 240668, 0, 3,
                                                                       237308, 154688, 237434,
                                                                       83288, 83388, 157088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 240878, 0, 3,
                                                                       237434, 154778, 237560,
                                                                       83388, 83488, 157238,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 241088, 0, 3,
                                                                       237560, 154868, 237686,
                                                                       83488, 83588, 157388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 241298, 0, 3,
                                                                       237686, 154958, 237812,
                                                                       83588, 83688, 157538,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 241508, 0, 3,
                                                                       237812, 155048, 237938,
                                                                       83688, 83788, 157688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 241718, 0, 3,
                                                                       237938, 155138, 238064,
                                                                       83788, 83888, 157838,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 241928, 0, 3,
                                                                       238064, 155228, 238190,
                                                                       83888, 83988, 157988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 242138, 0, 3,
                                                                       238190, 155318, 238316,
                                                                       83988, 84088, 158138,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 242348, 0, 3,
                                                                       238316, 155408, 238442,
                                                                       84088, 84188, 158288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 242558, 0, 3,
                                                                       238442, 155498, 238568,
                                                                       84188, 84288, 158438,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 242768, 0, 3,
                                                                       238568, 155588, 238694,
                                                                       84288, 84388, 158588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 242978, 0, 3,
                                                                       238820, 155768, 238946,
                                                                       84588, 84688, 158738,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 243188, 0, 3,
                                                                       238946, 155858, 239072,
                                                                       84688, 84788, 158888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 243398, 0, 3,
                                                                       239072, 155948, 239198,
                                                                       84788, 84888, 159038,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 243608, 0, 3,
                                                                       239198, 156038, 239324,
                                                                       84888, 84988, 159188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 243818, 0, 3,
                                                                       239324, 156128, 239450,
                                                                       84988, 85088, 159338,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 244028, 0, 3,
                                                                       239450, 156218, 239576,
                                                                       85088, 85188, 159488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 244238, 0, 3,
                                                                       239576, 156308, 239702,
                                                                       85188, 85288, 159638,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 244448, 0, 3,
                                                                       239702, 156398, 239828,
                                                                       85288, 85388, 159788,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 244658, 0, 3,
                                                                       239828, 156488, 239954,
                                                                       85388, 85488, 159938,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 244868, 0, 3,
                                                                       239954, 156578, 240080,
                                                                       85488, 85588, 160088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 245078, 0, 3,
                                                                       240080, 156668, 240206,
                                                                       85588, 85688, 160238,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 245288, 0, 3,
                                                                       240206, 156758, 240332,
                                                                       85688, 85788, 160388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 245498, 0, 3,
                                                                       240458, 156938, 240668,
                                                                       85988, 86138, 160538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 245813, 0, 3,
                                                                       240668, 157088, 240878,
                                                                       86138, 86288, 160763,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 246128, 0, 3,
                                                                       240878, 157238, 241088,
                                                                       86288, 86438, 160988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 246443, 0, 3,
                                                                       241088, 157388, 241298,
                                                                       86438, 86588, 161213,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 246758, 0, 3,
                                                                       241298, 157538, 241508,
                                                                       86588, 86738, 161438,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 247073, 0, 3,
                                                                       241508, 157688, 241718,
                                                                       86738, 86888, 161663,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 247388, 0, 3,
                                                                       241718, 157838, 241928,
                                                                       86888, 87038, 161888,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 247703, 0, 3,
                                                                       241928, 157988, 242138,
                                                                       87038, 87188, 162113,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 248018, 0, 3,
                                                                       242138, 158138, 242348,
                                                                       87188, 87338, 162338,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 248333, 0, 3,
                                                                       242348, 158288, 242558,
                                                                       87338, 87488, 162563,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 248648, 0, 3,
                                                                       242558, 158438, 242768,
                                                                       87488, 87638, 162788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 248963, 0, 3,
                                                                       242978, 158738, 243188,
                                                                       87938, 88088, 163013,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 249278, 0, 3,
                                                                       243188, 158888, 243398,
                                                                       88088, 88238, 163238,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 249593, 0, 3,
                                                                       243398, 159038, 243608,
                                                                       88238, 88388, 163463,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 249908, 0, 3,
                                                                       243608, 159188, 243818,
                                                                       88388, 88538, 163688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 250223, 0, 3,
                                                                       243818, 159338, 244028,
                                                                       88538, 88688, 163913,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 250538, 0, 3,
                                                                       244028, 159488, 244238,
                                                                       88688, 88838, 164138,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 250853, 0, 3,
                                                                       244238, 159638, 244448,
                                                                       88838, 88988, 164363,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 251168, 0, 3,
                                                                       244448, 159788, 244658,
                                                                       88988, 89138, 164588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 251483, 0, 3,
                                                                       244658, 159938, 244868,
                                                                       89138, 89288, 164813,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 251798, 0, 3,
                                                                       244868, 160088, 245078,
                                                                       89288, 89438, 165038,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 252113, 0, 3,
                                                                       245078, 160238, 245288,
                                                                       89438, 89588, 165263,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 252428, 0, 3,
                                                                       245498, 160538, 245813,
                                                                       89888, 90098, 165488,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 252869, 0, 3,
                                                                       245813, 160763, 246128,
                                                                       90098, 90308, 165803,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 253310, 0, 3,
                                                                       246128, 160988, 246443,
                                                                       90308, 90518, 166118,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 253751, 0, 3,
                                                                       246443, 161213, 246758,
                                                                       90518, 90728, 166433,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 254192, 0, 3,
                                                                       246758, 161438, 247073,
                                                                       90728, 90938, 166748,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 254633, 0, 3,
                                                                       247073, 161663, 247388,
                                                                       90938, 91148, 167063,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 255074, 0, 3,
                                                                       247388, 161888, 247703,
                                                                       91148, 91358, 167378,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 255515, 0, 3,
                                                                       247703, 162113, 248018,
                                                                       91358, 91568, 167693,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 255956, 0, 3,
                                                                       248018, 162338, 248333,
                                                                       91568, 91778, 168008,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 256397, 0, 3,
                                                                       248333, 162563, 248648,
                                                                       91778, 91988, 168323,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 256838, 0, 3,
                                                                       248963, 163013, 249278,
                                                                       92408, 92618, 168638,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 257279, 0, 3,
                                                                       249278, 163238, 249593,
                                                                       92618, 92828, 168953,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 257720, 0, 3,
                                                                       249593, 163463, 249908,
                                                                       92828, 93038, 169268,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 258161, 0, 3,
                                                                       249908, 163688, 250223,
                                                                       93038, 93248, 169583,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 258602, 0, 3,
                                                                       250223, 163913, 250538,
                                                                       93248, 93458, 169898,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 259043, 0, 3,
                                                                       250538, 164138, 250853,
                                                                       93458, 93668, 170213,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 259484, 0, 3,
                                                                       250853, 164363, 251168,
                                                                       93668, 93878, 170528,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 259925, 0, 3,
                                                                       251168, 164588, 251483,
                                                                       93878, 94088, 170843,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 260366, 0, 3,
                                                                       251483, 164813, 251798,
                                                                       94088, 94298, 171158,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 260807, 0, 3,
                                                                       251798, 165038, 252113,
                                                                       94298, 94508, 171473,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 261248, 0, 3,
                                                                       252428, 165488, 252869,
                                                                       94928, 95208, 171788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 261836, 0, 3,
                                                                       252869, 165803, 253310,
                                                                       95208, 95488, 172208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 262424, 0, 3,
                                                                       253310, 166118, 253751,
                                                                       95488, 95768, 172628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 263012, 0, 3,
                                                                       253751, 166433, 254192,
                                                                       95768, 96048, 173048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 263600, 0, 3,
                                                                       254192, 166748, 254633,
                                                                       96048, 96328, 173468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 264188, 0, 3,
                                                                       254633, 167063, 255074,
                                                                       96328, 96608, 173888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 264776, 0, 3,
                                                                       255074, 167378, 255515,
                                                                       96608, 96888, 174308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 265364, 0, 3,
                                                                       255515, 167693, 255956,
                                                                       96888, 97168, 174728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 265952, 0, 3,
                                                                       255956, 168008, 256397,
                                                                       97168, 97448, 175148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 266540, 0, 3,
                                                                       256838, 168638, 257279,
                                                                       98008, 98288, 175568,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 267128, 0, 3,
                                                                       257279, 168953, 257720,
                                                                       98288, 98568, 175988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 267716, 0, 3,
                                                                       257720, 169268, 258161,
                                                                       98568, 98848, 176408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 268304, 0, 3,
                                                                       258161, 169583, 258602,
                                                                       98848, 99128, 176828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 268892, 0, 3,
                                                                       258602, 169898, 259043,
                                                                       99128, 99408, 177248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 269480, 0, 3,
                                                                       259043, 170213, 259484,
                                                                       99408, 99688, 177668,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 270068, 0, 3,
                                                                       259484, 170528, 259925,
                                                                       99688, 99968, 178088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 270656, 0, 3,
                                                                       259925, 170843, 260366,
                                                                       99968, 100248, 178508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 271244, 0, 3,
                                                                       260366, 171158, 260807,
                                                                       100248, 100528, 178928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 271832, 0, 3,
                                                                       261248, 171788, 261836,
                                                                       101088, 101448, 179348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 272588, 0, 3,
                                                                       261836, 172208, 262424,
                                                                       101448, 101808, 179888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 273344, 0, 3,
                                                                       262424, 172628, 263012,
                                                                       101808, 102168, 180428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 274100, 0, 3,
                                                                       263012, 173048, 263600,
                                                                       102168, 102528, 180968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 274856, 0, 3,
                                                                       263600, 173468, 264188,
                                                                       102528, 102888, 181508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 275612, 0, 3,
                                                                       264188, 173888, 264776,
                                                                       102888, 103248, 182048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 276368, 0, 3,
                                                                       264776, 174308, 265364,
                                                                       103248, 103608, 182588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 277124, 0, 3,
                                                                       265364, 174728, 265952,
                                                                       103608, 103968, 183128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 277880, 0, 3,
                                                                       266540, 175568, 267128,
                                                                       104688, 105048, 183668,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 278636, 0, 3,
                                                                       267128, 175988, 267716,
                                                                       105048, 105408, 184208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 279392, 0, 3,
                                                                       267716, 176408, 268304,
                                                                       105408, 105768, 184748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 280148, 0, 3,
                                                                       268304, 176828, 268892,
                                                                       105768, 106128, 185288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 280904, 0, 3,
                                                                       268892, 177248, 269480,
                                                                       106128, 106488, 185828,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 281660, 0, 3,
                                                                       269480, 177668, 270068,
                                                                       106488, 106848, 186368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 282416, 0, 3,
                                                                       270068, 178088, 270656,
                                                                       106848, 107208, 186908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 283172, 0, 3,
                                                                       270656, 178508, 271244,
                                                                       107208, 107568, 187448,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 283928, 0, 3,
                                                                       271832, 179348, 272588,
                                                                       108288, 108738, 187988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 284873, 0, 3,
                                                                       272588, 179888, 273344,
                                                                       108738, 109188, 188663,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 285818, 0, 3,
                                                                       273344, 180428, 274100,
                                                                       109188, 109638, 189338,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 286763, 0, 3,
                                                                       274100, 180968, 274856,
                                                                       109638, 110088, 190013,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 287708, 0, 3,
                                                                       274856, 181508, 275612,
                                                                       110088, 110538, 190688,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 288653, 0, 3,
                                                                       275612, 182048, 276368,
                                                                       110538, 110988, 191363,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 289598, 0, 3,
                                                                       276368, 182588, 277124,
                                                                       110988, 111438, 192038,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 290543, 0, 3,
                                                                       277880, 183668, 278636,
                                                                       112338, 112788, 192713,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 291488, 0, 3,
                                                                       278636, 184208, 279392,
                                                                       112788, 113238, 193388,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 292433, 0, 3,
                                                                       279392, 184748, 280148,
                                                                       113238, 113688, 194063,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 293378, 0, 3,
                                                                       280148, 185288, 280904,
                                                                       113688, 114138, 194738,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 294323, 0, 3,
                                                                       280904, 185828, 281660,
                                                                       114138, 114588, 195413,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 295268, 0, 3,
                                                                       281660, 186368, 282416,
                                                                       114588, 115038, 196088,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 296213, 0, 3,
                                                                       282416, 186908, 283172,
                                                                       115038, 115488, 196763,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 297158, 0, 3,
                                                                       283928, 187988, 284873,
                                                                       116388, 116938, 197438,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 298313, 0, 3,
                                                                       284873, 188663, 285818,
                                                                       116938, 117488, 198263,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 299468, 0, 3,
                                                                       285818, 189338, 286763,
                                                                       117488, 118038, 199088,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 300623, 0, 3,
                                                                       286763, 190013, 287708,
                                                                       118038, 118588, 199913,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 301778, 0, 3,
                                                                       287708, 190688, 288653,
                                                                       118588, 119138, 200738,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 302933, 0, 3,
                                                                       288653, 191363, 289598,
                                                                       119138, 119688, 201563,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 304088, 0, 3,
                                                                       290543, 192713, 291488,
                                                                       120788, 121338, 202388,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 305243, 0, 3,
                                                                       291488, 193388, 292433,
                                                                       121338, 121888, 203213,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 306398, 0, 3,
                                                                       292433, 194063, 293378,
                                                                       121888, 122438, 204038,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 307553, 0, 3,
                                                                       293378, 194738, 294323,
                                                                       122438, 122988, 204863,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 308708, 0, 3,
                                                                       294323, 195413, 295268,
                                                                       122988, 123538, 205688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 309863, 0, 3,
                                                                       295268, 196088, 296213,
                                                                       123538, 124088, 206513,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 311018, 0, 3,
                                                                       297158, 197438, 298313,
                                                                       125188, 125848, 207338,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 312404, 0, 3,
                                                                       298313, 198263, 299468,
                                                                       125848, 126508, 208328,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 313790, 0, 3,
                                                                       299468, 199088, 300623,
                                                                       126508, 127168, 209318,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 315176, 0, 3,
                                                                       300623, 199913, 301778,
                                                                       127168, 127828, 210308,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 316562, 0, 3,
                                                                       301778, 200738, 302933,
                                                                       127828, 128488, 211298,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 317948, 0, 3,
                                                                       304088, 202388, 305243,
                                                                       129808, 130468, 212288,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 319334, 0, 3,
                                                                       305243, 203213, 306398,
                                                                       130468, 131128, 213278,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 320720, 0, 3,
                                                                       306398, 204038, 307553,
                                                                       131128, 131788, 214268,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 322106, 0, 3,
                                                                       307553, 204863, 308708,
                                                                       131788, 132448, 215258,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 323492, 0, 3,
                                                                       308708, 205688, 309863,
                                                                       132448, 133108, 216248,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 324878, 0, 3,
                                                                       311018, 207338, 312404,
                                                                       134428, 135208, 217238,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 326516, 0, 3,
                                                                       312404, 208328, 313790,
                                                                       135208, 135988, 218408,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 328154, 0, 3,
                                                                       313790, 209318, 315176,
                                                                       135988, 136768, 219578,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 329792, 0, 3,
                                                                       315176, 210308, 316562,
                                                                       136768, 137548, 220748,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 331430, 0, 3,
                                                                       317948, 212288, 319334,
                                                                       139108, 139888, 221918,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 333068, 0, 3,
                                                                       319334, 213278, 320720,
                                                                       139888, 140668, 223088,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 334706, 0, 3,
                                                                       320720, 214268, 322106,
                                                                       140668, 141448, 224258,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 336344, 0, 3,
                                                                       322106, 215258, 323492,
                                                                       141448, 142228, 225428,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 337982, 0, 3,
                                                                       324878, 217238, 326516,
                                                                       143788, 144698, 226598,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 339893, 0, 3,
                                                                       326516, 218408, 328154,
                                                                       144698, 145608, 227963,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 341804, 0, 3,
                                                                       328154, 219578, 329792,
                                                                       145608, 146518, 229328,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 343715, 0, 3,
                                                                       331430, 221918, 333068,
                                                                       148338, 149248, 230693,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 345626, 0, 3,
                                                                       333068, 223088, 334706,
                                                                       149248, 150158, 232058,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 347537, 0, 3,
                                                                       334706, 224258, 336344,
                                                                       150158, 151068, 233423,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349448, 3, 152888,
                                                                       152903, 234830, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349476, 3, 152903,
                                                                       152918, 234851, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349504, 3, 152918,
                                                                       152933, 234872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349532, 3, 152933,
                                                                       152948, 234893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349560, 3, 152948,
                                                                       152963, 234914, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349588, 3, 152963,
                                                                       152978, 234935, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349616, 3, 152978,
                                                                       152993, 234956, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349644, 3, 152993,
                                                                       153008, 234977, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349672, 3, 153008,
                                                                       153023, 234998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349700, 3, 153023,
                                                                       153038, 235019, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349728, 3, 153038,
                                                                       153053, 235040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349756, 3, 153053,
                                                                       153068, 235061, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349784, 3, 153068,
                                                                       153083, 235082, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349812, 3, 153113,
                                                                       153128, 235145, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349840, 3, 153128,
                                                                       153143, 235166, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349868, 3, 153143,
                                                                       153158, 235187, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349896, 3, 153158,
                                                                       153173, 235208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349924, 3, 153173,
                                                                       153188, 235229, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349952, 3, 153188,
                                                                       153203, 235250, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 349980, 3, 153203,
                                                                       153218, 235271, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 350008, 3, 153218,
                                                                       153233, 235292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 350036, 3, 153233,
                                                                       153248, 235313, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 350064, 3, 153248,
                                                                       153263, 235334, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 350092, 3, 153263,
                                                                       153278, 235355, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 350120, 3, 153278,
                                                                       153293, 235376, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 350148, 3, 153293,
                                                                       153308, 235397, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350176, 0, 3,
                                                                       349448, 234830, 349476,
                                                                       153338, 153383, 235544,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350260, 0, 3,
                                                                       349476, 234851, 349504,
                                                                       153383, 153428, 235607,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350344, 0, 3,
                                                                       349504, 234872, 349532,
                                                                       153428, 153473, 235670,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350428, 0, 3,
                                                                       349532, 234893, 349560,
                                                                       153473, 153518, 235733,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350512, 0, 3,
                                                                       349560, 234914, 349588,
                                                                       153518, 153563, 235796,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350596, 0, 3,
                                                                       349588, 234935, 349616,
                                                                       153563, 153608, 235859,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350680, 0, 3,
                                                                       349616, 234956, 349644,
                                                                       153608, 153653, 235922,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350764, 0, 3,
                                                                       349644, 234977, 349672,
                                                                       153653, 153698, 235985,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350848, 0, 3,
                                                                       349672, 234998, 349700,
                                                                       153698, 153743, 236048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 350932, 0, 3,
                                                                       349700, 235019, 349728,
                                                                       153743, 153788, 236111,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351016, 0, 3,
                                                                       349728, 235040, 349756,
                                                                       153788, 153833, 236174,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351100, 0, 3,
                                                                       349756, 235061, 349784,
                                                                       153833, 153878, 236237,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351184, 0, 3,
                                                                       349812, 235145, 349840,
                                                                       153968, 154013, 236426,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351268, 0, 3,
                                                                       349840, 235166, 349868,
                                                                       154013, 154058, 236489,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351352, 0, 3,
                                                                       349868, 235187, 349896,
                                                                       154058, 154103, 236552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351436, 0, 3,
                                                                       349896, 235208, 349924,
                                                                       154103, 154148, 236615,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351520, 0, 3,
                                                                       349924, 235229, 349952,
                                                                       154148, 154193, 236678,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351604, 0, 3,
                                                                       349952, 235250, 349980,
                                                                       154193, 154238, 236741,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351688, 0, 3,
                                                                       349980, 235271, 350008,
                                                                       154238, 154283, 236804,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351772, 0, 3,
                                                                       350008, 235292, 350036,
                                                                       154283, 154328, 236867,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351856, 0, 3,
                                                                       350036, 235313, 350064,
                                                                       154328, 154373, 236930,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 351940, 0, 3,
                                                                       350064, 235334, 350092,
                                                                       154373, 154418, 236993,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 352024, 0, 3,
                                                                       350092, 235355, 350120,
                                                                       154418, 154463, 237056,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 352108, 0, 3,
                                                                       350120, 235376, 350148,
                                                                       154463, 154508, 237119,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 352192, 0, 3,
                                                                       350176, 235544, 350260,
                                                                       154598, 154688, 237434,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 352360, 0, 3,
                                                                       350260, 235607, 350344,
                                                                       154688, 154778, 237560,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 352528, 0, 3,
                                                                       350344, 235670, 350428,
                                                                       154778, 154868, 237686,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 352696, 0, 3,
                                                                       350428, 235733, 350512,
                                                                       154868, 154958, 237812,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 352864, 0, 3,
                                                                       350512, 235796, 350596,
                                                                       154958, 155048, 237938,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 353032, 0, 3,
                                                                       350596, 235859, 350680,
                                                                       155048, 155138, 238064,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 353200, 0, 3,
                                                                       350680, 235922, 350764,
                                                                       155138, 155228, 238190,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 353368, 0, 3,
                                                                       350764, 235985, 350848,
                                                                       155228, 155318, 238316,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 353536, 0, 3,
                                                                       350848, 236048, 350932,
                                                                       155318, 155408, 238442,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 353704, 0, 3,
                                                                       350932, 236111, 351016,
                                                                       155408, 155498, 238568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 353872, 0, 3,
                                                                       351016, 236174, 351100,
                                                                       155498, 155588, 238694,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 354040, 0, 3,
                                                                       351184, 236426, 351268,
                                                                       155768, 155858, 239072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 354208, 0, 3,
                                                                       351268, 236489, 351352,
                                                                       155858, 155948, 239198,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 354376, 0, 3,
                                                                       351352, 236552, 351436,
                                                                       155948, 156038, 239324,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 354544, 0, 3,
                                                                       351436, 236615, 351520,
                                                                       156038, 156128, 239450,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 354712, 0, 3,
                                                                       351520, 236678, 351604,
                                                                       156128, 156218, 239576,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 354880, 0, 3,
                                                                       351604, 236741, 351688,
                                                                       156218, 156308, 239702,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 355048, 0, 3,
                                                                       351688, 236804, 351772,
                                                                       156308, 156398, 239828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 355216, 0, 3,
                                                                       351772, 236867, 351856,
                                                                       156398, 156488, 239954,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 355384, 0, 3,
                                                                       351856, 236930, 351940,
                                                                       156488, 156578, 240080,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 355552, 0, 3,
                                                                       351940, 236993, 352024,
                                                                       156578, 156668, 240206,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 355720, 0, 3,
                                                                       352024, 237056, 352108,
                                                                       156668, 156758, 240332,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 355888, 0, 3,
                                                                       352192, 237434, 352360,
                                                                       156938, 157088, 240878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 356168, 0, 3,
                                                                       352360, 237560, 352528,
                                                                       157088, 157238, 241088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 356448, 0, 3,
                                                                       352528, 237686, 352696,
                                                                       157238, 157388, 241298,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 356728, 0, 3,
                                                                       352696, 237812, 352864,
                                                                       157388, 157538, 241508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 357008, 0, 3,
                                                                       352864, 237938, 353032,
                                                                       157538, 157688, 241718,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 357288, 0, 3,
                                                                       353032, 238064, 353200,
                                                                       157688, 157838, 241928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 357568, 0, 3,
                                                                       353200, 238190, 353368,
                                                                       157838, 157988, 242138,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 357848, 0, 3,
                                                                       353368, 238316, 353536,
                                                                       157988, 158138, 242348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 358128, 0, 3,
                                                                       353536, 238442, 353704,
                                                                       158138, 158288, 242558,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 358408, 0, 3,
                                                                       353704, 238568, 353872,
                                                                       158288, 158438, 242768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 358688, 0, 3,
                                                                       354040, 239072, 354208,
                                                                       158738, 158888, 243398,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 358968, 0, 3,
                                                                       354208, 239198, 354376,
                                                                       158888, 159038, 243608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 359248, 0, 3,
                                                                       354376, 239324, 354544,
                                                                       159038, 159188, 243818,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 359528, 0, 3,
                                                                       354544, 239450, 354712,
                                                                       159188, 159338, 244028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 359808, 0, 3,
                                                                       354712, 239576, 354880,
                                                                       159338, 159488, 244238,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 360088, 0, 3,
                                                                       354880, 239702, 355048,
                                                                       159488, 159638, 244448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 360368, 0, 3,
                                                                       355048, 239828, 355216,
                                                                       159638, 159788, 244658,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 360648, 0, 3,
                                                                       355216, 239954, 355384,
                                                                       159788, 159938, 244868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 360928, 0, 3,
                                                                       355384, 240080, 355552,
                                                                       159938, 160088, 245078,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 361208, 0, 3,
                                                                       355552, 240206, 355720,
                                                                       160088, 160238, 245288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 361488, 0, 3,
                                                                       355888, 240878, 356168,
                                                                       160538, 160763, 246128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 361908, 0, 3,
                                                                       356168, 241088, 356448,
                                                                       160763, 160988, 246443,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 362328, 0, 3,
                                                                       356448, 241298, 356728,
                                                                       160988, 161213, 246758,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 362748, 0, 3,
                                                                       356728, 241508, 357008,
                                                                       161213, 161438, 247073,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 363168, 0, 3,
                                                                       357008, 241718, 357288,
                                                                       161438, 161663, 247388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 363588, 0, 3,
                                                                       357288, 241928, 357568,
                                                                       161663, 161888, 247703,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 364008, 0, 3,
                                                                       357568, 242138, 357848,
                                                                       161888, 162113, 248018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 364428, 0, 3,
                                                                       357848, 242348, 358128,
                                                                       162113, 162338, 248333,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 364848, 0, 3,
                                                                       358128, 242558, 358408,
                                                                       162338, 162563, 248648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 365268, 0, 3,
                                                                       358688, 243398, 358968,
                                                                       163013, 163238, 249593,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 365688, 0, 3,
                                                                       358968, 243608, 359248,
                                                                       163238, 163463, 249908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 366108, 0, 3,
                                                                       359248, 243818, 359528,
                                                                       163463, 163688, 250223,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 366528, 0, 3,
                                                                       359528, 244028, 359808,
                                                                       163688, 163913, 250538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 366948, 0, 3,
                                                                       359808, 244238, 360088,
                                                                       163913, 164138, 250853,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 367368, 0, 3,
                                                                       360088, 244448, 360368,
                                                                       164138, 164363, 251168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 367788, 0, 3,
                                                                       360368, 244658, 360648,
                                                                       164363, 164588, 251483,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 368208, 0, 3,
                                                                       360648, 244868, 360928,
                                                                       164588, 164813, 251798,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 368628, 0, 3,
                                                                       360928, 245078, 361208,
                                                                       164813, 165038, 252113,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 369048, 0, 3,
                                                                       361488, 246128, 361908,
                                                                       165488, 165803, 253310,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 369636, 0, 3,
                                                                       361908, 246443, 362328,
                                                                       165803, 166118, 253751,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 370224, 0, 3,
                                                                       362328, 246758, 362748,
                                                                       166118, 166433, 254192,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 370812, 0, 3,
                                                                       362748, 247073, 363168,
                                                                       166433, 166748, 254633,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 371400, 0, 3,
                                                                       363168, 247388, 363588,
                                                                       166748, 167063, 255074,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 371988, 0, 3,
                                                                       363588, 247703, 364008,
                                                                       167063, 167378, 255515,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 372576, 0, 3,
                                                                       364008, 248018, 364428,
                                                                       167378, 167693, 255956,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 373164, 0, 3,
                                                                       364428, 248333, 364848,
                                                                       167693, 168008, 256397,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 373752, 0, 3,
                                                                       365268, 249593, 365688,
                                                                       168638, 168953, 257720,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 374340, 0, 3,
                                                                       365688, 249908, 366108,
                                                                       168953, 169268, 258161,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 374928, 0, 3,
                                                                       366108, 250223, 366528,
                                                                       169268, 169583, 258602,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 375516, 0, 3,
                                                                       366528, 250538, 366948,
                                                                       169583, 169898, 259043,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 376104, 0, 3,
                                                                       366948, 250853, 367368,
                                                                       169898, 170213, 259484,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 376692, 0, 3,
                                                                       367368, 251168, 367788,
                                                                       170213, 170528, 259925,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 377280, 0, 3,
                                                                       367788, 251483, 368208,
                                                                       170528, 170843, 260366,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 377868, 0, 3,
                                                                       368208, 251798, 368628,
                                                                       170843, 171158, 260807,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 378456, 0, 3,
                                                                       369048, 253310, 369636,
                                                                       171788, 172208, 262424,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 379240, 0, 3,
                                                                       369636, 253751, 370224,
                                                                       172208, 172628, 263012,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 380024, 0, 3,
                                                                       370224, 254192, 370812,
                                                                       172628, 173048, 263600,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 380808, 0, 3,
                                                                       370812, 254633, 371400,
                                                                       173048, 173468, 264188,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 381592, 0, 3,
                                                                       371400, 255074, 371988,
                                                                       173468, 173888, 264776,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 382376, 0, 3,
                                                                       371988, 255515, 372576,
                                                                       173888, 174308, 265364,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 383160, 0, 3,
                                                                       372576, 255956, 373164,
                                                                       174308, 174728, 265952,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 383944, 0, 3,
                                                                       373752, 257720, 374340,
                                                                       175568, 175988, 267716,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 384728, 0, 3,
                                                                       374340, 258161, 374928,
                                                                       175988, 176408, 268304,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 385512, 0, 3,
                                                                       374928, 258602, 375516,
                                                                       176408, 176828, 268892,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 386296, 0, 3,
                                                                       375516, 259043, 376104,
                                                                       176828, 177248, 269480,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 387080, 0, 3,
                                                                       376104, 259484, 376692,
                                                                       177248, 177668, 270068,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 387864, 0, 3,
                                                                       376692, 259925, 377280,
                                                                       177668, 178088, 270656,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 388648, 0, 3,
                                                                       377280, 260366, 377868,
                                                                       178088, 178508, 271244,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 389432, 0, 3,
                                                                       378456, 262424, 379240,
                                                                       179348, 179888, 273344,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 390440, 0, 3,
                                                                       379240, 263012, 380024,
                                                                       179888, 180428, 274100,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 391448, 0, 3,
                                                                       380024, 263600, 380808,
                                                                       180428, 180968, 274856,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 392456, 0, 3,
                                                                       380808, 264188, 381592,
                                                                       180968, 181508, 275612,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 393464, 0, 3,
                                                                       381592, 264776, 382376,
                                                                       181508, 182048, 276368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 394472, 0, 3,
                                                                       382376, 265364, 383160,
                                                                       182048, 182588, 277124,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 395480, 0, 3,
                                                                       383944, 267716, 384728,
                                                                       183668, 184208, 279392,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 396488, 0, 3,
                                                                       384728, 268304, 385512,
                                                                       184208, 184748, 280148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 397496, 0, 3,
                                                                       385512, 268892, 386296,
                                                                       184748, 185288, 280904,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 398504, 0, 3,
                                                                       386296, 269480, 387080,
                                                                       185288, 185828, 281660,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 399512, 0, 3,
                                                                       387080, 270068, 387864,
                                                                       185828, 186368, 282416,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 400520, 0, 3,
                                                                       387864, 270656, 388648,
                                                                       186368, 186908, 283172,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 401528, 0, 3,
                                                                       389432, 273344, 390440,
                                                                       187988, 188663, 285818,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 402788, 0, 3,
                                                                       390440, 274100, 391448,
                                                                       188663, 189338, 286763,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 404048, 0, 3,
                                                                       391448, 274856, 392456,
                                                                       189338, 190013, 287708,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 405308, 0, 3,
                                                                       392456, 275612, 393464,
                                                                       190013, 190688, 288653,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 406568, 0, 3,
                                                                       393464, 276368, 394472,
                                                                       190688, 191363, 289598,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 407828, 0, 3,
                                                                       395480, 279392, 396488,
                                                                       192713, 193388, 292433,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 409088, 0, 3,
                                                                       396488, 280148, 397496,
                                                                       193388, 194063, 293378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 410348, 0, 3,
                                                                       397496, 280904, 398504,
                                                                       194063, 194738, 294323,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 411608, 0, 3,
                                                                       398504, 281660, 399512,
                                                                       194738, 195413, 295268,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 412868, 0, 3,
                                                                       399512, 282416, 400520,
                                                                       195413, 196088, 296213,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 414128, 0, 3,
                                                                       401528, 285818, 402788,
                                                                       197438, 198263, 299468,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 415668, 0, 3,
                                                                       402788, 286763, 404048,
                                                                       198263, 199088, 300623,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 417208, 0, 3,
                                                                       404048, 287708, 405308,
                                                                       199088, 199913, 301778,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 418748, 0, 3,
                                                                       405308, 288653, 406568,
                                                                       199913, 200738, 302933,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 420288, 0, 3,
                                                                       407828, 292433, 409088,
                                                                       202388, 203213, 306398,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 421828, 0, 3,
                                                                       409088, 293378, 410348,
                                                                       203213, 204038, 307553,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 423368, 0, 3,
                                                                       410348, 294323, 411608,
                                                                       204038, 204863, 308708,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 424908, 0, 3,
                                                                       411608, 295268, 412868,
                                                                       204863, 205688, 309863,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 426448, 0, 3,
                                                                       414128, 299468, 415668,
                                                                       207338, 208328, 313790,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 428296, 0, 3,
                                                                       415668, 300623, 417208,
                                                                       208328, 209318, 315176,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 430144, 0, 3,
                                                                       417208, 301778, 418748,
                                                                       209318, 210308, 316562,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 431992, 0, 3,
                                                                       420288, 306398, 421828,
                                                                       212288, 213278, 320720,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 433840, 0, 3,
                                                                       421828, 307553, 423368,
                                                                       213278, 214268, 322106,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 435688, 0, 3,
                                                                       423368, 308708, 424908,
                                                                       214268, 215258, 323492,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 437536, 0, 3,
                                                                       426448, 313790, 428296,
                                                                       217238, 218408, 328154,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 439720, 0, 3,
                                                                       428296, 315176, 430144,
                                                                       218408, 219578, 329792,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 441904, 0, 3,
                                                                       431992, 320720, 433840,
                                                                       221918, 223088, 334706,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 444088, 0, 3,
                                                                       433840, 322106, 435688,
                                                                       223088, 224258, 336344,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 446272, 0, 3,
                                                                       437536, 328154, 439720,
                                                                       226598, 227963, 341804,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 448820, 0, 3,
                                                                       441904, 334706, 444088,
                                                                       230693, 232058, 347537,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451368, 3, 234788,
                                                                       234809, 349448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451404, 3, 234809,
                                                                       234830, 349476, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451440, 3, 234830,
                                                                       234851, 349504, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451476, 3, 234851,
                                                                       234872, 349532, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451512, 3, 234872,
                                                                       234893, 349560, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451548, 3, 234893,
                                                                       234914, 349588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451584, 3, 234914,
                                                                       234935, 349616, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451620, 3, 234935,
                                                                       234956, 349644, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451656, 3, 234956,
                                                                       234977, 349672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451692, 3, 234977,
                                                                       234998, 349700, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451728, 3, 234998,
                                                                       235019, 349728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451764, 3, 235019,
                                                                       235040, 349756, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451800, 3, 235040,
                                                                       235061, 349784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451836, 3, 235103,
                                                                       235124, 349812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451872, 3, 235124,
                                                                       235145, 349840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451908, 3, 235145,
                                                                       235166, 349868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451944, 3, 235166,
                                                                       235187, 349896, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 451980, 3, 235187,
                                                                       235208, 349924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 452016, 3, 235208,
                                                                       235229, 349952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 452052, 3, 235229,
                                                                       235250, 349980, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 452088, 3, 235250,
                                                                       235271, 350008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 452124, 3, 235271,
                                                                       235292, 350036, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 452160, 3, 235292,
                                                                       235313, 350064, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 452196, 3, 235313,
                                                                       235334, 350092, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 452232, 3, 235334,
                                                                       235355, 350120, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 452268, 3, 235355,
                                                                       235376, 350148, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 452304, 0, 3,
                                                                       451368, 349448, 451404,
                                                                       235418, 235481, 350176,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 452412, 0, 3,
                                                                       451404, 349476, 451440,
                                                                       235481, 235544, 350260,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 452520, 0, 3,
                                                                       451440, 349504, 451476,
                                                                       235544, 235607, 350344,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 452628, 0, 3,
                                                                       451476, 349532, 451512,
                                                                       235607, 235670, 350428,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 452736, 0, 3,
                                                                       451512, 349560, 451548,
                                                                       235670, 235733, 350512,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 452844, 0, 3,
                                                                       451548, 349588, 451584,
                                                                       235733, 235796, 350596,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 452952, 0, 3,
                                                                       451584, 349616, 451620,
                                                                       235796, 235859, 350680,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453060, 0, 3,
                                                                       451620, 349644, 451656,
                                                                       235859, 235922, 350764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453168, 0, 3,
                                                                       451656, 349672, 451692,
                                                                       235922, 235985, 350848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453276, 0, 3,
                                                                       451692, 349700, 451728,
                                                                       235985, 236048, 350932,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453384, 0, 3,
                                                                       451728, 349728, 451764,
                                                                       236048, 236111, 351016,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453492, 0, 3,
                                                                       451764, 349756, 451800,
                                                                       236111, 236174, 351100,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453600, 0, 3,
                                                                       451836, 349812, 451872,
                                                                       236300, 236363, 351184,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453708, 0, 3,
                                                                       451872, 349840, 451908,
                                                                       236363, 236426, 351268,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453816, 0, 3,
                                                                       451908, 349868, 451944,
                                                                       236426, 236489, 351352,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 453924, 0, 3,
                                                                       451944, 349896, 451980,
                                                                       236489, 236552, 351436,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 454032, 0, 3,
                                                                       451980, 349924, 452016,
                                                                       236552, 236615, 351520,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 454140, 0, 3,
                                                                       452016, 349952, 452052,
                                                                       236615, 236678, 351604,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 454248, 0, 3,
                                                                       452052, 349980, 452088,
                                                                       236678, 236741, 351688,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 454356, 0, 3,
                                                                       452088, 350008, 452124,
                                                                       236741, 236804, 351772,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 454464, 0, 3,
                                                                       452124, 350036, 452160,
                                                                       236804, 236867, 351856,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 454572, 0, 3,
                                                                       452160, 350064, 452196,
                                                                       236867, 236930, 351940,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 454680, 0, 3,
                                                                       452196, 350092, 452232,
                                                                       236930, 236993, 352024,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 454788, 0, 3,
                                                                       452232, 350120, 452268,
                                                                       236993, 237056, 352108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 454896, 0, 3,
                                                                       452304, 350176, 452412,
                                                                       237182, 237308, 352192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 455112, 0, 3,
                                                                       452412, 350260, 452520,
                                                                       237308, 237434, 352360,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 455328, 0, 3,
                                                                       452520, 350344, 452628,
                                                                       237434, 237560, 352528,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 455544, 0, 3,
                                                                       452628, 350428, 452736,
                                                                       237560, 237686, 352696,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 455760, 0, 3,
                                                                       452736, 350512, 452844,
                                                                       237686, 237812, 352864,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 455976, 0, 3,
                                                                       452844, 350596, 452952,
                                                                       237812, 237938, 353032,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 456192, 0, 3,
                                                                       452952, 350680, 453060,
                                                                       237938, 238064, 353200,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 456408, 0, 3,
                                                                       453060, 350764, 453168,
                                                                       238064, 238190, 353368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 456624, 0, 3,
                                                                       453168, 350848, 453276,
                                                                       238190, 238316, 353536,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 456840, 0, 3,
                                                                       453276, 350932, 453384,
                                                                       238316, 238442, 353704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 457056, 0, 3,
                                                                       453384, 351016, 453492,
                                                                       238442, 238568, 353872,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 457272, 0, 3,
                                                                       453600, 351184, 453708,
                                                                       238820, 238946, 354040,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 457488, 0, 3,
                                                                       453708, 351268, 453816,
                                                                       238946, 239072, 354208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 457704, 0, 3,
                                                                       453816, 351352, 453924,
                                                                       239072, 239198, 354376,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 457920, 0, 3,
                                                                       453924, 351436, 454032,
                                                                       239198, 239324, 354544,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 458136, 0, 3,
                                                                       454032, 351520, 454140,
                                                                       239324, 239450, 354712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 458352, 0, 3,
                                                                       454140, 351604, 454248,
                                                                       239450, 239576, 354880,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 458568, 0, 3,
                                                                       454248, 351688, 454356,
                                                                       239576, 239702, 355048,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 458784, 0, 3,
                                                                       454356, 351772, 454464,
                                                                       239702, 239828, 355216,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 459000, 0, 3,
                                                                       454464, 351856, 454572,
                                                                       239828, 239954, 355384,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 459216, 0, 3,
                                                                       454572, 351940, 454680,
                                                                       239954, 240080, 355552,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 459432, 0, 3,
                                                                       454680, 352024, 454788,
                                                                       240080, 240206, 355720,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 459648, 0, 3,
                                                                       454896, 352192, 455112,
                                                                       240458, 240668, 355888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 460008, 0, 3,
                                                                       455112, 352360, 455328,
                                                                       240668, 240878, 356168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 460368, 0, 3,
                                                                       455328, 352528, 455544,
                                                                       240878, 241088, 356448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 460728, 0, 3,
                                                                       455544, 352696, 455760,
                                                                       241088, 241298, 356728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 461088, 0, 3,
                                                                       455760, 352864, 455976,
                                                                       241298, 241508, 357008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 461448, 0, 3,
                                                                       455976, 353032, 456192,
                                                                       241508, 241718, 357288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 461808, 0, 3,
                                                                       456192, 353200, 456408,
                                                                       241718, 241928, 357568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 462168, 0, 3,
                                                                       456408, 353368, 456624,
                                                                       241928, 242138, 357848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 462528, 0, 3,
                                                                       456624, 353536, 456840,
                                                                       242138, 242348, 358128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 462888, 0, 3,
                                                                       456840, 353704, 457056,
                                                                       242348, 242558, 358408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 463248, 0, 3,
                                                                       457272, 354040, 457488,
                                                                       242978, 243188, 358688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 463608, 0, 3,
                                                                       457488, 354208, 457704,
                                                                       243188, 243398, 358968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 463968, 0, 3,
                                                                       457704, 354376, 457920,
                                                                       243398, 243608, 359248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 464328, 0, 3,
                                                                       457920, 354544, 458136,
                                                                       243608, 243818, 359528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 464688, 0, 3,
                                                                       458136, 354712, 458352,
                                                                       243818, 244028, 359808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 465048, 0, 3,
                                                                       458352, 354880, 458568,
                                                                       244028, 244238, 360088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 465408, 0, 3,
                                                                       458568, 355048, 458784,
                                                                       244238, 244448, 360368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 465768, 0, 3,
                                                                       458784, 355216, 459000,
                                                                       244448, 244658, 360648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 466128, 0, 3,
                                                                       459000, 355384, 459216,
                                                                       244658, 244868, 360928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 466488, 0, 3,
                                                                       459216, 355552, 459432,
                                                                       244868, 245078, 361208,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 466848, 0, 3,
                                                                       459648, 355888, 460008,
                                                                       245498, 245813, 361488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 467388, 0, 3,
                                                                       460008, 356168, 460368,
                                                                       245813, 246128, 361908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 467928, 0, 3,
                                                                       460368, 356448, 460728,
                                                                       246128, 246443, 362328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 468468, 0, 3,
                                                                       460728, 356728, 461088,
                                                                       246443, 246758, 362748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 469008, 0, 3,
                                                                       461088, 357008, 461448,
                                                                       246758, 247073, 363168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 469548, 0, 3,
                                                                       461448, 357288, 461808,
                                                                       247073, 247388, 363588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 470088, 0, 3,
                                                                       461808, 357568, 462168,
                                                                       247388, 247703, 364008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 470628, 0, 3,
                                                                       462168, 357848, 462528,
                                                                       247703, 248018, 364428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 471168, 0, 3,
                                                                       462528, 358128, 462888,
                                                                       248018, 248333, 364848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 471708, 0, 3,
                                                                       463248, 358688, 463608,
                                                                       248963, 249278, 365268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 472248, 0, 3,
                                                                       463608, 358968, 463968,
                                                                       249278, 249593, 365688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 472788, 0, 3,
                                                                       463968, 359248, 464328,
                                                                       249593, 249908, 366108,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 473328, 0, 3,
                                                                       464328, 359528, 464688,
                                                                       249908, 250223, 366528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 473868, 0, 3,
                                                                       464688, 359808, 465048,
                                                                       250223, 250538, 366948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 474408, 0, 3,
                                                                       465048, 360088, 465408,
                                                                       250538, 250853, 367368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 474948, 0, 3,
                                                                       465408, 360368, 465768,
                                                                       250853, 251168, 367788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 475488, 0, 3,
                                                                       465768, 360648, 466128,
                                                                       251168, 251483, 368208,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 476028, 0, 3,
                                                                       466128, 360928, 466488,
                                                                       251483, 251798, 368628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 476568, 0, 3,
                                                                       466848, 361488, 467388,
                                                                       252428, 252869, 369048,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 477324, 0, 3,
                                                                       467388, 361908, 467928,
                                                                       252869, 253310, 369636,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 478080, 0, 3,
                                                                       467928, 362328, 468468,
                                                                       253310, 253751, 370224,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 478836, 0, 3,
                                                                       468468, 362748, 469008,
                                                                       253751, 254192, 370812,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 479592, 0, 3,
                                                                       469008, 363168, 469548,
                                                                       254192, 254633, 371400,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 480348, 0, 3,
                                                                       469548, 363588, 470088,
                                                                       254633, 255074, 371988,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 481104, 0, 3,
                                                                       470088, 364008, 470628,
                                                                       255074, 255515, 372576,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 481860, 0, 3,
                                                                       470628, 364428, 471168,
                                                                       255515, 255956, 373164,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 482616, 0, 3,
                                                                       471708, 365268, 472248,
                                                                       256838, 257279, 373752,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 483372, 0, 3,
                                                                       472248, 365688, 472788,
                                                                       257279, 257720, 374340,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 484128, 0, 3,
                                                                       472788, 366108, 473328,
                                                                       257720, 258161, 374928,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 484884, 0, 3,
                                                                       473328, 366528, 473868,
                                                                       258161, 258602, 375516,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 485640, 0, 3,
                                                                       473868, 366948, 474408,
                                                                       258602, 259043, 376104,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 486396, 0, 3,
                                                                       474408, 367368, 474948,
                                                                       259043, 259484, 376692,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 487152, 0, 3,
                                                                       474948, 367788, 475488,
                                                                       259484, 259925, 377280,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 487908, 0, 3,
                                                                       475488, 368208, 476028,
                                                                       259925, 260366, 377868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 488664, 0, 3,
                                                                       476568, 369048, 477324,
                                                                       261248, 261836, 378456,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 489672, 0, 3,
                                                                       477324, 369636, 478080,
                                                                       261836, 262424, 379240,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 490680, 0, 3,
                                                                       478080, 370224, 478836,
                                                                       262424, 263012, 380024,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 491688, 0, 3,
                                                                       478836, 370812, 479592,
                                                                       263012, 263600, 380808,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 492696, 0, 3,
                                                                       479592, 371400, 480348,
                                                                       263600, 264188, 381592,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 493704, 0, 3,
                                                                       480348, 371988, 481104,
                                                                       264188, 264776, 382376,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 494712, 0, 3,
                                                                       481104, 372576, 481860,
                                                                       264776, 265364, 383160,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 495720, 0, 3,
                                                                       482616, 373752, 483372,
                                                                       266540, 267128, 383944,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 496728, 0, 3,
                                                                       483372, 374340, 484128,
                                                                       267128, 267716, 384728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 497736, 0, 3,
                                                                       484128, 374928, 484884,
                                                                       267716, 268304, 385512,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 498744, 0, 3,
                                                                       484884, 375516, 485640,
                                                                       268304, 268892, 386296,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 499752, 0, 3,
                                                                       485640, 376104, 486396,
                                                                       268892, 269480, 387080,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 500760, 0, 3,
                                                                       486396, 376692, 487152,
                                                                       269480, 270068, 387864,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 501768, 0, 3,
                                                                       487152, 377280, 487908,
                                                                       270068, 270656, 388648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 502776, 0, 3,
                                                                       488664, 378456, 489672,
                                                                       271832, 272588, 389432,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 504072, 0, 3,
                                                                       489672, 379240, 490680,
                                                                       272588, 273344, 390440,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 505368, 0, 3,
                                                                       490680, 380024, 491688,
                                                                       273344, 274100, 391448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 506664, 0, 3,
                                                                       491688, 380808, 492696,
                                                                       274100, 274856, 392456,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 507960, 0, 3,
                                                                       492696, 381592, 493704,
                                                                       274856, 275612, 393464,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 509256, 0, 3,
                                                                       493704, 382376, 494712,
                                                                       275612, 276368, 394472,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 510552, 0, 3,
                                                                       495720, 383944, 496728,
                                                                       277880, 278636, 395480,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 511848, 0, 3,
                                                                       496728, 384728, 497736,
                                                                       278636, 279392, 396488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 513144, 0, 3,
                                                                       497736, 385512, 498744,
                                                                       279392, 280148, 397496,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 514440, 0, 3,
                                                                       498744, 386296, 499752,
                                                                       280148, 280904, 398504,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 515736, 0, 3,
                                                                       499752, 387080, 500760,
                                                                       280904, 281660, 399512,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 517032, 0, 3,
                                                                       500760, 387864, 501768,
                                                                       281660, 282416, 400520,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 518328, 0, 3,
                                                                       502776, 389432, 504072,
                                                                       283928, 284873, 401528,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 519948, 0, 3,
                                                                       504072, 390440, 505368,
                                                                       284873, 285818, 402788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 521568, 0, 3,
                                                                       505368, 391448, 506664,
                                                                       285818, 286763, 404048,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 523188, 0, 3,
                                                                       506664, 392456, 507960,
                                                                       286763, 287708, 405308,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 524808, 0, 3,
                                                                       507960, 393464, 509256,
                                                                       287708, 288653, 406568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 526428, 0, 3,
                                                                       510552, 395480, 511848,
                                                                       290543, 291488, 407828,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 528048, 0, 3,
                                                                       511848, 396488, 513144,
                                                                       291488, 292433, 409088,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 529668, 0, 3,
                                                                       513144, 397496, 514440,
                                                                       292433, 293378, 410348,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 531288, 0, 3,
                                                                       514440, 398504, 515736,
                                                                       293378, 294323, 411608,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 532908, 0, 3,
                                                                       515736, 399512, 517032,
                                                                       294323, 295268, 412868,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 534528, 0, 3,
                                                                       518328, 401528, 519948,
                                                                       297158, 298313, 414128,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 536508, 0, 3,
                                                                       519948, 402788, 521568,
                                                                       298313, 299468, 415668,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 538488, 0, 3,
                                                                       521568, 404048, 523188,
                                                                       299468, 300623, 417208,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 540468, 0, 3,
                                                                       523188, 405308, 524808,
                                                                       300623, 301778, 418748,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 542448, 0, 3,
                                                                       526428, 407828, 528048,
                                                                       304088, 305243, 420288,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 544428, 0, 3,
                                                                       528048, 409088, 529668,
                                                                       305243, 306398, 421828,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 546408, 0, 3,
                                                                       529668, 410348, 531288,
                                                                       306398, 307553, 423368,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 548388, 0, 3,
                                                                       531288, 411608, 532908,
                                                                       307553, 308708, 424908,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 550368, 0, 3,
                                                                       534528, 414128, 536508,
                                                                       311018, 312404, 426448,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 552744, 0, 3,
                                                                       536508, 415668, 538488,
                                                                       312404, 313790, 428296,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 555120, 0, 3,
                                                                       538488, 417208, 540468,
                                                                       313790, 315176, 430144,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 557496, 0, 3,
                                                                       542448, 420288, 544428,
                                                                       317948, 319334, 431992,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 559872, 0, 3,
                                                                       544428, 421828, 546408,
                                                                       319334, 320720, 433840,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 562248, 0, 3,
                                                                       546408, 423368, 548388,
                                                                       320720, 322106, 435688,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 564624, 0, 3,
                                                                       550368, 426448, 552744,
                                                                       324878, 326516, 437536,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 567432, 0, 3,
                                                                       552744, 428296, 555120,
                                                                       326516, 328154, 439720,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 570240, 0, 3,
                                                                       557496, 431992, 559872,
                                                                       331430, 333068, 441904,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 573048, 0, 3,
                                                                       559872, 433840, 562248,
                                                                       333068, 334706, 444088,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsk_three_center_electron_repulsion_0(buffer, 575856, 0, 3,
                                                                       564624, 437536, 567432,
                                                                       337982, 339893, 446272,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsk_three_center_electron_repulsion_0(buffer, 579132, 0, 3,
                                                                       570240, 441904, 573048,
                                                                       343715, 345626, 448820,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 582408, 488664, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 583836, 495720, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 585264, 502776, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 587100, 510552, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 588936, 518328, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 591231, 526428, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 593526, 534528, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 596331, 542448, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 599136, 550368, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 602502, 557496, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 605868, 564624, 2808, ncols);

                    simdfunc::contract_primitives(buffer, 609846, 570240, 2808, ncols);

                    simdfunc::contract_primitives(buffer, 613824, 575856, 3276, ncols);

                    simdfunc::contract_primitives(buffer, 618465, 579132, 3276, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 583416, 582408, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 584844, 583836, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 586560, 585264, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 588396, 587100, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 590556, 588936, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 592851, 591231, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 595506, 593526, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 598311, 596331, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 601512, 599136, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 604878, 602502, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 608676, 605868, 78, 1, nmax);

        simdtrf::transform_k_inner(buffer, 612654, 609846, 78, 1, nmax);

        simdtrf::transform_k_inner(buffer, 617100, 613824, 91, 1, nmax);

        simdtrf::transform_k_inner(buffer, 621741, 618465, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 623106, 583416, 586560, 15,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 624366, 584844, 588396, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 625626, 586560, 590556, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 627246, 588396, 592851, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 628866, 590556, 595506, 15,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 630891, 592851, 598311, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 632916, 595506, 601512, 15,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 635391, 598311, 604878, 15,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 637866, 601512, 608676, 15,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 640836, 604878, 612654, 15,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 643806, 608676, 617100, 15,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 647316, 612654, 621741, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 650826, 623106, 625626, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 653346, 624366, 627246, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 655866, 625626, 628866, 15,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 659106, 627246, 630891, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 662346, 628866, 632916, 15,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 666396, 630891, 635391, 15,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 670446, 632916, 637866, 15,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 675396, 635391, 640836, 15,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 680346, 637866, 643806, 15,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 686286, 640836, 647316, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 692226, 650826, 655866, 15,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 696426, 653346, 659106, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 700626, 655866, 662346, 15,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 706026, 659106, 666396, 15,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 711426, 662346, 670446, 15,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 718176, 666396, 675396, 15,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 724926, 670446, 680346, 15,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 733176, 675396, 686286, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 741426, 692226, 700626, 15,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 747726, 696426, 706026, 15,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 754026, 700626, 711426, 15,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 762126, 706026, 718176, 15,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 770226, 711426, 724926, 15,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 780351, 718176, 733176, 15,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 790476, 741426, 754026, 15,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 799296, 747726, 762126, 15,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 808116, 754026, 770226, 15,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 819456, 762126, 780351, 15,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 830796, 790476, 808116, 15,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 842556, 799296, 819456, 15,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 854316, 842556, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 854316, 195, nmax);

        simdtrf::transform_i_inner(buffer, 854316, 830796, 28, 15, nmax);

        simdtrf::transform_i_outer(values + 2535 * nvalues + n * npairs, nvalues, buffer, 854316,
                                   195, nmax);
    }

    for (size_t m = 0; m < 5070; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
