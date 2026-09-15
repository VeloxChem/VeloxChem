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


#include "SimdThreeCenterElectronRepulsionRecIIL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
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
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iil_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iil_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 580553, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2873 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 580553, 418608, 23191, dimensions);

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

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 20,
                                                             ncols, fj, 6, fq);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 8, 9,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 9, 10,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 10, 11,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 11, 12,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 12, 13,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 13, 14,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 125, 0, 3, 14, 15,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 131, 0, 3, 15, 16,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 137, 0, 3, 16, 17,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 17, 18,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 149, 0, 3, 18, 19,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 155, 0, 3, 19, 20,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 161, 0, 3, 20, 21,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 167, 0, 3, 21, 22,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 22, 23,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 179, 0, 3, 23, 24,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 185, 0, 3, 24, 25,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 191, 0, 3, 25, 26,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 197, 0, 3, 26, 27,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 29, 32,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 32, 35,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 35, 38,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 38, 41,
                                                                       107, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 41, 44,
                                                                       113, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 44, 47,
                                                                       119, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 47, 50,
                                                                       125, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 50, 53,
                                                                       131, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 53, 56,
                                                                       137, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 56, 59,
                                                                       143, 149, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 303, 0, 3, 59, 62,
                                                                       149, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 62, 65,
                                                                       155, 161, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 323, 0, 3, 65, 68,
                                                                       161, 167, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 333, 0, 3, 68, 71,
                                                                       167, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 343, 0, 3, 71, 74,
                                                                       173, 179, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 353, 0, 3, 74, 77,
                                                                       179, 185, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 363, 0, 3, 77, 80,
                                                                       185, 191, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 373, 0, 3, 80, 83,
                                                                       191, 197, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 383, 0, 3, 89, 95,
                                                                       203, 213, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 95,
                                                                       101, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 101,
                                                                       107, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 107,
                                                                       113, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 443, 0, 3, 113,
                                                                       119, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 119,
                                                                       125, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 473, 0, 3, 125,
                                                                       131, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 131,
                                                                       137, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 503, 0, 3, 137,
                                                                       143, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 143,
                                                                       149, 293, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 533, 0, 3, 149,
                                                                       155, 303, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 155,
                                                                       161, 313, 323, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 563, 0, 3, 161,
                                                                       167, 323, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 167,
                                                                       173, 333, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 593, 0, 3, 173,
                                                                       179, 343, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 608, 0, 3, 179,
                                                                       185, 353, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 623, 0, 3, 185,
                                                                       191, 363, 373, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 203,
                                                                       213, 383, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 659, 0, 3, 213,
                                                                       223, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 680, 0, 3, 223,
                                                                       233, 413, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 701, 0, 3, 233,
                                                                       243, 428, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 243,
                                                                       253, 443, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 253,
                                                                       263, 458, 473, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 764, 0, 3, 263,
                                                                       273, 473, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 785, 0, 3, 273,
                                                                       283, 488, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 283,
                                                                       293, 503, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 827, 0, 3, 293,
                                                                       303, 518, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 848, 0, 3, 303,
                                                                       313, 533, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 869, 0, 3, 313,
                                                                       323, 548, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 323,
                                                                       333, 563, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 911, 0, 3, 333,
                                                                       343, 578, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 932, 0, 3, 343,
                                                                       353, 593, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 953, 0, 3, 353,
                                                                       363, 608, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 383,
                                                                       398, 638, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 398,
                                                                       413, 659, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 413,
                                                                       428, 680, 701, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 428,
                                                                       443, 701, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 443,
                                                                       458, 722, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 458,
                                                                       473, 743, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 473,
                                                                       488, 764, 785, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 488,
                                                                       503, 785, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1198, 0, 3, 503,
                                                                       518, 806, 827, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 518,
                                                                       533, 827, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 533,
                                                                       548, 848, 869, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 548,
                                                                       563, 869, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 563,
                                                                       578, 890, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1338, 0, 3, 578,
                                                                       593, 911, 932, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1366, 0, 3, 593,
                                                                       608, 932, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1394, 0, 3, 638,
                                                                       659, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1430, 0, 3, 659,
                                                                       680, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1466, 0, 3, 680,
                                                                       701, 1030, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1502, 0, 3, 701,
                                                                       722, 1058, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1538, 0, 3, 722,
                                                                       743, 1086, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1574, 0, 3, 743,
                                                                       764, 1114, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1610, 0, 3, 764,
                                                                       785, 1142, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1646, 0, 3, 785,
                                                                       806, 1170, 1198, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1682, 0, 3, 806,
                                                                       827, 1198, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1718, 0, 3, 827,
                                                                       848, 1226, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1754, 0, 3, 848,
                                                                       869, 1254, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1790, 0, 3, 869,
                                                                       890, 1282, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1826, 0, 3, 890,
                                                                       911, 1310, 1338, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1862, 0, 3, 911,
                                                                       932, 1338, 1366, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 974,
                                                                       1002, 1394, 1430, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1943, 0, 3, 1002,
                                                                       1030, 1430, 1466, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1030,
                                                                       1058, 1466, 1502, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2033, 0, 3, 1058,
                                                                       1086, 1502, 1538, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2078, 0, 3, 1086,
                                                                       1114, 1538, 1574, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2123, 0, 3, 1114,
                                                                       1142, 1574, 1610, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 1142,
                                                                       1170, 1610, 1646, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2213, 0, 3, 1170,
                                                                       1198, 1646, 1682, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2258, 0, 3, 1198,
                                                                       1226, 1682, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2303, 0, 3, 1226,
                                                                       1254, 1718, 1754, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1254,
                                                                       1282, 1754, 1790, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2393, 0, 3, 1282,
                                                                       1310, 1790, 1826, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2438, 0, 3, 1310,
                                                                       1338, 1826, 1862, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1394,
                                                                       1430, 1898, 1943, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2538, 0, 3, 1430,
                                                                       1466, 1943, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2593, 0, 3, 1466,
                                                                       1502, 1988, 2033, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1502,
                                                                       1538, 2033, 2078, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2703, 0, 3, 1538,
                                                                       1574, 2078, 2123, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2758, 0, 3, 1574,
                                                                       1610, 2123, 2168, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1610,
                                                                       1646, 2168, 2213, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2868, 0, 3, 1646,
                                                                       1682, 2213, 2258, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2923, 0, 3, 1682,
                                                                       1718, 2258, 2303, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2978, 0, 3, 1718,
                                                                       1754, 2303, 2348, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3033, 0, 3, 1754,
                                                                       1790, 2348, 2393, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 1790,
                                                                       1826, 2393, 2438, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 1898,
                                                                       1943, 2483, 2538, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3209, 0, 3, 1943,
                                                                       1988, 2538, 2593, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3275, 0, 3, 1988,
                                                                       2033, 2593, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3341, 0, 3, 2033,
                                                                       2078, 2648, 2703, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3407, 0, 3, 2078,
                                                                       2123, 2703, 2758, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 2123,
                                                                       2168, 2758, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3539, 0, 3, 2168,
                                                                       2213, 2813, 2868, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3605, 0, 3, 2213,
                                                                       2258, 2868, 2923, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3671, 0, 3, 2258,
                                                                       2303, 2923, 2978, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3737, 0, 3, 2303,
                                                                       2348, 2978, 3033, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 3803, 0, 3, 2348,
                                                                       2393, 3033, 3088, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3869, 0, 3, 2483,
                                                                       2538, 3143, 3209, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 3947, 0, 3, 2538,
                                                                       2593, 3209, 3275, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4025, 0, 3, 2593,
                                                                       2648, 3275, 3341, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4103, 0, 3, 2648,
                                                                       2703, 3341, 3407, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4181, 0, 3, 2703,
                                                                       2758, 3407, 3473, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4259, 0, 3, 2758,
                                                                       2813, 3473, 3539, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4337, 0, 3, 2813,
                                                                       2868, 3539, 3605, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4415, 0, 3, 2868,
                                                                       2923, 3605, 3671, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4493, 0, 3, 2923,
                                                                       2978, 3671, 3737, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 4571, 0, 3, 2978,
                                                                       3033, 3737, 3803, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4649, 0, 3, 3143,
                                                                       3209, 3869, 3947, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4740, 0, 3, 3209,
                                                                       3275, 3947, 4025, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4831, 0, 3, 3275,
                                                                       3341, 4025, 4103, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 4922, 0, 3, 3341,
                                                                       3407, 4103, 4181, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5013, 0, 3, 3407,
                                                                       3473, 4181, 4259, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5104, 0, 3, 3473,
                                                                       3539, 4259, 4337, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5195, 0, 3, 3539,
                                                                       3605, 4337, 4415, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5286, 0, 3, 3605,
                                                                       3671, 4415, 4493, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 5377, 0, 3, 3671,
                                                                       3737, 4493, 4571, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5468, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5471, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5474, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5477, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5480, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5483, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5486, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5489, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5492, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5495, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5498, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5501, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5504, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5507, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5510, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5513, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5516, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5519, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5522, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5525, 3, 10, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5534, 3, 11, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5543, 3, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5552, 3, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5561, 3, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5570, 3, 15, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5579, 3, 16, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5588, 3, 17, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5597, 3, 18, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5606, 3, 19, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5615, 3, 20, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5624, 3, 21, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5633, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5642, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5651, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5660, 3, 25, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5669, 3, 26, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 5678, 3, 27, 86,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5687, 3, 35, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5705, 3, 38, 107,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5723, 3, 41, 113,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5741, 3, 44, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5759, 3, 47, 125,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5777, 3, 50, 131,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5795, 3, 53, 137,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5813, 3, 56, 143,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5831, 3, 59, 149,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5849, 3, 62, 155,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5867, 3, 65, 161,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5885, 3, 68, 167,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5903, 3, 71, 173,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5921, 3, 74, 179,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5939, 3, 77, 185,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5957, 3, 80, 191,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 5975, 3, 83, 197,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 5993, 3, 101, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6023, 3, 107, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6053, 3, 113, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6083, 3, 119, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6113, 3, 125, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6143, 3, 131, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6173, 3, 137, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6203, 3, 143, 293,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6233, 3, 149, 303,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6263, 3, 155, 313,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6293, 3, 161, 323,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6323, 3, 167, 333,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6353, 3, 173, 343,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6383, 3, 179, 353,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6413, 3, 185, 363,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 6443, 3, 191, 373,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6473, 3, 223, 413,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6518, 3, 233, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6563, 3, 243, 443,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6608, 3, 253, 458,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6653, 3, 263, 473,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6698, 3, 273, 488,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6743, 3, 283, 503,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6788, 3, 293, 518,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6833, 3, 303, 533,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6878, 3, 313, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6923, 3, 323, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6968, 3, 333, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7013, 3, 343, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7058, 3, 353, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 7103, 3, 363, 623,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7148, 3, 413, 680,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7211, 3, 428, 701,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7274, 3, 443, 722,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7337, 3, 458, 743,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7400, 3, 473, 764,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7463, 3, 488, 785,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7526, 3, 503, 806,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7589, 3, 518, 827,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7652, 3, 533, 848,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7715, 3, 548, 869,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7778, 3, 563, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7841, 3, 578, 911,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7904, 3, 593, 932,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 7967, 3, 608, 953,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8030, 3, 680,
                                                                       1030, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8114, 3, 701,
                                                                       1058, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8198, 3, 722,
                                                                       1086, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8282, 3, 743,
                                                                       1114, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8366, 3, 764,
                                                                       1142, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8450, 3, 785,
                                                                       1170, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8534, 3, 806,
                                                                       1198, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8618, 3, 827,
                                                                       1226, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8702, 3, 848,
                                                                       1254, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8786, 3, 869,
                                                                       1282, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8870, 3, 890,
                                                                       1310, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 8954, 3, 911,
                                                                       1338, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 9038, 3, 932,
                                                                       1366, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9122, 3, 1030,
                                                                       1466, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9230, 3, 1058,
                                                                       1502, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9338, 3, 1086,
                                                                       1538, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9446, 3, 1114,
                                                                       1574, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9554, 3, 1142,
                                                                       1610, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9662, 3, 1170,
                                                                       1646, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9770, 3, 1198,
                                                                       1682, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9878, 3, 1226,
                                                                       1718, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 9986, 3, 1254,
                                                                       1754, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10094, 3, 1282,
                                                                       1790, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10202, 3, 1310,
                                                                       1826, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 10310, 3, 1338,
                                                                       1862, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10418, 3, 1466,
                                                                       1988, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10553, 3, 1502,
                                                                       2033, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10688, 3, 1538,
                                                                       2078, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10823, 3, 1574,
                                                                       2123, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 10958, 3, 1610,
                                                                       2168, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11093, 3, 1646,
                                                                       2213, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11228, 3, 1682,
                                                                       2258, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11363, 3, 1718,
                                                                       2303, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11498, 3, 1754,
                                                                       2348, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11633, 3, 1790,
                                                                       2393, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 11768, 3, 1826,
                                                                       2438, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 11903, 3, 1988,
                                                                       2593, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12068, 3, 2033,
                                                                       2648, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12233, 3, 2078,
                                                                       2703, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12398, 3, 2123,
                                                                       2758, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12563, 3, 2168,
                                                                       2813, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12728, 3, 2213,
                                                                       2868, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 12893, 3, 2258,
                                                                       2923, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13058, 3, 2303,
                                                                       2978, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13223, 3, 2348,
                                                                       3033, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 13388, 3, 2393,
                                                                       3088, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13553, 3, 2593,
                                                                       3275, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13751, 3, 2648,
                                                                       3341, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 13949, 3, 2703,
                                                                       3407, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14147, 3, 2758,
                                                                       3473, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14345, 3, 2813,
                                                                       3539, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14543, 3, 2868,
                                                                       3605, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14741, 3, 2923,
                                                                       3671, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 14939, 3, 2978,
                                                                       3737, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 15137, 3, 3033,
                                                                       3803, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15335, 3, 3275,
                                                                       4025, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15569, 3, 3341,
                                                                       4103, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 15803, 3, 3407,
                                                                       4181, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16037, 3, 3473,
                                                                       4259, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16271, 3, 3539,
                                                                       4337, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16505, 3, 3605,
                                                                       4415, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16739, 3, 3671,
                                                                       4493, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 16973, 3, 3737,
                                                                       4571, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17207, 3, 4025,
                                                                       4831, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17480, 3, 4103,
                                                                       4922, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 17753, 3, 4181,
                                                                       5013, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 18026, 3, 4259,
                                                                       5104, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 18299, 3, 4337,
                                                                       5195, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 18572, 3, 4415,
                                                                       5286, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 18845, 3, 4493,
                                                                       5377, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19118, 3, 8, 9,
                                                                       5468, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19124, 3, 9, 10,
                                                                       5471, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19130, 3, 10, 11,
                                                                       5474, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19136, 3, 11, 12,
                                                                       5477, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19142, 3, 12, 13,
                                                                       5480, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19148, 3, 13, 14,
                                                                       5483, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19154, 3, 14, 15,
                                                                       5486, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19160, 3, 15, 16,
                                                                       5489, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19166, 3, 16, 17,
                                                                       5492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19172, 3, 17, 18,
                                                                       5495, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19178, 3, 18, 19,
                                                                       5498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19184, 3, 19, 20,
                                                                       5501, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19190, 3, 20, 21,
                                                                       5504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19196, 3, 21, 22,
                                                                       5507, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19202, 3, 22, 23,
                                                                       5510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19208, 3, 23, 24,
                                                                       5513, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19214, 3, 24, 25,
                                                                       5516, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19220, 3, 25, 26,
                                                                       5519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 19226, 3, 26, 27,
                                                                       5522, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19232, 0, 3,
                                                                       19118, 5468, 19124, 5525,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19250, 0, 3,
                                                                       19124, 5471, 19130, 5534,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19268, 0, 3,
                                                                       19130, 5474, 19136, 5543,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19286, 0, 3,
                                                                       19136, 5477, 19142, 5552,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19304, 0, 3,
                                                                       19142, 5480, 19148, 5561,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19322, 0, 3,
                                                                       19148, 5483, 19154, 5570,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19340, 0, 3,
                                                                       19154, 5486, 19160, 5579,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19358, 0, 3,
                                                                       19160, 5489, 19166, 5588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19376, 0, 3,
                                                                       19166, 5492, 19172, 5597,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19394, 0, 3,
                                                                       19172, 5495, 19178, 5606,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19412, 0, 3,
                                                                       19178, 5498, 19184, 5615,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19430, 0, 3,
                                                                       19184, 5501, 19190, 5624,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19448, 0, 3,
                                                                       19190, 5504, 19196, 5633,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19466, 0, 3,
                                                                       19196, 5507, 19202, 5642,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19484, 0, 3,
                                                                       19202, 5510, 19208, 5651,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19502, 0, 3,
                                                                       19208, 5513, 19214, 5660,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19520, 0, 3,
                                                                       19214, 5516, 19220, 5669,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 19538, 0, 3,
                                                                       19220, 5519, 19226, 5678,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19556, 0, 3,
                                                                       19232, 5525, 19250, 89,
                                                                       95, 5687, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19592, 0, 3,
                                                                       19250, 5534, 19268, 95,
                                                                       101, 5705, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19628, 0, 3,
                                                                       19268, 5543, 19286, 101,
                                                                       107, 5723, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19664, 0, 3,
                                                                       19286, 5552, 19304, 107,
                                                                       113, 5741, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19700, 0, 3,
                                                                       19304, 5561, 19322, 113,
                                                                       119, 5759, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19736, 0, 3,
                                                                       19322, 5570, 19340, 119,
                                                                       125, 5777, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19772, 0, 3,
                                                                       19340, 5579, 19358, 125,
                                                                       131, 5795, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19808, 0, 3,
                                                                       19358, 5588, 19376, 131,
                                                                       137, 5813, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19844, 0, 3,
                                                                       19376, 5597, 19394, 137,
                                                                       143, 5831, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19880, 0, 3,
                                                                       19394, 5606, 19412, 143,
                                                                       149, 5849, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19916, 0, 3,
                                                                       19412, 5615, 19430, 149,
                                                                       155, 5867, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19952, 0, 3,
                                                                       19430, 5624, 19448, 155,
                                                                       161, 5885, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 19988, 0, 3,
                                                                       19448, 5633, 19466, 161,
                                                                       167, 5903, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 20024, 0, 3,
                                                                       19466, 5642, 19484, 167,
                                                                       173, 5921, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 20060, 0, 3,
                                                                       19484, 5651, 19502, 173,
                                                                       179, 5939, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 20096, 0, 3,
                                                                       19502, 5660, 19520, 179,
                                                                       185, 5957, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 20132, 0, 3,
                                                                       19520, 5669, 19538, 185,
                                                                       191, 5975, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20168, 0, 3,
                                                                       19556, 5687, 19592, 203,
                                                                       213, 5993, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20228, 0, 3,
                                                                       19592, 5705, 19628, 213,
                                                                       223, 6023, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20288, 0, 3,
                                                                       19628, 5723, 19664, 223,
                                                                       233, 6053, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20348, 0, 3,
                                                                       19664, 5741, 19700, 233,
                                                                       243, 6083, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20408, 0, 3,
                                                                       19700, 5759, 19736, 243,
                                                                       253, 6113, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20468, 0, 3,
                                                                       19736, 5777, 19772, 253,
                                                                       263, 6143, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20528, 0, 3,
                                                                       19772, 5795, 19808, 263,
                                                                       273, 6173, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20588, 0, 3,
                                                                       19808, 5813, 19844, 273,
                                                                       283, 6203, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20648, 0, 3,
                                                                       19844, 5831, 19880, 283,
                                                                       293, 6233, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20708, 0, 3,
                                                                       19880, 5849, 19916, 293,
                                                                       303, 6263, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20768, 0, 3,
                                                                       19916, 5867, 19952, 303,
                                                                       313, 6293, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20828, 0, 3,
                                                                       19952, 5885, 19988, 313,
                                                                       323, 6323, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20888, 0, 3,
                                                                       19988, 5903, 20024, 323,
                                                                       333, 6353, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 20948, 0, 3,
                                                                       20024, 5921, 20060, 333,
                                                                       343, 6383, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21008, 0, 3,
                                                                       20060, 5939, 20096, 343,
                                                                       353, 6413, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 21068, 0, 3,
                                                                       20096, 5957, 20132, 353,
                                                                       363, 6443, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21128, 0, 3,
                                                                       20168, 5993, 20228, 383,
                                                                       398, 6473, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21218, 0, 3,
                                                                       20228, 6023, 20288, 398,
                                                                       413, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21308, 0, 3,
                                                                       20288, 6053, 20348, 413,
                                                                       428, 6563, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21398, 0, 3,
                                                                       20348, 6083, 20408, 428,
                                                                       443, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21488, 0, 3,
                                                                       20408, 6113, 20468, 443,
                                                                       458, 6653, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21578, 0, 3,
                                                                       20468, 6143, 20528, 458,
                                                                       473, 6698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21668, 0, 3,
                                                                       20528, 6173, 20588, 473,
                                                                       488, 6743, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21758, 0, 3,
                                                                       20588, 6203, 20648, 488,
                                                                       503, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21848, 0, 3,
                                                                       20648, 6233, 20708, 503,
                                                                       518, 6833, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 21938, 0, 3,
                                                                       20708, 6263, 20768, 518,
                                                                       533, 6878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22028, 0, 3,
                                                                       20768, 6293, 20828, 533,
                                                                       548, 6923, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22118, 0, 3,
                                                                       20828, 6323, 20888, 548,
                                                                       563, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22208, 0, 3,
                                                                       20888, 6353, 20948, 563,
                                                                       578, 7013, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22298, 0, 3,
                                                                       20948, 6383, 21008, 578,
                                                                       593, 7058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 22388, 0, 3,
                                                                       21008, 6413, 21068, 593,
                                                                       608, 7103, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22478, 0, 3,
                                                                       21128, 6473, 21218, 638,
                                                                       659, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22604, 0, 3,
                                                                       21218, 6518, 21308, 659,
                                                                       680, 7211, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22730, 0, 3,
                                                                       21308, 6563, 21398, 680,
                                                                       701, 7274, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22856, 0, 3,
                                                                       21398, 6608, 21488, 701,
                                                                       722, 7337, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 22982, 0, 3,
                                                                       21488, 6653, 21578, 722,
                                                                       743, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23108, 0, 3,
                                                                       21578, 6698, 21668, 743,
                                                                       764, 7463, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23234, 0, 3,
                                                                       21668, 6743, 21758, 764,
                                                                       785, 7526, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23360, 0, 3,
                                                                       21758, 6788, 21848, 785,
                                                                       806, 7589, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23486, 0, 3,
                                                                       21848, 6833, 21938, 806,
                                                                       827, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23612, 0, 3,
                                                                       21938, 6878, 22028, 827,
                                                                       848, 7715, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23738, 0, 3,
                                                                       22028, 6923, 22118, 848,
                                                                       869, 7778, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23864, 0, 3,
                                                                       22118, 6968, 22208, 869,
                                                                       890, 7841, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 23990, 0, 3,
                                                                       22208, 7013, 22298, 890,
                                                                       911, 7904, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 24116, 0, 3,
                                                                       22298, 7058, 22388, 911,
                                                                       932, 7967, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24242, 0, 3,
                                                                       22478, 7148, 22604, 974,
                                                                       1002, 8030, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24410, 0, 3,
                                                                       22604, 7211, 22730, 1002,
                                                                       1030, 8114, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24578, 0, 3,
                                                                       22730, 7274, 22856, 1030,
                                                                       1058, 8198, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24746, 0, 3,
                                                                       22856, 7337, 22982, 1058,
                                                                       1086, 8282, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 24914, 0, 3,
                                                                       22982, 7400, 23108, 1086,
                                                                       1114, 8366, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25082, 0, 3,
                                                                       23108, 7463, 23234, 1114,
                                                                       1142, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25250, 0, 3,
                                                                       23234, 7526, 23360, 1142,
                                                                       1170, 8534, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25418, 0, 3,
                                                                       23360, 7589, 23486, 1170,
                                                                       1198, 8618, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25586, 0, 3,
                                                                       23486, 7652, 23612, 1198,
                                                                       1226, 8702, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25754, 0, 3,
                                                                       23612, 7715, 23738, 1226,
                                                                       1254, 8786, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 25922, 0, 3,
                                                                       23738, 7778, 23864, 1254,
                                                                       1282, 8870, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26090, 0, 3,
                                                                       23864, 7841, 23990, 1282,
                                                                       1310, 8954, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 26258, 0, 3,
                                                                       23990, 7904, 24116, 1310,
                                                                       1338, 9038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26426, 0, 3,
                                                                       24242, 8030, 24410, 1394,
                                                                       1430, 9122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26642, 0, 3,
                                                                       24410, 8114, 24578, 1430,
                                                                       1466, 9230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 26858, 0, 3,
                                                                       24578, 8198, 24746, 1466,
                                                                       1502, 9338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27074, 0, 3,
                                                                       24746, 8282, 24914, 1502,
                                                                       1538, 9446, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27290, 0, 3,
                                                                       24914, 8366, 25082, 1538,
                                                                       1574, 9554, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27506, 0, 3,
                                                                       25082, 8450, 25250, 1574,
                                                                       1610, 9662, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27722, 0, 3,
                                                                       25250, 8534, 25418, 1610,
                                                                       1646, 9770, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 27938, 0, 3,
                                                                       25418, 8618, 25586, 1646,
                                                                       1682, 9878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28154, 0, 3,
                                                                       25586, 8702, 25754, 1682,
                                                                       1718, 9986, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28370, 0, 3,
                                                                       25754, 8786, 25922, 1718,
                                                                       1754, 10094, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28586, 0, 3,
                                                                       25922, 8870, 26090, 1754,
                                                                       1790, 10202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 28802, 0, 3,
                                                                       26090, 8954, 26258, 1790,
                                                                       1826, 10310, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29018, 0, 3,
                                                                       26426, 9122, 26642, 1898,
                                                                       1943, 10418, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29288, 0, 3,
                                                                       26642, 9230, 26858, 1943,
                                                                       1988, 10553, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29558, 0, 3,
                                                                       26858, 9338, 27074, 1988,
                                                                       2033, 10688, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 29828, 0, 3,
                                                                       27074, 9446, 27290, 2033,
                                                                       2078, 10823, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30098, 0, 3,
                                                                       27290, 9554, 27506, 2078,
                                                                       2123, 10958, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30368, 0, 3,
                                                                       27506, 9662, 27722, 2123,
                                                                       2168, 11093, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30638, 0, 3,
                                                                       27722, 9770, 27938, 2168,
                                                                       2213, 11228, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 30908, 0, 3,
                                                                       27938, 9878, 28154, 2213,
                                                                       2258, 11363, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31178, 0, 3,
                                                                       28154, 9986, 28370, 2258,
                                                                       2303, 11498, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31448, 0, 3,
                                                                       28370, 10094, 28586, 2303,
                                                                       2348, 11633, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 31718, 0, 3,
                                                                       28586, 10202, 28802, 2348,
                                                                       2393, 11768, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 31988, 0, 3,
                                                                       29018, 10418, 29288, 2483,
                                                                       2538, 11903, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32318, 0, 3,
                                                                       29288, 10553, 29558, 2538,
                                                                       2593, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32648, 0, 3,
                                                                       29558, 10688, 29828, 2593,
                                                                       2648, 12233, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 32978, 0, 3,
                                                                       29828, 10823, 30098, 2648,
                                                                       2703, 12398, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33308, 0, 3,
                                                                       30098, 10958, 30368, 2703,
                                                                       2758, 12563, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33638, 0, 3,
                                                                       30368, 11093, 30638, 2758,
                                                                       2813, 12728, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 33968, 0, 3,
                                                                       30638, 11228, 30908, 2813,
                                                                       2868, 12893, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34298, 0, 3,
                                                                       30908, 11363, 31178, 2868,
                                                                       2923, 13058, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34628, 0, 3,
                                                                       31178, 11498, 31448, 2923,
                                                                       2978, 13223, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 34958, 0, 3,
                                                                       31448, 11633, 31718, 2978,
                                                                       3033, 13388, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 35288, 0, 3,
                                                                       31988, 11903, 32318, 3143,
                                                                       3209, 13553, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 35684, 0, 3,
                                                                       32318, 12068, 32648, 3209,
                                                                       3275, 13751, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 36080, 0, 3,
                                                                       32648, 12233, 32978, 3275,
                                                                       3341, 13949, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 36476, 0, 3,
                                                                       32978, 12398, 33308, 3341,
                                                                       3407, 14147, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 36872, 0, 3,
                                                                       33308, 12563, 33638, 3407,
                                                                       3473, 14345, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 37268, 0, 3,
                                                                       33638, 12728, 33968, 3473,
                                                                       3539, 14543, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 37664, 0, 3,
                                                                       33968, 12893, 34298, 3539,
                                                                       3605, 14741, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38060, 0, 3,
                                                                       34298, 13058, 34628, 3605,
                                                                       3671, 14939, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 38456, 0, 3,
                                                                       34628, 13223, 34958, 3671,
                                                                       3737, 15137, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 38852, 0, 3,
                                                                       35288, 13553, 35684, 3869,
                                                                       3947, 15335, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 39320, 0, 3,
                                                                       35684, 13751, 36080, 3947,
                                                                       4025, 15569, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 39788, 0, 3,
                                                                       36080, 13949, 36476, 4025,
                                                                       4103, 15803, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 40256, 0, 3,
                                                                       36476, 14147, 36872, 4103,
                                                                       4181, 16037, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 40724, 0, 3,
                                                                       36872, 14345, 37268, 4181,
                                                                       4259, 16271, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 41192, 0, 3,
                                                                       37268, 14543, 37664, 4259,
                                                                       4337, 16505, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 41660, 0, 3,
                                                                       37664, 14741, 38060, 4337,
                                                                       4415, 16739, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 42128, 0, 3,
                                                                       38060, 14939, 38456, 4415,
                                                                       4493, 16973, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 42596, 0, 3,
                                                                       38852, 15335, 39320, 4649,
                                                                       4740, 17207, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 43142, 0, 3,
                                                                       39320, 15569, 39788, 4740,
                                                                       4831, 17480, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 43688, 0, 3,
                                                                       39788, 15803, 40256, 4831,
                                                                       4922, 17753, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 44234, 0, 3,
                                                                       40256, 16037, 40724, 4922,
                                                                       5013, 18026, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 44780, 0, 3,
                                                                       40724, 16271, 41192, 5013,
                                                                       5104, 18299, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 45326, 0, 3,
                                                                       41192, 16505, 41660, 5104,
                                                                       5195, 18572, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 45872, 0, 3,
                                                                       41660, 16739, 42128, 5195,
                                                                       5286, 18845, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46418, 3, 5468,
                                                                       5471, 19130, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46428, 3, 5471,
                                                                       5474, 19136, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46438, 3, 5474,
                                                                       5477, 19142, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46448, 3, 5477,
                                                                       5480, 19148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46458, 3, 5480,
                                                                       5483, 19154, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46468, 3, 5483,
                                                                       5486, 19160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46478, 3, 5486,
                                                                       5489, 19166, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46488, 3, 5489,
                                                                       5492, 19172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46498, 3, 5492,
                                                                       5495, 19178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46508, 3, 5495,
                                                                       5498, 19184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46518, 3, 5498,
                                                                       5501, 19190, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46528, 3, 5501,
                                                                       5504, 19196, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46538, 3, 5504,
                                                                       5507, 19202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46548, 3, 5507,
                                                                       5510, 19208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46558, 3, 5510,
                                                                       5513, 19214, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46568, 3, 5513,
                                                                       5516, 19220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 46578, 3, 5516,
                                                                       5519, 19226, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46588, 0, 3,
                                                                       46418, 19130, 46428, 5525,
                                                                       5534, 19268, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46618, 0, 3,
                                                                       46428, 19136, 46438, 5534,
                                                                       5543, 19286, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46648, 0, 3,
                                                                       46438, 19142, 46448, 5543,
                                                                       5552, 19304, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46678, 0, 3,
                                                                       46448, 19148, 46458, 5552,
                                                                       5561, 19322, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46708, 0, 3,
                                                                       46458, 19154, 46468, 5561,
                                                                       5570, 19340, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46738, 0, 3,
                                                                       46468, 19160, 46478, 5570,
                                                                       5579, 19358, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46768, 0, 3,
                                                                       46478, 19166, 46488, 5579,
                                                                       5588, 19376, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46798, 0, 3,
                                                                       46488, 19172, 46498, 5588,
                                                                       5597, 19394, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46828, 0, 3,
                                                                       46498, 19178, 46508, 5597,
                                                                       5606, 19412, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46858, 0, 3,
                                                                       46508, 19184, 46518, 5606,
                                                                       5615, 19430, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46888, 0, 3,
                                                                       46518, 19190, 46528, 5615,
                                                                       5624, 19448, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46918, 0, 3,
                                                                       46528, 19196, 46538, 5624,
                                                                       5633, 19466, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46948, 0, 3,
                                                                       46538, 19202, 46548, 5633,
                                                                       5642, 19484, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 46978, 0, 3,
                                                                       46548, 19208, 46558, 5642,
                                                                       5651, 19502, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47008, 0, 3,
                                                                       46558, 19214, 46568, 5651,
                                                                       5660, 19520, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 47038, 0, 3,
                                                                       46568, 19220, 46578, 5660,
                                                                       5669, 19538, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47068, 0, 3,
                                                                       46588, 19268, 46618, 5687,
                                                                       5705, 19628, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47128, 0, 3,
                                                                       46618, 19286, 46648, 5705,
                                                                       5723, 19664, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47188, 0, 3,
                                                                       46648, 19304, 46678, 5723,
                                                                       5741, 19700, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47248, 0, 3,
                                                                       46678, 19322, 46708, 5741,
                                                                       5759, 19736, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47308, 0, 3,
                                                                       46708, 19340, 46738, 5759,
                                                                       5777, 19772, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47368, 0, 3,
                                                                       46738, 19358, 46768, 5777,
                                                                       5795, 19808, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47428, 0, 3,
                                                                       46768, 19376, 46798, 5795,
                                                                       5813, 19844, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47488, 0, 3,
                                                                       46798, 19394, 46828, 5813,
                                                                       5831, 19880, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47548, 0, 3,
                                                                       46828, 19412, 46858, 5831,
                                                                       5849, 19916, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47608, 0, 3,
                                                                       46858, 19430, 46888, 5849,
                                                                       5867, 19952, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47668, 0, 3,
                                                                       46888, 19448, 46918, 5867,
                                                                       5885, 19988, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47728, 0, 3,
                                                                       46918, 19466, 46948, 5885,
                                                                       5903, 20024, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47788, 0, 3,
                                                                       46948, 19484, 46978, 5903,
                                                                       5921, 20060, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47848, 0, 3,
                                                                       46978, 19502, 47008, 5921,
                                                                       5939, 20096, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 47908, 0, 3,
                                                                       47008, 19520, 47038, 5939,
                                                                       5957, 20132, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 47968, 0, 3,
                                                                       47068, 19628, 47128, 5993,
                                                                       6023, 20288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48068, 0, 3,
                                                                       47128, 19664, 47188, 6023,
                                                                       6053, 20348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48168, 0, 3,
                                                                       47188, 19700, 47248, 6053,
                                                                       6083, 20408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48268, 0, 3,
                                                                       47248, 19736, 47308, 6083,
                                                                       6113, 20468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48368, 0, 3,
                                                                       47308, 19772, 47368, 6113,
                                                                       6143, 20528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48468, 0, 3,
                                                                       47368, 19808, 47428, 6143,
                                                                       6173, 20588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48568, 0, 3,
                                                                       47428, 19844, 47488, 6173,
                                                                       6203, 20648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48668, 0, 3,
                                                                       47488, 19880, 47548, 6203,
                                                                       6233, 20708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48768, 0, 3,
                                                                       47548, 19916, 47608, 6233,
                                                                       6263, 20768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48868, 0, 3,
                                                                       47608, 19952, 47668, 6263,
                                                                       6293, 20828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 48968, 0, 3,
                                                                       47668, 19988, 47728, 6293,
                                                                       6323, 20888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 49068, 0, 3,
                                                                       47728, 20024, 47788, 6323,
                                                                       6353, 20948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 49168, 0, 3,
                                                                       47788, 20060, 47848, 6353,
                                                                       6383, 21008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 49268, 0, 3,
                                                                       47848, 20096, 47908, 6383,
                                                                       6413, 21068, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49368, 0, 3,
                                                                       47968, 20288, 48068, 6473,
                                                                       6518, 21308, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49518, 0, 3,
                                                                       48068, 20348, 48168, 6518,
                                                                       6563, 21398, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49668, 0, 3,
                                                                       48168, 20408, 48268, 6563,
                                                                       6608, 21488, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49818, 0, 3,
                                                                       48268, 20468, 48368, 6608,
                                                                       6653, 21578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 49968, 0, 3,
                                                                       48368, 20528, 48468, 6653,
                                                                       6698, 21668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50118, 0, 3,
                                                                       48468, 20588, 48568, 6698,
                                                                       6743, 21758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50268, 0, 3,
                                                                       48568, 20648, 48668, 6743,
                                                                       6788, 21848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50418, 0, 3,
                                                                       48668, 20708, 48768, 6788,
                                                                       6833, 21938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50568, 0, 3,
                                                                       48768, 20768, 48868, 6833,
                                                                       6878, 22028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50718, 0, 3,
                                                                       48868, 20828, 48968, 6878,
                                                                       6923, 22118, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 50868, 0, 3,
                                                                       48968, 20888, 49068, 6923,
                                                                       6968, 22208, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 51018, 0, 3,
                                                                       49068, 20948, 49168, 6968,
                                                                       7013, 22298, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 51168, 0, 3,
                                                                       49168, 21008, 49268, 7013,
                                                                       7058, 22388, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51318, 0, 3,
                                                                       49368, 21308, 49518, 7148,
                                                                       7211, 22730, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51528, 0, 3,
                                                                       49518, 21398, 49668, 7211,
                                                                       7274, 22856, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51738, 0, 3,
                                                                       49668, 21488, 49818, 7274,
                                                                       7337, 22982, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 51948, 0, 3,
                                                                       49818, 21578, 49968, 7337,
                                                                       7400, 23108, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52158, 0, 3,
                                                                       49968, 21668, 50118, 7400,
                                                                       7463, 23234, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52368, 0, 3,
                                                                       50118, 21758, 50268, 7463,
                                                                       7526, 23360, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52578, 0, 3,
                                                                       50268, 21848, 50418, 7526,
                                                                       7589, 23486, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52788, 0, 3,
                                                                       50418, 21938, 50568, 7589,
                                                                       7652, 23612, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 52998, 0, 3,
                                                                       50568, 22028, 50718, 7652,
                                                                       7715, 23738, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53208, 0, 3,
                                                                       50718, 22118, 50868, 7715,
                                                                       7778, 23864, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53418, 0, 3,
                                                                       50868, 22208, 51018, 7778,
                                                                       7841, 23990, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 53628, 0, 3,
                                                                       51018, 22298, 51168, 7841,
                                                                       7904, 24116, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 53838, 0, 3,
                                                                       51318, 22730, 51528, 8030,
                                                                       8114, 24578, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54118, 0, 3,
                                                                       51528, 22856, 51738, 8114,
                                                                       8198, 24746, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54398, 0, 3,
                                                                       51738, 22982, 51948, 8198,
                                                                       8282, 24914, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54678, 0, 3,
                                                                       51948, 23108, 52158, 8282,
                                                                       8366, 25082, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 54958, 0, 3,
                                                                       52158, 23234, 52368, 8366,
                                                                       8450, 25250, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55238, 0, 3,
                                                                       52368, 23360, 52578, 8450,
                                                                       8534, 25418, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55518, 0, 3,
                                                                       52578, 23486, 52788, 8534,
                                                                       8618, 25586, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 55798, 0, 3,
                                                                       52788, 23612, 52998, 8618,
                                                                       8702, 25754, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56078, 0, 3,
                                                                       52998, 23738, 53208, 8702,
                                                                       8786, 25922, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56358, 0, 3,
                                                                       53208, 23864, 53418, 8786,
                                                                       8870, 26090, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 56638, 0, 3,
                                                                       53418, 23990, 53628, 8870,
                                                                       8954, 26258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 56918, 0, 3,
                                                                       53838, 24578, 54118, 9122,
                                                                       9230, 26858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57278, 0, 3,
                                                                       54118, 24746, 54398, 9230,
                                                                       9338, 27074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57638, 0, 3,
                                                                       54398, 24914, 54678, 9338,
                                                                       9446, 27290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 57998, 0, 3,
                                                                       54678, 25082, 54958, 9446,
                                                                       9554, 27506, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58358, 0, 3,
                                                                       54958, 25250, 55238, 9554,
                                                                       9662, 27722, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 58718, 0, 3,
                                                                       55238, 25418, 55518, 9662,
                                                                       9770, 27938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59078, 0, 3,
                                                                       55518, 25586, 55798, 9770,
                                                                       9878, 28154, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59438, 0, 3,
                                                                       55798, 25754, 56078, 9878,
                                                                       9986, 28370, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 59798, 0, 3,
                                                                       56078, 25922, 56358, 9986,
                                                                       10094, 28586, ncols,
                                                                       gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 60158, 0, 3,
                                                                       56358, 26090, 56638,
                                                                       10094, 10202, 28802,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60518, 0, 3,
                                                                       56918, 26858, 57278,
                                                                       10418, 10553, 29558,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 60968, 0, 3,
                                                                       57278, 27074, 57638,
                                                                       10553, 10688, 29828,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 61418, 0, 3,
                                                                       57638, 27290, 57998,
                                                                       10688, 10823, 30098,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 61868, 0, 3,
                                                                       57998, 27506, 58358,
                                                                       10823, 10958, 30368,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62318, 0, 3,
                                                                       58358, 27722, 58718,
                                                                       10958, 11093, 30638,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 62768, 0, 3,
                                                                       58718, 27938, 59078,
                                                                       11093, 11228, 30908,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63218, 0, 3,
                                                                       59078, 28154, 59438,
                                                                       11228, 11363, 31178,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 63668, 0, 3,
                                                                       59438, 28370, 59798,
                                                                       11363, 11498, 31448,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 64118, 0, 3,
                                                                       59798, 28586, 60158,
                                                                       11498, 11633, 31718,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 64568, 0, 3,
                                                                       60518, 29558, 60968,
                                                                       11903, 12068, 32648,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 65118, 0, 3,
                                                                       60968, 29828, 61418,
                                                                       12068, 12233, 32978,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 65668, 0, 3,
                                                                       61418, 30098, 61868,
                                                                       12233, 12398, 33308,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 66218, 0, 3,
                                                                       61868, 30368, 62318,
                                                                       12398, 12563, 33638,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 66768, 0, 3,
                                                                       62318, 30638, 62768,
                                                                       12563, 12728, 33968,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67318, 0, 3,
                                                                       62768, 30908, 63218,
                                                                       12728, 12893, 34298,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 67868, 0, 3,
                                                                       63218, 31178, 63668,
                                                                       12893, 13058, 34628,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 68418, 0, 3,
                                                                       63668, 31448, 64118,
                                                                       13058, 13223, 34958,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 68968, 0, 3,
                                                                       64568, 32648, 65118,
                                                                       13553, 13751, 36080,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 69628, 0, 3,
                                                                       65118, 32978, 65668,
                                                                       13751, 13949, 36476,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 70288, 0, 3,
                                                                       65668, 33308, 66218,
                                                                       13949, 14147, 36872,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 70948, 0, 3,
                                                                       66218, 33638, 66768,
                                                                       14147, 14345, 37268,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 71608, 0, 3,
                                                                       66768, 33968, 67318,
                                                                       14345, 14543, 37664,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 72268, 0, 3,
                                                                       67318, 34298, 67868,
                                                                       14543, 14741, 38060,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 72928, 0, 3,
                                                                       67868, 34628, 68418,
                                                                       14741, 14939, 38456,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 73588, 0, 3,
                                                                       68968, 36080, 69628,
                                                                       15335, 15569, 39788,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 74368, 0, 3,
                                                                       69628, 36476, 70288,
                                                                       15569, 15803, 40256,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 75148, 0, 3,
                                                                       70288, 36872, 70948,
                                                                       15803, 16037, 40724,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 75928, 0, 3,
                                                                       70948, 37268, 71608,
                                                                       16037, 16271, 41192,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 76708, 0, 3,
                                                                       71608, 37664, 72268,
                                                                       16271, 16505, 41660,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 77488, 0, 3,
                                                                       72268, 38060, 72928,
                                                                       16505, 16739, 42128,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 78268, 0, 3,
                                                                       73588, 39788, 74368,
                                                                       17207, 17480, 43688,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 79178, 0, 3,
                                                                       74368, 40256, 75148,
                                                                       17480, 17753, 44234,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 80088, 0, 3,
                                                                       75148, 40724, 75928,
                                                                       17753, 18026, 44780,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 80998, 0, 3,
                                                                       75928, 41192, 76708,
                                                                       18026, 18299, 45326,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 81908, 0, 3,
                                                                       76708, 41660, 77488,
                                                                       18299, 18572, 45872,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82818, 3, 19118,
                                                                       19124, 46418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82833, 3, 19124,
                                                                       19130, 46428, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82848, 3, 19130,
                                                                       19136, 46438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82863, 3, 19136,
                                                                       19142, 46448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82878, 3, 19142,
                                                                       19148, 46458, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82893, 3, 19148,
                                                                       19154, 46468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82908, 3, 19154,
                                                                       19160, 46478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82923, 3, 19160,
                                                                       19166, 46488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82938, 3, 19166,
                                                                       19172, 46498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82953, 3, 19172,
                                                                       19178, 46508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82968, 3, 19178,
                                                                       19184, 46518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82983, 3, 19184,
                                                                       19190, 46528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 82998, 3, 19190,
                                                                       19196, 46538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 83013, 3, 19196,
                                                                       19202, 46548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 83028, 3, 19202,
                                                                       19208, 46558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 83043, 3, 19208,
                                                                       19214, 46568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 83058, 3, 19214,
                                                                       19220, 46578, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83073, 0, 3,
                                                                       82818, 46418, 82833,
                                                                       19232, 19250, 46588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83118, 0, 3,
                                                                       82833, 46428, 82848,
                                                                       19250, 19268, 46618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83163, 0, 3,
                                                                       82848, 46438, 82863,
                                                                       19268, 19286, 46648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83208, 0, 3,
                                                                       82863, 46448, 82878,
                                                                       19286, 19304, 46678,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83253, 0, 3,
                                                                       82878, 46458, 82893,
                                                                       19304, 19322, 46708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83298, 0, 3,
                                                                       82893, 46468, 82908,
                                                                       19322, 19340, 46738,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83343, 0, 3,
                                                                       82908, 46478, 82923,
                                                                       19340, 19358, 46768,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83388, 0, 3,
                                                                       82923, 46488, 82938,
                                                                       19358, 19376, 46798,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83433, 0, 3,
                                                                       82938, 46498, 82953,
                                                                       19376, 19394, 46828,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83478, 0, 3,
                                                                       82953, 46508, 82968,
                                                                       19394, 19412, 46858,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83523, 0, 3,
                                                                       82968, 46518, 82983,
                                                                       19412, 19430, 46888,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83568, 0, 3,
                                                                       82983, 46528, 82998,
                                                                       19430, 19448, 46918,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83613, 0, 3,
                                                                       82998, 46538, 83013,
                                                                       19448, 19466, 46948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83658, 0, 3,
                                                                       83013, 46548, 83028,
                                                                       19466, 19484, 46978,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83703, 0, 3,
                                                                       83028, 46558, 83043,
                                                                       19484, 19502, 47008,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 83748, 0, 3,
                                                                       83043, 46568, 83058,
                                                                       19502, 19520, 47038,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 83793, 0, 3,
                                                                       83073, 46588, 83118,
                                                                       19556, 19592, 47068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 83883, 0, 3,
                                                                       83118, 46618, 83163,
                                                                       19592, 19628, 47128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 83973, 0, 3,
                                                                       83163, 46648, 83208,
                                                                       19628, 19664, 47188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84063, 0, 3,
                                                                       83208, 46678, 83253,
                                                                       19664, 19700, 47248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84153, 0, 3,
                                                                       83253, 46708, 83298,
                                                                       19700, 19736, 47308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84243, 0, 3,
                                                                       83298, 46738, 83343,
                                                                       19736, 19772, 47368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84333, 0, 3,
                                                                       83343, 46768, 83388,
                                                                       19772, 19808, 47428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84423, 0, 3,
                                                                       83388, 46798, 83433,
                                                                       19808, 19844, 47488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84513, 0, 3,
                                                                       83433, 46828, 83478,
                                                                       19844, 19880, 47548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84603, 0, 3,
                                                                       83478, 46858, 83523,
                                                                       19880, 19916, 47608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84693, 0, 3,
                                                                       83523, 46888, 83568,
                                                                       19916, 19952, 47668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84783, 0, 3,
                                                                       83568, 46918, 83613,
                                                                       19952, 19988, 47728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84873, 0, 3,
                                                                       83613, 46948, 83658,
                                                                       19988, 20024, 47788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 84963, 0, 3,
                                                                       83658, 46978, 83703,
                                                                       20024, 20060, 47848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 85053, 0, 3,
                                                                       83703, 47008, 83748,
                                                                       20060, 20096, 47908,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85143, 0, 3,
                                                                       83793, 47068, 83883,
                                                                       20168, 20228, 47968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85293, 0, 3,
                                                                       83883, 47128, 83973,
                                                                       20228, 20288, 48068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85443, 0, 3,
                                                                       83973, 47188, 84063,
                                                                       20288, 20348, 48168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85593, 0, 3,
                                                                       84063, 47248, 84153,
                                                                       20348, 20408, 48268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85743, 0, 3,
                                                                       84153, 47308, 84243,
                                                                       20408, 20468, 48368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 85893, 0, 3,
                                                                       84243, 47368, 84333,
                                                                       20468, 20528, 48468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86043, 0, 3,
                                                                       84333, 47428, 84423,
                                                                       20528, 20588, 48568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86193, 0, 3,
                                                                       84423, 47488, 84513,
                                                                       20588, 20648, 48668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86343, 0, 3,
                                                                       84513, 47548, 84603,
                                                                       20648, 20708, 48768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86493, 0, 3,
                                                                       84603, 47608, 84693,
                                                                       20708, 20768, 48868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86643, 0, 3,
                                                                       84693, 47668, 84783,
                                                                       20768, 20828, 48968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86793, 0, 3,
                                                                       84783, 47728, 84873,
                                                                       20828, 20888, 49068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 86943, 0, 3,
                                                                       84873, 47788, 84963,
                                                                       20888, 20948, 49168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 87093, 0, 3,
                                                                       84963, 47848, 85053,
                                                                       20948, 21008, 49268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 87243, 0, 3,
                                                                       85143, 47968, 85293,
                                                                       21128, 21218, 49368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 87468, 0, 3,
                                                                       85293, 48068, 85443,
                                                                       21218, 21308, 49518,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 87693, 0, 3,
                                                                       85443, 48168, 85593,
                                                                       21308, 21398, 49668,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 87918, 0, 3,
                                                                       85593, 48268, 85743,
                                                                       21398, 21488, 49818,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 88143, 0, 3,
                                                                       85743, 48368, 85893,
                                                                       21488, 21578, 49968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 88368, 0, 3,
                                                                       85893, 48468, 86043,
                                                                       21578, 21668, 50118,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 88593, 0, 3,
                                                                       86043, 48568, 86193,
                                                                       21668, 21758, 50268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 88818, 0, 3,
                                                                       86193, 48668, 86343,
                                                                       21758, 21848, 50418,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89043, 0, 3,
                                                                       86343, 48768, 86493,
                                                                       21848, 21938, 50568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89268, 0, 3,
                                                                       86493, 48868, 86643,
                                                                       21938, 22028, 50718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89493, 0, 3,
                                                                       86643, 48968, 86793,
                                                                       22028, 22118, 50868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89718, 0, 3,
                                                                       86793, 49068, 86943,
                                                                       22118, 22208, 51018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 89943, 0, 3,
                                                                       86943, 49168, 87093,
                                                                       22208, 22298, 51168,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90168, 0, 3,
                                                                       87243, 49368, 87468,
                                                                       22478, 22604, 51318,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90483, 0, 3,
                                                                       87468, 49518, 87693,
                                                                       22604, 22730, 51528,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 90798, 0, 3,
                                                                       87693, 49668, 87918,
                                                                       22730, 22856, 51738,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 91113, 0, 3,
                                                                       87918, 49818, 88143,
                                                                       22856, 22982, 51948,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 91428, 0, 3,
                                                                       88143, 49968, 88368,
                                                                       22982, 23108, 52158,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 91743, 0, 3,
                                                                       88368, 50118, 88593,
                                                                       23108, 23234, 52368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 92058, 0, 3,
                                                                       88593, 50268, 88818,
                                                                       23234, 23360, 52578,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 92373, 0, 3,
                                                                       88818, 50418, 89043,
                                                                       23360, 23486, 52788,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 92688, 0, 3,
                                                                       89043, 50568, 89268,
                                                                       23486, 23612, 52998,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 93003, 0, 3,
                                                                       89268, 50718, 89493,
                                                                       23612, 23738, 53208,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 93318, 0, 3,
                                                                       89493, 50868, 89718,
                                                                       23738, 23864, 53418,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 93633, 0, 3,
                                                                       89718, 51018, 89943,
                                                                       23864, 23990, 53628,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 93948, 0, 3,
                                                                       90168, 51318, 90483,
                                                                       24242, 24410, 53838,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94368, 0, 3,
                                                                       90483, 51528, 90798,
                                                                       24410, 24578, 54118,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 94788, 0, 3,
                                                                       90798, 51738, 91113,
                                                                       24578, 24746, 54398,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 95208, 0, 3,
                                                                       91113, 51948, 91428,
                                                                       24746, 24914, 54678,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 95628, 0, 3,
                                                                       91428, 52158, 91743,
                                                                       24914, 25082, 54958,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96048, 0, 3,
                                                                       91743, 52368, 92058,
                                                                       25082, 25250, 55238,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96468, 0, 3,
                                                                       92058, 52578, 92373,
                                                                       25250, 25418, 55518,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 96888, 0, 3,
                                                                       92373, 52788, 92688,
                                                                       25418, 25586, 55798,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 97308, 0, 3,
                                                                       92688, 52998, 93003,
                                                                       25586, 25754, 56078,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 97728, 0, 3,
                                                                       93003, 53208, 93318,
                                                                       25754, 25922, 56358,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 98148, 0, 3,
                                                                       93318, 53418, 93633,
                                                                       25922, 26090, 56638,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 98568, 0, 3,
                                                                       93948, 53838, 94368,
                                                                       26426, 26642, 56918,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 99108, 0, 3,
                                                                       94368, 54118, 94788,
                                                                       26642, 26858, 57278,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 99648, 0, 3,
                                                                       94788, 54398, 95208,
                                                                       26858, 27074, 57638,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 100188, 0, 3,
                                                                       95208, 54678, 95628,
                                                                       27074, 27290, 57998,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 100728, 0, 3,
                                                                       95628, 54958, 96048,
                                                                       27290, 27506, 58358,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 101268, 0, 3,
                                                                       96048, 55238, 96468,
                                                                       27506, 27722, 58718,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 101808, 0, 3,
                                                                       96468, 55518, 96888,
                                                                       27722, 27938, 59078,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 102348, 0, 3,
                                                                       96888, 55798, 97308,
                                                                       27938, 28154, 59438,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 102888, 0, 3,
                                                                       97308, 56078, 97728,
                                                                       28154, 28370, 59798,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 103428, 0, 3,
                                                                       97728, 56358, 98148,
                                                                       28370, 28586, 60158,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 103968, 0, 3,
                                                                       98568, 56918, 99108,
                                                                       29018, 29288, 60518,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 104643, 0, 3,
                                                                       99108, 57278, 99648,
                                                                       29288, 29558, 60968,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 105318, 0, 3,
                                                                       99648, 57638, 100188,
                                                                       29558, 29828, 61418,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 105993, 0, 3,
                                                                       100188, 57998, 100728,
                                                                       29828, 30098, 61868,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 106668, 0, 3,
                                                                       100728, 58358, 101268,
                                                                       30098, 30368, 62318,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 107343, 0, 3,
                                                                       101268, 58718, 101808,
                                                                       30368, 30638, 62768,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 108018, 0, 3,
                                                                       101808, 59078, 102348,
                                                                       30638, 30908, 63218,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 108693, 0, 3,
                                                                       102348, 59438, 102888,
                                                                       30908, 31178, 63668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 109368, 0, 3,
                                                                       102888, 59798, 103428,
                                                                       31178, 31448, 64118,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 110043, 0, 3,
                                                                       103968, 60518, 104643,
                                                                       31988, 32318, 64568,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 110868, 0, 3,
                                                                       104643, 60968, 105318,
                                                                       32318, 32648, 65118,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 111693, 0, 3,
                                                                       105318, 61418, 105993,
                                                                       32648, 32978, 65668,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 112518, 0, 3,
                                                                       105993, 61868, 106668,
                                                                       32978, 33308, 66218,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 113343, 0, 3,
                                                                       106668, 62318, 107343,
                                                                       33308, 33638, 66768,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 114168, 0, 3,
                                                                       107343, 62768, 108018,
                                                                       33638, 33968, 67318,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 114993, 0, 3,
                                                                       108018, 63218, 108693,
                                                                       33968, 34298, 67868,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 115818, 0, 3,
                                                                       108693, 63668, 109368,
                                                                       34298, 34628, 68418,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 116643, 0, 3,
                                                                       110043, 64568, 110868,
                                                                       35288, 35684, 68968,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 117633, 0, 3,
                                                                       110868, 65118, 111693,
                                                                       35684, 36080, 69628,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 118623, 0, 3,
                                                                       111693, 65668, 112518,
                                                                       36080, 36476, 70288,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 119613, 0, 3,
                                                                       112518, 66218, 113343,
                                                                       36476, 36872, 70948,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 120603, 0, 3,
                                                                       113343, 66768, 114168,
                                                                       36872, 37268, 71608,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 121593, 0, 3,
                                                                       114168, 67318, 114993,
                                                                       37268, 37664, 72268,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 122583, 0, 3,
                                                                       114993, 67868, 115818,
                                                                       37664, 38060, 72928,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 123573, 0, 3,
                                                                       116643, 68968, 117633,
                                                                       38852, 39320, 73588,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 124743, 0, 3,
                                                                       117633, 69628, 118623,
                                                                       39320, 39788, 74368,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 125913, 0, 3,
                                                                       118623, 70288, 119613,
                                                                       39788, 40256, 75148,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 127083, 0, 3,
                                                                       119613, 70948, 120603,
                                                                       40256, 40724, 75928,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 128253, 0, 3,
                                                                       120603, 71608, 121593,
                                                                       40724, 41192, 76708,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 129423, 0, 3,
                                                                       121593, 72268, 122583,
                                                                       41192, 41660, 77488,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 130593, 0, 3,
                                                                       123573, 73588, 124743,
                                                                       42596, 43142, 78268,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 131958, 0, 3,
                                                                       124743, 74368, 125913,
                                                                       43142, 43688, 79178,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 133323, 0, 3,
                                                                       125913, 75148, 127083,
                                                                       43688, 44234, 80088,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 134688, 0, 3,
                                                                       127083, 75928, 128253,
                                                                       44234, 44780, 80998,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 136053, 0, 3,
                                                                       128253, 76708, 129423,
                                                                       44780, 45326, 81908,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137418, 3, 46418,
                                                                       46428, 82848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137439, 3, 46428,
                                                                       46438, 82863, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137460, 3, 46438,
                                                                       46448, 82878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137481, 3, 46448,
                                                                       46458, 82893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137502, 3, 46458,
                                                                       46468, 82908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137523, 3, 46468,
                                                                       46478, 82923, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137544, 3, 46478,
                                                                       46488, 82938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137565, 3, 46488,
                                                                       46498, 82953, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137586, 3, 46498,
                                                                       46508, 82968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137607, 3, 46508,
                                                                       46518, 82983, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137628, 3, 46518,
                                                                       46528, 82998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137649, 3, 46528,
                                                                       46538, 83013, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137670, 3, 46538,
                                                                       46548, 83028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137691, 3, 46548,
                                                                       46558, 83043, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137712, 3, 46558,
                                                                       46568, 83058, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137733, 0, 3,
                                                                       137418, 82848, 137439,
                                                                       46588, 46618, 83163,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137796, 0, 3,
                                                                       137439, 82863, 137460,
                                                                       46618, 46648, 83208,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137859, 0, 3,
                                                                       137460, 82878, 137481,
                                                                       46648, 46678, 83253,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137922, 0, 3,
                                                                       137481, 82893, 137502,
                                                                       46678, 46708, 83298,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 137985, 0, 3,
                                                                       137502, 82908, 137523,
                                                                       46708, 46738, 83343,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138048, 0, 3,
                                                                       137523, 82923, 137544,
                                                                       46738, 46768, 83388,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138111, 0, 3,
                                                                       137544, 82938, 137565,
                                                                       46768, 46798, 83433,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138174, 0, 3,
                                                                       137565, 82953, 137586,
                                                                       46798, 46828, 83478,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138237, 0, 3,
                                                                       137586, 82968, 137607,
                                                                       46828, 46858, 83523,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138300, 0, 3,
                                                                       137607, 82983, 137628,
                                                                       46858, 46888, 83568,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138363, 0, 3,
                                                                       137628, 82998, 137649,
                                                                       46888, 46918, 83613,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138426, 0, 3,
                                                                       137649, 83013, 137670,
                                                                       46918, 46948, 83658,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138489, 0, 3,
                                                                       137670, 83028, 137691,
                                                                       46948, 46978, 83703,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 138552, 0, 3,
                                                                       137691, 83043, 137712,
                                                                       46978, 47008, 83748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 138615, 0, 3,
                                                                       137733, 83163, 137796,
                                                                       47068, 47128, 83973,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 138741, 0, 3,
                                                                       137796, 83208, 137859,
                                                                       47128, 47188, 84063,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 138867, 0, 3,
                                                                       137859, 83253, 137922,
                                                                       47188, 47248, 84153,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 138993, 0, 3,
                                                                       137922, 83298, 137985,
                                                                       47248, 47308, 84243,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139119, 0, 3,
                                                                       137985, 83343, 138048,
                                                                       47308, 47368, 84333,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139245, 0, 3,
                                                                       138048, 83388, 138111,
                                                                       47368, 47428, 84423,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139371, 0, 3,
                                                                       138111, 83433, 138174,
                                                                       47428, 47488, 84513,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139497, 0, 3,
                                                                       138174, 83478, 138237,
                                                                       47488, 47548, 84603,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139623, 0, 3,
                                                                       138237, 83523, 138300,
                                                                       47548, 47608, 84693,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139749, 0, 3,
                                                                       138300, 83568, 138363,
                                                                       47608, 47668, 84783,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 139875, 0, 3,
                                                                       138363, 83613, 138426,
                                                                       47668, 47728, 84873,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140001, 0, 3,
                                                                       138426, 83658, 138489,
                                                                       47728, 47788, 84963,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 140127, 0, 3,
                                                                       138489, 83703, 138552,
                                                                       47788, 47848, 85053,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 140253, 0, 3,
                                                                       138615, 83973, 138741,
                                                                       47968, 48068, 85443,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 140463, 0, 3,
                                                                       138741, 84063, 138867,
                                                                       48068, 48168, 85593,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 140673, 0, 3,
                                                                       138867, 84153, 138993,
                                                                       48168, 48268, 85743,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 140883, 0, 3,
                                                                       138993, 84243, 139119,
                                                                       48268, 48368, 85893,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141093, 0, 3,
                                                                       139119, 84333, 139245,
                                                                       48368, 48468, 86043,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141303, 0, 3,
                                                                       139245, 84423, 139371,
                                                                       48468, 48568, 86193,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141513, 0, 3,
                                                                       139371, 84513, 139497,
                                                                       48568, 48668, 86343,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141723, 0, 3,
                                                                       139497, 84603, 139623,
                                                                       48668, 48768, 86493,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 141933, 0, 3,
                                                                       139623, 84693, 139749,
                                                                       48768, 48868, 86643,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142143, 0, 3,
                                                                       139749, 84783, 139875,
                                                                       48868, 48968, 86793,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142353, 0, 3,
                                                                       139875, 84873, 140001,
                                                                       48968, 49068, 86943,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 142563, 0, 3,
                                                                       140001, 84963, 140127,
                                                                       49068, 49168, 87093,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 142773, 0, 3,
                                                                       140253, 85443, 140463,
                                                                       49368, 49518, 87693,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 143088, 0, 3,
                                                                       140463, 85593, 140673,
                                                                       49518, 49668, 87918,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 143403, 0, 3,
                                                                       140673, 85743, 140883,
                                                                       49668, 49818, 88143,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 143718, 0, 3,
                                                                       140883, 85893, 141093,
                                                                       49818, 49968, 88368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 144033, 0, 3,
                                                                       141093, 86043, 141303,
                                                                       49968, 50118, 88593,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 144348, 0, 3,
                                                                       141303, 86193, 141513,
                                                                       50118, 50268, 88818,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 144663, 0, 3,
                                                                       141513, 86343, 141723,
                                                                       50268, 50418, 89043,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 144978, 0, 3,
                                                                       141723, 86493, 141933,
                                                                       50418, 50568, 89268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 145293, 0, 3,
                                                                       141933, 86643, 142143,
                                                                       50568, 50718, 89493,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 145608, 0, 3,
                                                                       142143, 86793, 142353,
                                                                       50718, 50868, 89718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 145923, 0, 3,
                                                                       142353, 86943, 142563,
                                                                       50868, 51018, 89943,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 146238, 0, 3,
                                                                       142773, 87693, 143088,
                                                                       51318, 51528, 90798,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 146679, 0, 3,
                                                                       143088, 87918, 143403,
                                                                       51528, 51738, 91113,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 147120, 0, 3,
                                                                       143403, 88143, 143718,
                                                                       51738, 51948, 91428,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 147561, 0, 3,
                                                                       143718, 88368, 144033,
                                                                       51948, 52158, 91743,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 148002, 0, 3,
                                                                       144033, 88593, 144348,
                                                                       52158, 52368, 92058,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 148443, 0, 3,
                                                                       144348, 88818, 144663,
                                                                       52368, 52578, 92373,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 148884, 0, 3,
                                                                       144663, 89043, 144978,
                                                                       52578, 52788, 92688,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 149325, 0, 3,
                                                                       144978, 89268, 145293,
                                                                       52788, 52998, 93003,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 149766, 0, 3,
                                                                       145293, 89493, 145608,
                                                                       52998, 53208, 93318,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 150207, 0, 3,
                                                                       145608, 89718, 145923,
                                                                       53208, 53418, 93633,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 150648, 0, 3,
                                                                       146238, 90798, 146679,
                                                                       53838, 54118, 94788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 151236, 0, 3,
                                                                       146679, 91113, 147120,
                                                                       54118, 54398, 95208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 151824, 0, 3,
                                                                       147120, 91428, 147561,
                                                                       54398, 54678, 95628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 152412, 0, 3,
                                                                       147561, 91743, 148002,
                                                                       54678, 54958, 96048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 153000, 0, 3,
                                                                       148002, 92058, 148443,
                                                                       54958, 55238, 96468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 153588, 0, 3,
                                                                       148443, 92373, 148884,
                                                                       55238, 55518, 96888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 154176, 0, 3,
                                                                       148884, 92688, 149325,
                                                                       55518, 55798, 97308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 154764, 0, 3,
                                                                       149325, 93003, 149766,
                                                                       55798, 56078, 97728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 155352, 0, 3,
                                                                       149766, 93318, 150207,
                                                                       56078, 56358, 98148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 155940, 0, 3,
                                                                       150648, 94788, 151236,
                                                                       56918, 57278, 99648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 156696, 0, 3,
                                                                       151236, 95208, 151824,
                                                                       57278, 57638, 100188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 157452, 0, 3,
                                                                       151824, 95628, 152412,
                                                                       57638, 57998, 100728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 158208, 0, 3,
                                                                       152412, 96048, 153000,
                                                                       57998, 58358, 101268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 158964, 0, 3,
                                                                       153000, 96468, 153588,
                                                                       58358, 58718, 101808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 159720, 0, 3,
                                                                       153588, 96888, 154176,
                                                                       58718, 59078, 102348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 160476, 0, 3,
                                                                       154176, 97308, 154764,
                                                                       59078, 59438, 102888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 161232, 0, 3,
                                                                       154764, 97728, 155352,
                                                                       59438, 59798, 103428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 161988, 0, 3,
                                                                       155940, 99648, 156696,
                                                                       60518, 60968, 105318,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 162933, 0, 3,
                                                                       156696, 100188, 157452,
                                                                       60968, 61418, 105993,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 163878, 0, 3,
                                                                       157452, 100728, 158208,
                                                                       61418, 61868, 106668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 164823, 0, 3,
                                                                       158208, 101268, 158964,
                                                                       61868, 62318, 107343,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 165768, 0, 3,
                                                                       158964, 101808, 159720,
                                                                       62318, 62768, 108018,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 166713, 0, 3,
                                                                       159720, 102348, 160476,
                                                                       62768, 63218, 108693,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 167658, 0, 3,
                                                                       160476, 102888, 161232,
                                                                       63218, 63668, 109368,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 168603, 0, 3,
                                                                       161988, 105318, 162933,
                                                                       64568, 65118, 111693,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 169758, 0, 3,
                                                                       162933, 105993, 163878,
                                                                       65118, 65668, 112518,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 170913, 0, 3,
                                                                       163878, 106668, 164823,
                                                                       65668, 66218, 113343,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 172068, 0, 3,
                                                                       164823, 107343, 165768,
                                                                       66218, 66768, 114168,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 173223, 0, 3,
                                                                       165768, 108018, 166713,
                                                                       66768, 67318, 114993,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 174378, 0, 3,
                                                                       166713, 108693, 167658,
                                                                       67318, 67868, 115818,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 175533, 0, 3,
                                                                       168603, 111693, 169758,
                                                                       68968, 69628, 118623,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 176919, 0, 3,
                                                                       169758, 112518, 170913,
                                                                       69628, 70288, 119613,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 178305, 0, 3,
                                                                       170913, 113343, 172068,
                                                                       70288, 70948, 120603,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 179691, 0, 3,
                                                                       172068, 114168, 173223,
                                                                       70948, 71608, 121593,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 181077, 0, 3,
                                                                       173223, 114993, 174378,
                                                                       71608, 72268, 122583,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 182463, 0, 3,
                                                                       175533, 118623, 176919,
                                                                       73588, 74368, 125913,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 184101, 0, 3,
                                                                       176919, 119613, 178305,
                                                                       74368, 75148, 127083,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 185739, 0, 3,
                                                                       178305, 120603, 179691,
                                                                       75148, 75928, 128253,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 187377, 0, 3,
                                                                       179691, 121593, 181077,
                                                                       75928, 76708, 129423,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 189015, 0, 3,
                                                                       182463, 125913, 184101,
                                                                       78268, 79178, 133323,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 190926, 0, 3,
                                                                       184101, 127083, 185739,
                                                                       79178, 80088, 134688,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 192837, 0, 3,
                                                                       185739, 128253, 187377,
                                                                       80088, 80998, 136053,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194748, 3, 82818,
                                                                       82833, 137418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194776, 3, 82833,
                                                                       82848, 137439, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194804, 3, 82848,
                                                                       82863, 137460, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194832, 3, 82863,
                                                                       82878, 137481, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194860, 3, 82878,
                                                                       82893, 137502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194888, 3, 82893,
                                                                       82908, 137523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194916, 3, 82908,
                                                                       82923, 137544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194944, 3, 82923,
                                                                       82938, 137565, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 194972, 3, 82938,
                                                                       82953, 137586, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195000, 3, 82953,
                                                                       82968, 137607, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195028, 3, 82968,
                                                                       82983, 137628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195056, 3, 82983,
                                                                       82998, 137649, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195084, 3, 82998,
                                                                       83013, 137670, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195112, 3, 83013,
                                                                       83028, 137691, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 195140, 3, 83028,
                                                                       83043, 137712, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195168, 0, 3,
                                                                       194748, 137418, 194776,
                                                                       83073, 83118, 137733,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195252, 0, 3,
                                                                       194776, 137439, 194804,
                                                                       83118, 83163, 137796,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195336, 0, 3,
                                                                       194804, 137460, 194832,
                                                                       83163, 83208, 137859,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195420, 0, 3,
                                                                       194832, 137481, 194860,
                                                                       83208, 83253, 137922,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195504, 0, 3,
                                                                       194860, 137502, 194888,
                                                                       83253, 83298, 137985,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195588, 0, 3,
                                                                       194888, 137523, 194916,
                                                                       83298, 83343, 138048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195672, 0, 3,
                                                                       194916, 137544, 194944,
                                                                       83343, 83388, 138111,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195756, 0, 3,
                                                                       194944, 137565, 194972,
                                                                       83388, 83433, 138174,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195840, 0, 3,
                                                                       194972, 137586, 195000,
                                                                       83433, 83478, 138237,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 195924, 0, 3,
                                                                       195000, 137607, 195028,
                                                                       83478, 83523, 138300,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 196008, 0, 3,
                                                                       195028, 137628, 195056,
                                                                       83523, 83568, 138363,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 196092, 0, 3,
                                                                       195056, 137649, 195084,
                                                                       83568, 83613, 138426,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 196176, 0, 3,
                                                                       195084, 137670, 195112,
                                                                       83613, 83658, 138489,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 196260, 0, 3,
                                                                       195112, 137691, 195140,
                                                                       83658, 83703, 138552,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 196344, 0, 3,
                                                                       195168, 137733, 195252,
                                                                       83793, 83883, 138615,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 196512, 0, 3,
                                                                       195252, 137796, 195336,
                                                                       83883, 83973, 138741,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 196680, 0, 3,
                                                                       195336, 137859, 195420,
                                                                       83973, 84063, 138867,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 196848, 0, 3,
                                                                       195420, 137922, 195504,
                                                                       84063, 84153, 138993,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197016, 0, 3,
                                                                       195504, 137985, 195588,
                                                                       84153, 84243, 139119,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197184, 0, 3,
                                                                       195588, 138048, 195672,
                                                                       84243, 84333, 139245,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197352, 0, 3,
                                                                       195672, 138111, 195756,
                                                                       84333, 84423, 139371,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197520, 0, 3,
                                                                       195756, 138174, 195840,
                                                                       84423, 84513, 139497,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197688, 0, 3,
                                                                       195840, 138237, 195924,
                                                                       84513, 84603, 139623,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 197856, 0, 3,
                                                                       195924, 138300, 196008,
                                                                       84603, 84693, 139749,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 198024, 0, 3,
                                                                       196008, 138363, 196092,
                                                                       84693, 84783, 139875,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 198192, 0, 3,
                                                                       196092, 138426, 196176,
                                                                       84783, 84873, 140001,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 198360, 0, 3,
                                                                       196176, 138489, 196260,
                                                                       84873, 84963, 140127,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 198528, 0, 3,
                                                                       196344, 138615, 196512,
                                                                       85143, 85293, 140253,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 198808, 0, 3,
                                                                       196512, 138741, 196680,
                                                                       85293, 85443, 140463,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 199088, 0, 3,
                                                                       196680, 138867, 196848,
                                                                       85443, 85593, 140673,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 199368, 0, 3,
                                                                       196848, 138993, 197016,
                                                                       85593, 85743, 140883,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 199648, 0, 3,
                                                                       197016, 139119, 197184,
                                                                       85743, 85893, 141093,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 199928, 0, 3,
                                                                       197184, 139245, 197352,
                                                                       85893, 86043, 141303,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 200208, 0, 3,
                                                                       197352, 139371, 197520,
                                                                       86043, 86193, 141513,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 200488, 0, 3,
                                                                       197520, 139497, 197688,
                                                                       86193, 86343, 141723,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 200768, 0, 3,
                                                                       197688, 139623, 197856,
                                                                       86343, 86493, 141933,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 201048, 0, 3,
                                                                       197856, 139749, 198024,
                                                                       86493, 86643, 142143,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 201328, 0, 3,
                                                                       198024, 139875, 198192,
                                                                       86643, 86793, 142353,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 201608, 0, 3,
                                                                       198192, 140001, 198360,
                                                                       86793, 86943, 142563,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 201888, 0, 3,
                                                                       198528, 140253, 198808,
                                                                       87243, 87468, 142773,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 202308, 0, 3,
                                                                       198808, 140463, 199088,
                                                                       87468, 87693, 143088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 202728, 0, 3,
                                                                       199088, 140673, 199368,
                                                                       87693, 87918, 143403,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 203148, 0, 3,
                                                                       199368, 140883, 199648,
                                                                       87918, 88143, 143718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 203568, 0, 3,
                                                                       199648, 141093, 199928,
                                                                       88143, 88368, 144033,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 203988, 0, 3,
                                                                       199928, 141303, 200208,
                                                                       88368, 88593, 144348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 204408, 0, 3,
                                                                       200208, 141513, 200488,
                                                                       88593, 88818, 144663,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 204828, 0, 3,
                                                                       200488, 141723, 200768,
                                                                       88818, 89043, 144978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 205248, 0, 3,
                                                                       200768, 141933, 201048,
                                                                       89043, 89268, 145293,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 205668, 0, 3,
                                                                       201048, 142143, 201328,
                                                                       89268, 89493, 145608,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 206088, 0, 3,
                                                                       201328, 142353, 201608,
                                                                       89493, 89718, 145923,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 206508, 0, 3,
                                                                       201888, 142773, 202308,
                                                                       90168, 90483, 146238,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 207096, 0, 3,
                                                                       202308, 143088, 202728,
                                                                       90483, 90798, 146679,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 207684, 0, 3,
                                                                       202728, 143403, 203148,
                                                                       90798, 91113, 147120,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 208272, 0, 3,
                                                                       203148, 143718, 203568,
                                                                       91113, 91428, 147561,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 208860, 0, 3,
                                                                       203568, 144033, 203988,
                                                                       91428, 91743, 148002,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 209448, 0, 3,
                                                                       203988, 144348, 204408,
                                                                       91743, 92058, 148443,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 210036, 0, 3,
                                                                       204408, 144663, 204828,
                                                                       92058, 92373, 148884,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 210624, 0, 3,
                                                                       204828, 144978, 205248,
                                                                       92373, 92688, 149325,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 211212, 0, 3,
                                                                       205248, 145293, 205668,
                                                                       92688, 93003, 149766,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 211800, 0, 3,
                                                                       205668, 145608, 206088,
                                                                       93003, 93318, 150207,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 212388, 0, 3,
                                                                       206508, 146238, 207096,
                                                                       93948, 94368, 150648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 213172, 0, 3,
                                                                       207096, 146679, 207684,
                                                                       94368, 94788, 151236,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 213956, 0, 3,
                                                                       207684, 147120, 208272,
                                                                       94788, 95208, 151824,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 214740, 0, 3,
                                                                       208272, 147561, 208860,
                                                                       95208, 95628, 152412,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 215524, 0, 3,
                                                                       208860, 148002, 209448,
                                                                       95628, 96048, 153000,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 216308, 0, 3,
                                                                       209448, 148443, 210036,
                                                                       96048, 96468, 153588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 217092, 0, 3,
                                                                       210036, 148884, 210624,
                                                                       96468, 96888, 154176,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 217876, 0, 3,
                                                                       210624, 149325, 211212,
                                                                       96888, 97308, 154764,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 218660, 0, 3,
                                                                       211212, 149766, 211800,
                                                                       97308, 97728, 155352,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 219444, 0, 3,
                                                                       212388, 150648, 213172,
                                                                       98568, 99108, 155940,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 220452, 0, 3,
                                                                       213172, 151236, 213956,
                                                                       99108, 99648, 156696,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 221460, 0, 3,
                                                                       213956, 151824, 214740,
                                                                       99648, 100188, 157452,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 222468, 0, 3,
                                                                       214740, 152412, 215524,
                                                                       100188, 100728, 158208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 223476, 0, 3,
                                                                       215524, 153000, 216308,
                                                                       100728, 101268, 158964,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 224484, 0, 3,
                                                                       216308, 153588, 217092,
                                                                       101268, 101808, 159720,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 225492, 0, 3,
                                                                       217092, 154176, 217876,
                                                                       101808, 102348, 160476,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 226500, 0, 3,
                                                                       217876, 154764, 218660,
                                                                       102348, 102888, 161232,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 227508, 0, 3,
                                                                       219444, 155940, 220452,
                                                                       103968, 104643, 161988,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 228768, 0, 3,
                                                                       220452, 156696, 221460,
                                                                       104643, 105318, 162933,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 230028, 0, 3,
                                                                       221460, 157452, 222468,
                                                                       105318, 105993, 163878,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 231288, 0, 3,
                                                                       222468, 158208, 223476,
                                                                       105993, 106668, 164823,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 232548, 0, 3,
                                                                       223476, 158964, 224484,
                                                                       106668, 107343, 165768,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 233808, 0, 3,
                                                                       224484, 159720, 225492,
                                                                       107343, 108018, 166713,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 235068, 0, 3,
                                                                       225492, 160476, 226500,
                                                                       108018, 108693, 167658,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 236328, 0, 3,
                                                                       227508, 161988, 228768,
                                                                       110043, 110868, 168603,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 237868, 0, 3,
                                                                       228768, 162933, 230028,
                                                                       110868, 111693, 169758,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 239408, 0, 3,
                                                                       230028, 163878, 231288,
                                                                       111693, 112518, 170913,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 240948, 0, 3,
                                                                       231288, 164823, 232548,
                                                                       112518, 113343, 172068,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 242488, 0, 3,
                                                                       232548, 165768, 233808,
                                                                       113343, 114168, 173223,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 244028, 0, 3,
                                                                       233808, 166713, 235068,
                                                                       114168, 114993, 174378,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 245568, 0, 3,
                                                                       236328, 168603, 237868,
                                                                       116643, 117633, 175533,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 247416, 0, 3,
                                                                       237868, 169758, 239408,
                                                                       117633, 118623, 176919,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 249264, 0, 3,
                                                                       239408, 170913, 240948,
                                                                       118623, 119613, 178305,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 251112, 0, 3,
                                                                       240948, 172068, 242488,
                                                                       119613, 120603, 179691,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 252960, 0, 3,
                                                                       242488, 173223, 244028,
                                                                       120603, 121593, 181077,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 254808, 0, 3,
                                                                       245568, 175533, 247416,
                                                                       123573, 124743, 182463,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 256992, 0, 3,
                                                                       247416, 176919, 249264,
                                                                       124743, 125913, 184101,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 259176, 0, 3,
                                                                       249264, 178305, 251112,
                                                                       125913, 127083, 185739,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 261360, 0, 3,
                                                                       251112, 179691, 252960,
                                                                       127083, 128253, 187377,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 263544, 0, 3,
                                                                       254808, 182463, 256992,
                                                                       130593, 131958, 189015,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 266092, 0, 3,
                                                                       256992, 184101, 259176,
                                                                       131958, 133323, 190926,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 268640, 0, 3,
                                                                       259176, 185739, 261360,
                                                                       133323, 134688, 192837,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271188, 3, 137418,
                                                                       137439, 194804, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271224, 3, 137439,
                                                                       137460, 194832, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271260, 3, 137460,
                                                                       137481, 194860, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271296, 3, 137481,
                                                                       137502, 194888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271332, 3, 137502,
                                                                       137523, 194916, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271368, 3, 137523,
                                                                       137544, 194944, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271404, 3, 137544,
                                                                       137565, 194972, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271440, 3, 137565,
                                                                       137586, 195000, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271476, 3, 137586,
                                                                       137607, 195028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271512, 3, 137607,
                                                                       137628, 195056, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271548, 3, 137628,
                                                                       137649, 195084, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271584, 3, 137649,
                                                                       137670, 195112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 271620, 3, 137670,
                                                                       137691, 195140, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 271656, 0, 3,
                                                                       271188, 194804, 271224,
                                                                       137733, 137796, 195336,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 271764, 0, 3,
                                                                       271224, 194832, 271260,
                                                                       137796, 137859, 195420,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 271872, 0, 3,
                                                                       271260, 194860, 271296,
                                                                       137859, 137922, 195504,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 271980, 0, 3,
                                                                       271296, 194888, 271332,
                                                                       137922, 137985, 195588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272088, 0, 3,
                                                                       271332, 194916, 271368,
                                                                       137985, 138048, 195672,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272196, 0, 3,
                                                                       271368, 194944, 271404,
                                                                       138048, 138111, 195756,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272304, 0, 3,
                                                                       271404, 194972, 271440,
                                                                       138111, 138174, 195840,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272412, 0, 3,
                                                                       271440, 195000, 271476,
                                                                       138174, 138237, 195924,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272520, 0, 3,
                                                                       271476, 195028, 271512,
                                                                       138237, 138300, 196008,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272628, 0, 3,
                                                                       271512, 195056, 271548,
                                                                       138300, 138363, 196092,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272736, 0, 3,
                                                                       271548, 195084, 271584,
                                                                       138363, 138426, 196176,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 272844, 0, 3,
                                                                       271584, 195112, 271620,
                                                                       138426, 138489, 196260,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 272952, 0, 3,
                                                                       271656, 195336, 271764,
                                                                       138615, 138741, 196680,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 273168, 0, 3,
                                                                       271764, 195420, 271872,
                                                                       138741, 138867, 196848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 273384, 0, 3,
                                                                       271872, 195504, 271980,
                                                                       138867, 138993, 197016,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 273600, 0, 3,
                                                                       271980, 195588, 272088,
                                                                       138993, 139119, 197184,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 273816, 0, 3,
                                                                       272088, 195672, 272196,
                                                                       139119, 139245, 197352,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274032, 0, 3,
                                                                       272196, 195756, 272304,
                                                                       139245, 139371, 197520,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274248, 0, 3,
                                                                       272304, 195840, 272412,
                                                                       139371, 139497, 197688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274464, 0, 3,
                                                                       272412, 195924, 272520,
                                                                       139497, 139623, 197856,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274680, 0, 3,
                                                                       272520, 196008, 272628,
                                                                       139623, 139749, 198024,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 274896, 0, 3,
                                                                       272628, 196092, 272736,
                                                                       139749, 139875, 198192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 275112, 0, 3,
                                                                       272736, 196176, 272844,
                                                                       139875, 140001, 198360,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 275328, 0, 3,
                                                                       272952, 196680, 273168,
                                                                       140253, 140463, 199088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 275688, 0, 3,
                                                                       273168, 196848, 273384,
                                                                       140463, 140673, 199368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 276048, 0, 3,
                                                                       273384, 197016, 273600,
                                                                       140673, 140883, 199648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 276408, 0, 3,
                                                                       273600, 197184, 273816,
                                                                       140883, 141093, 199928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 276768, 0, 3,
                                                                       273816, 197352, 274032,
                                                                       141093, 141303, 200208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 277128, 0, 3,
                                                                       274032, 197520, 274248,
                                                                       141303, 141513, 200488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 277488, 0, 3,
                                                                       274248, 197688, 274464,
                                                                       141513, 141723, 200768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 277848, 0, 3,
                                                                       274464, 197856, 274680,
                                                                       141723, 141933, 201048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 278208, 0, 3,
                                                                       274680, 198024, 274896,
                                                                       141933, 142143, 201328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 278568, 0, 3,
                                                                       274896, 198192, 275112,
                                                                       142143, 142353, 201608,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 278928, 0, 3,
                                                                       275328, 199088, 275688,
                                                                       142773, 143088, 202728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 279468, 0, 3,
                                                                       275688, 199368, 276048,
                                                                       143088, 143403, 203148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 280008, 0, 3,
                                                                       276048, 199648, 276408,
                                                                       143403, 143718, 203568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 280548, 0, 3,
                                                                       276408, 199928, 276768,
                                                                       143718, 144033, 203988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 281088, 0, 3,
                                                                       276768, 200208, 277128,
                                                                       144033, 144348, 204408,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 281628, 0, 3,
                                                                       277128, 200488, 277488,
                                                                       144348, 144663, 204828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 282168, 0, 3,
                                                                       277488, 200768, 277848,
                                                                       144663, 144978, 205248,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 282708, 0, 3,
                                                                       277848, 201048, 278208,
                                                                       144978, 145293, 205668,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 283248, 0, 3,
                                                                       278208, 201328, 278568,
                                                                       145293, 145608, 206088,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 283788, 0, 3,
                                                                       278928, 202728, 279468,
                                                                       146238, 146679, 207684,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 284544, 0, 3,
                                                                       279468, 203148, 280008,
                                                                       146679, 147120, 208272,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 285300, 0, 3,
                                                                       280008, 203568, 280548,
                                                                       147120, 147561, 208860,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 286056, 0, 3,
                                                                       280548, 203988, 281088,
                                                                       147561, 148002, 209448,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 286812, 0, 3,
                                                                       281088, 204408, 281628,
                                                                       148002, 148443, 210036,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 287568, 0, 3,
                                                                       281628, 204828, 282168,
                                                                       148443, 148884, 210624,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 288324, 0, 3,
                                                                       282168, 205248, 282708,
                                                                       148884, 149325, 211212,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 289080, 0, 3,
                                                                       282708, 205668, 283248,
                                                                       149325, 149766, 211800,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 289836, 0, 3,
                                                                       283788, 207684, 284544,
                                                                       150648, 151236, 213956,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 290844, 0, 3,
                                                                       284544, 208272, 285300,
                                                                       151236, 151824, 214740,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 291852, 0, 3,
                                                                       285300, 208860, 286056,
                                                                       151824, 152412, 215524,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 292860, 0, 3,
                                                                       286056, 209448, 286812,
                                                                       152412, 153000, 216308,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 293868, 0, 3,
                                                                       286812, 210036, 287568,
                                                                       153000, 153588, 217092,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 294876, 0, 3,
                                                                       287568, 210624, 288324,
                                                                       153588, 154176, 217876,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 295884, 0, 3,
                                                                       288324, 211212, 289080,
                                                                       154176, 154764, 218660,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 296892, 0, 3,
                                                                       289836, 213956, 290844,
                                                                       155940, 156696, 221460,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 298188, 0, 3,
                                                                       290844, 214740, 291852,
                                                                       156696, 157452, 222468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 299484, 0, 3,
                                                                       291852, 215524, 292860,
                                                                       157452, 158208, 223476,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 300780, 0, 3,
                                                                       292860, 216308, 293868,
                                                                       158208, 158964, 224484,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 302076, 0, 3,
                                                                       293868, 217092, 294876,
                                                                       158964, 159720, 225492,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 303372, 0, 3,
                                                                       294876, 217876, 295884,
                                                                       159720, 160476, 226500,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 304668, 0, 3,
                                                                       296892, 221460, 298188,
                                                                       161988, 162933, 230028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 306288, 0, 3,
                                                                       298188, 222468, 299484,
                                                                       162933, 163878, 231288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 307908, 0, 3,
                                                                       299484, 223476, 300780,
                                                                       163878, 164823, 232548,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 309528, 0, 3,
                                                                       300780, 224484, 302076,
                                                                       164823, 165768, 233808,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 311148, 0, 3,
                                                                       302076, 225492, 303372,
                                                                       165768, 166713, 235068,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 312768, 0, 3,
                                                                       304668, 230028, 306288,
                                                                       168603, 169758, 239408,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 314748, 0, 3,
                                                                       306288, 231288, 307908,
                                                                       169758, 170913, 240948,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 316728, 0, 3,
                                                                       307908, 232548, 309528,
                                                                       170913, 172068, 242488,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 318708, 0, 3,
                                                                       309528, 233808, 311148,
                                                                       172068, 173223, 244028,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 320688, 0, 3,
                                                                       312768, 239408, 314748,
                                                                       175533, 176919, 249264,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 323064, 0, 3,
                                                                       314748, 240948, 316728,
                                                                       176919, 178305, 251112,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsk_three_center_electron_repulsion_0(buffer, 325440, 0, 3,
                                                                       316728, 242488, 318708,
                                                                       178305, 179691, 252960,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 327816, 0, 3,
                                                                       320688, 249264, 323064,
                                                                       182463, 184101, 259176,
                                                                       ncols, gamma, p, q);

                    compute_prim_osk_three_center_electron_repulsion_0(buffer, 330624, 0, 3,
                                                                       323064, 251112, 325440,
                                                                       184101, 185739, 261360,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsk_three_center_electron_repulsion_0(buffer, 333432, 0, 3,
                                                                       327816, 259176, 330624,
                                                                       189015, 190926, 268640,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336708, 3, 194748,
                                                                       194776, 271188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336753, 3, 194776,
                                                                       194804, 271224, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336798, 3, 194804,
                                                                       194832, 271260, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336843, 3, 194832,
                                                                       194860, 271296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336888, 3, 194860,
                                                                       194888, 271332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336933, 3, 194888,
                                                                       194916, 271368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 336978, 3, 194916,
                                                                       194944, 271404, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337023, 3, 194944,
                                                                       194972, 271440, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337068, 3, 194972,
                                                                       195000, 271476, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337113, 3, 195000,
                                                                       195028, 271512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337158, 3, 195028,
                                                                       195056, 271548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337203, 3, 195056,
                                                                       195084, 271584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 337248, 3, 195084,
                                                                       195112, 271620, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337293, 0, 3,
                                                                       336708, 271188, 336753,
                                                                       195168, 195252, 271656,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337428, 0, 3,
                                                                       336753, 271224, 336798,
                                                                       195252, 195336, 271764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337563, 0, 3,
                                                                       336798, 271260, 336843,
                                                                       195336, 195420, 271872,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337698, 0, 3,
                                                                       336843, 271296, 336888,
                                                                       195420, 195504, 271980,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337833, 0, 3,
                                                                       336888, 271332, 336933,
                                                                       195504, 195588, 272088,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 337968, 0, 3,
                                                                       336933, 271368, 336978,
                                                                       195588, 195672, 272196,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338103, 0, 3,
                                                                       336978, 271404, 337023,
                                                                       195672, 195756, 272304,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338238, 0, 3,
                                                                       337023, 271440, 337068,
                                                                       195756, 195840, 272412,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338373, 0, 3,
                                                                       337068, 271476, 337113,
                                                                       195840, 195924, 272520,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338508, 0, 3,
                                                                       337113, 271512, 337158,
                                                                       195924, 196008, 272628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338643, 0, 3,
                                                                       337158, 271548, 337203,
                                                                       196008, 196092, 272736,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 338778, 0, 3,
                                                                       337203, 271584, 337248,
                                                                       196092, 196176, 272844,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 338913, 0, 3,
                                                                       337293, 271656, 337428,
                                                                       196344, 196512, 272952,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 339183, 0, 3,
                                                                       337428, 271764, 337563,
                                                                       196512, 196680, 273168,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 339453, 0, 3,
                                                                       337563, 271872, 337698,
                                                                       196680, 196848, 273384,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 339723, 0, 3,
                                                                       337698, 271980, 337833,
                                                                       196848, 197016, 273600,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 339993, 0, 3,
                                                                       337833, 272088, 337968,
                                                                       197016, 197184, 273816,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 340263, 0, 3,
                                                                       337968, 272196, 338103,
                                                                       197184, 197352, 274032,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 340533, 0, 3,
                                                                       338103, 272304, 338238,
                                                                       197352, 197520, 274248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 340803, 0, 3,
                                                                       338238, 272412, 338373,
                                                                       197520, 197688, 274464,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 341073, 0, 3,
                                                                       338373, 272520, 338508,
                                                                       197688, 197856, 274680,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 341343, 0, 3,
                                                                       338508, 272628, 338643,
                                                                       197856, 198024, 274896,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 341613, 0, 3,
                                                                       338643, 272736, 338778,
                                                                       198024, 198192, 275112,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 341883, 0, 3,
                                                                       338913, 272952, 339183,
                                                                       198528, 198808, 275328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 342333, 0, 3,
                                                                       339183, 273168, 339453,
                                                                       198808, 199088, 275688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 342783, 0, 3,
                                                                       339453, 273384, 339723,
                                                                       199088, 199368, 276048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 343233, 0, 3,
                                                                       339723, 273600, 339993,
                                                                       199368, 199648, 276408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 343683, 0, 3,
                                                                       339993, 273816, 340263,
                                                                       199648, 199928, 276768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 344133, 0, 3,
                                                                       340263, 274032, 340533,
                                                                       199928, 200208, 277128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 344583, 0, 3,
                                                                       340533, 274248, 340803,
                                                                       200208, 200488, 277488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 345033, 0, 3,
                                                                       340803, 274464, 341073,
                                                                       200488, 200768, 277848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 345483, 0, 3,
                                                                       341073, 274680, 341343,
                                                                       200768, 201048, 278208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 345933, 0, 3,
                                                                       341343, 274896, 341613,
                                                                       201048, 201328, 278568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 346383, 0, 3,
                                                                       341883, 275328, 342333,
                                                                       201888, 202308, 278928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 347058, 0, 3,
                                                                       342333, 275688, 342783,
                                                                       202308, 202728, 279468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 347733, 0, 3,
                                                                       342783, 276048, 343233,
                                                                       202728, 203148, 280008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 348408, 0, 3,
                                                                       343233, 276408, 343683,
                                                                       203148, 203568, 280548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 349083, 0, 3,
                                                                       343683, 276768, 344133,
                                                                       203568, 203988, 281088,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 349758, 0, 3,
                                                                       344133, 277128, 344583,
                                                                       203988, 204408, 281628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 350433, 0, 3,
                                                                       344583, 277488, 345033,
                                                                       204408, 204828, 282168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 351108, 0, 3,
                                                                       345033, 277848, 345483,
                                                                       204828, 205248, 282708,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 351783, 0, 3,
                                                                       345483, 278208, 345933,
                                                                       205248, 205668, 283248,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 352458, 0, 3,
                                                                       346383, 278928, 347058,
                                                                       206508, 207096, 283788,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 353403, 0, 3,
                                                                       347058, 279468, 347733,
                                                                       207096, 207684, 284544,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 354348, 0, 3,
                                                                       347733, 280008, 348408,
                                                                       207684, 208272, 285300,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 355293, 0, 3,
                                                                       348408, 280548, 349083,
                                                                       208272, 208860, 286056,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 356238, 0, 3,
                                                                       349083, 281088, 349758,
                                                                       208860, 209448, 286812,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 357183, 0, 3,
                                                                       349758, 281628, 350433,
                                                                       209448, 210036, 287568,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 358128, 0, 3,
                                                                       350433, 282168, 351108,
                                                                       210036, 210624, 288324,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 359073, 0, 3,
                                                                       351108, 282708, 351783,
                                                                       210624, 211212, 289080,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 360018, 0, 3,
                                                                       352458, 283788, 353403,
                                                                       212388, 213172, 289836,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 361278, 0, 3,
                                                                       353403, 284544, 354348,
                                                                       213172, 213956, 290844,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 362538, 0, 3,
                                                                       354348, 285300, 355293,
                                                                       213956, 214740, 291852,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 363798, 0, 3,
                                                                       355293, 286056, 356238,
                                                                       214740, 215524, 292860,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 365058, 0, 3,
                                                                       356238, 286812, 357183,
                                                                       215524, 216308, 293868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 366318, 0, 3,
                                                                       357183, 287568, 358128,
                                                                       216308, 217092, 294876,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 367578, 0, 3,
                                                                       358128, 288324, 359073,
                                                                       217092, 217876, 295884,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 368838, 0, 3,
                                                                       360018, 289836, 361278,
                                                                       219444, 220452, 296892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 370458, 0, 3,
                                                                       361278, 290844, 362538,
                                                                       220452, 221460, 298188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 372078, 0, 3,
                                                                       362538, 291852, 363798,
                                                                       221460, 222468, 299484,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 373698, 0, 3,
                                                                       363798, 292860, 365058,
                                                                       222468, 223476, 300780,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 375318, 0, 3,
                                                                       365058, 293868, 366318,
                                                                       223476, 224484, 302076,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 376938, 0, 3,
                                                                       366318, 294876, 367578,
                                                                       224484, 225492, 303372,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 378558, 0, 3,
                                                                       368838, 296892, 370458,
                                                                       227508, 228768, 304668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 380583, 0, 3,
                                                                       370458, 298188, 372078,
                                                                       228768, 230028, 306288,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 382608, 0, 3,
                                                                       372078, 299484, 373698,
                                                                       230028, 231288, 307908,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 384633, 0, 3,
                                                                       373698, 300780, 375318,
                                                                       231288, 232548, 309528,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 386658, 0, 3,
                                                                       375318, 302076, 376938,
                                                                       232548, 233808, 311148,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 388683, 0, 3,
                                                                       378558, 304668, 380583,
                                                                       236328, 237868, 312768,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 391158, 0, 3,
                                                                       380583, 306288, 382608,
                                                                       237868, 239408, 314748,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 393633, 0, 3,
                                                                       382608, 307908, 384633,
                                                                       239408, 240948, 316728,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 396108, 0, 3,
                                                                       384633, 309528, 386658,
                                                                       240948, 242488, 318708,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 398583, 0, 3,
                                                                       388683, 312768, 391158,
                                                                       245568, 247416, 320688,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 401553, 0, 3,
                                                                       391158, 314748, 393633,
                                                                       247416, 249264, 323064,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsl_three_center_electron_repulsion_0(buffer, 404523, 0, 3,
                                                                       393633, 316728, 396108,
                                                                       249264, 251112, 325440,
                                                                       ncols, gamma, p, q);

                    compute_prim_osl_three_center_electron_repulsion_0(buffer, 407493, 0, 3,
                                                                       398583, 320688, 401553,
                                                                       254808, 256992, 327816,
                                                                       ncols, gamma, p, q);

                    compute_prim_osl_three_center_electron_repulsion_0(buffer, 411003, 0, 3,
                                                                       401553, 323064, 404523,
                                                                       256992, 259176, 330624,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsl_three_center_electron_repulsion_0(buffer, 414513, 0, 3,
                                                                       407493, 327816, 411003,
                                                                       263544, 266092, 333432,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 418608, 360018, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 420344, 368838, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 422576, 378558, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 425366, 388683, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 428776, 398583, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 432868, 407493, 3510, ncols);

                    simdfunc::contract_primitives(buffer, 437704, 414513, 4095, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 419868, 418608, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 421964, 420344, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 424601, 422576, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 427841, 425366, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 431746, 428776, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 436378, 432868, 78, 1, nmax);

        simdtrf::transform_l_inner(buffer, 441799, 437704, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 443346, 419868, 421964, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 444774, 421964, 424601, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 446610, 424601, 427841, 17,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 448905, 427841, 431746, 17,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 451710, 431746, 436378, 17,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 455076, 436378, 441799, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 459054, 443346, 444774, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 461910, 444774, 446610, 17,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 465582, 446610, 448905, 17,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 470172, 448905, 451710, 17,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 475782, 451710, 455076, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 482514, 459054, 461910, 17,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 487274, 461910, 465582, 17,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 493394, 465582, 470172, 17,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 501044, 470172, 475782, 17,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 510394, 482514, 487274, 17,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 517534, 487274, 493394, 17,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 526714, 493394, 501044, 17,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 538189, 510394, 517534, 17,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 548185, 517534, 526714, 17,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 561037, 538189, 548185, 17,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 574365, 561037, 28, 17, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 574365, 221, nmax);
    }

    for (size_t m = 0; m < 2873; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
