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


#include "SimdThreeCenterElectronRepulsionRsRecIII.hpp"

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
compute_rs_iii_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_iii_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 627320, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 4394 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 627320, 389488, 31535, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 18,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 26, 3, 18,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9108, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9111, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9114, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9117, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9120, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9123, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9126, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9129, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9132, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9135, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9138, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9141, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9144, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9147, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9150, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9153, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9156, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9159, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9162, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9165, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9168, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9171, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9174, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9177, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9180, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9183, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9186, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9189, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9192, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9195, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9198, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9201, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9204, 3, 44,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 9207, 3, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9210, 3, 9, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9219, 3, 10, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9228, 3, 11, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9237, 3, 12, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9246, 3, 13, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9255, 3, 14, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9264, 3, 15, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9273, 3, 16, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9282, 3, 17, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9291, 3, 18, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9300, 3, 19, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9309, 3, 20, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9318, 3, 21, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9327, 3, 22, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9336, 3, 23, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9345, 3, 24, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9354, 3, 29, 106,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9363, 3, 30, 109,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9372, 3, 31, 112,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9381, 3, 32, 115,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9390, 3, 33, 118,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9399, 3, 34, 121,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9408, 3, 35, 124,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9417, 3, 36, 127,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9426, 3, 37, 130,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9435, 3, 38, 133,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9444, 3, 39, 136,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9453, 3, 40, 139,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9462, 3, 41, 142,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9471, 3, 42, 145,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9480, 3, 43, 148,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 9489, 3, 44, 151,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9498, 3, 52, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9516, 3, 55, 172,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9534, 3, 58, 178,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9552, 3, 61, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9570, 3, 64, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9588, 3, 67, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9606, 3, 70, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9624, 3, 73, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9642, 3, 76, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9660, 3, 79, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9678, 3, 82, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9696, 3, 85, 232,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9714, 3, 88, 238,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9732, 3, 91, 244,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9750, 3, 94, 250,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9768, 3, 106, 268,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9786, 3, 109, 274,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9804, 3, 112, 280,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9822, 3, 115, 286,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9840, 3, 118, 292,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9858, 3, 121, 298,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9876, 3, 124, 304,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9894, 3, 127, 310,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9912, 3, 130, 316,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9930, 3, 133, 322,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9948, 3, 136, 328,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9966, 3, 139, 334,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 9984, 3, 142, 340,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10002, 3, 145,
                                                                       346, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 10020, 3, 148,
                                                                       352, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10038, 3, 166,
                                                                       378, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10068, 3, 172,
                                                                       388, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10098, 3, 178,
                                                                       398, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10128, 3, 184,
                                                                       408, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10158, 3, 190,
                                                                       418, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10188, 3, 196,
                                                                       428, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10218, 3, 202,
                                                                       438, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10248, 3, 208,
                                                                       448, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10278, 3, 214,
                                                                       458, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10308, 3, 220,
                                                                       468, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10338, 3, 226,
                                                                       478, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10368, 3, 232,
                                                                       488, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10398, 3, 238,
                                                                       498, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10428, 3, 244,
                                                                       508, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10458, 3, 268,
                                                                       538, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10488, 3, 274,
                                                                       548, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10518, 3, 280,
                                                                       558, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10548, 3, 286,
                                                                       568, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10578, 3, 292,
                                                                       578, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10608, 3, 298,
                                                                       588, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10638, 3, 304,
                                                                       598, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10668, 3, 310,
                                                                       608, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10698, 3, 316,
                                                                       618, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10728, 3, 322,
                                                                       628, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10758, 3, 328,
                                                                       638, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10788, 3, 334,
                                                                       648, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10818, 3, 340,
                                                                       658, ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 10848, 3, 346,
                                                                       668, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10878, 3, 378,
                                                                       708, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10923, 3, 388,
                                                                       723, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 10968, 3, 398,
                                                                       738, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11013, 3, 408,
                                                                       753, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11058, 3, 418,
                                                                       768, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11103, 3, 428,
                                                                       783, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11148, 3, 438,
                                                                       798, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11193, 3, 448,
                                                                       813, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11238, 3, 458,
                                                                       828, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11283, 3, 468,
                                                                       843, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11328, 3, 478,
                                                                       858, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11373, 3, 488,
                                                                       873, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11418, 3, 498,
                                                                       888, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11463, 3, 538,
                                                                       933, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11508, 3, 548,
                                                                       948, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11553, 3, 558,
                                                                       963, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11598, 3, 568,
                                                                       978, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11643, 3, 578,
                                                                       993, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11688, 3, 588,
                                                                       1008, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11733, 3, 598,
                                                                       1023, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11778, 3, 608,
                                                                       1038, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11823, 3, 618,
                                                                       1053, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11868, 3, 628,
                                                                       1068, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11913, 3, 638,
                                                                       1083, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 11958, 3, 648,
                                                                       1098, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 12003, 3, 658,
                                                                       1113, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12048, 3, 708,
                                                                       1170, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12111, 3, 723,
                                                                       1191, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12174, 3, 738,
                                                                       1212, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12237, 3, 753,
                                                                       1233, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12300, 3, 768,
                                                                       1254, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12363, 3, 783,
                                                                       1275, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12426, 3, 798,
                                                                       1296, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12489, 3, 813,
                                                                       1317, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12552, 3, 828,
                                                                       1338, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12615, 3, 843,
                                                                       1359, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12678, 3, 858,
                                                                       1380, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12741, 3, 873,
                                                                       1401, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12804, 3, 933,
                                                                       1464, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12867, 3, 948,
                                                                       1485, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12930, 3, 963,
                                                                       1506, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 12993, 3, 978,
                                                                       1527, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13056, 3, 993,
                                                                       1548, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13119, 3, 1008,
                                                                       1569, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13182, 3, 1023,
                                                                       1590, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13245, 3, 1038,
                                                                       1611, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13308, 3, 1053,
                                                                       1632, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13371, 3, 1068,
                                                                       1653, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13434, 3, 1083,
                                                                       1674, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 13497, 3, 1098,
                                                                       1695, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13560, 3, 1170,
                                                                       1772, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13644, 3, 1191,
                                                                       1800, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13728, 3, 1212,
                                                                       1828, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13812, 3, 1233,
                                                                       1856, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13896, 3, 1254,
                                                                       1884, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 13980, 3, 1275,
                                                                       1912, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14064, 3, 1296,
                                                                       1940, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14148, 3, 1317,
                                                                       1968, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14232, 3, 1338,
                                                                       1996, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14316, 3, 1359,
                                                                       2024, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14400, 3, 1380,
                                                                       2052, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14484, 3, 1464,
                                                                       2136, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14568, 3, 1485,
                                                                       2164, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14652, 3, 1506,
                                                                       2192, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14736, 3, 1527,
                                                                       2220, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14820, 3, 1548,
                                                                       2248, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14904, 3, 1569,
                                                                       2276, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 14988, 3, 1590,
                                                                       2304, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15072, 3, 1611,
                                                                       2332, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15156, 3, 1632,
                                                                       2360, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15240, 3, 1653,
                                                                       2388, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 15324, 3, 1674,
                                                                       2416, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15408, 3, 1772,
                                                                       2516, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15516, 3, 1800,
                                                                       2552, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15624, 3, 1828,
                                                                       2588, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15732, 3, 1856,
                                                                       2624, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15840, 3, 1884,
                                                                       2660, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 15948, 3, 1912,
                                                                       2696, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16056, 3, 1940,
                                                                       2732, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16164, 3, 1968,
                                                                       2768, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16272, 3, 1996,
                                                                       2804, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16380, 3, 2024,
                                                                       2840, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16488, 3, 2136,
                                                                       2948, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16596, 3, 2164,
                                                                       2984, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16704, 3, 2192,
                                                                       3020, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16812, 3, 2220,
                                                                       3056, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 16920, 3, 2248,
                                                                       3092, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17028, 3, 2276,
                                                                       3128, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17136, 3, 2304,
                                                                       3164, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17244, 3, 2332,
                                                                       3200, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17352, 3, 2360,
                                                                       3236, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 17460, 3, 2388,
                                                                       3272, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17568, 3, 2516,
                                                                       3398, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17703, 3, 2552,
                                                                       3443, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17838, 3, 2588,
                                                                       3488, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 17973, 3, 2624,
                                                                       3533, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18108, 3, 2660,
                                                                       3578, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18243, 3, 2696,
                                                                       3623, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18378, 3, 2732,
                                                                       3668, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18513, 3, 2768,
                                                                       3713, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18648, 3, 2804,
                                                                       3758, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18783, 3, 2948,
                                                                       3893, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 18918, 3, 2984,
                                                                       3938, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19053, 3, 3020,
                                                                       3983, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19188, 3, 3056,
                                                                       4028, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19323, 3, 3092,
                                                                       4073, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19458, 3, 3128,
                                                                       4118, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19593, 3, 3164,
                                                                       4163, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19728, 3, 3200,
                                                                       4208, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 19863, 3, 3236,
                                                                       4253, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 19998, 3, 3398,
                                                                       4408, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20163, 3, 3443,
                                                                       4463, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20328, 3, 3488,
                                                                       4518, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20493, 3, 3533,
                                                                       4573, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20658, 3, 3578,
                                                                       4628, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20823, 3, 3623,
                                                                       4683, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 20988, 3, 3668,
                                                                       4738, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21153, 3, 3713,
                                                                       4793, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21318, 3, 3893,
                                                                       4958, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21483, 3, 3938,
                                                                       5013, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21648, 3, 3983,
                                                                       5068, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21813, 3, 4028,
                                                                       5123, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 21978, 3, 4073,
                                                                       5178, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22143, 3, 4118,
                                                                       5233, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22308, 3, 4163,
                                                                       5288, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 22473, 3, 4208,
                                                                       5343, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22638, 3, 4408,
                                                                       5530, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 22836, 3, 4463,
                                                                       5596, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23034, 3, 4518,
                                                                       5662, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23232, 3, 4573,
                                                                       5728, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23430, 3, 4628,
                                                                       5794, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23628, 3, 4683,
                                                                       5860, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 23826, 3, 4738,
                                                                       5926, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24024, 3, 4958,
                                                                       6124, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24222, 3, 5013,
                                                                       6190, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24420, 3, 5068,
                                                                       6256, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24618, 3, 5123,
                                                                       6322, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 24816, 3, 5178,
                                                                       6388, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 25014, 3, 5233,
                                                                       6454, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 25212, 3, 5288,
                                                                       6520, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25410, 3, 5530,
                                                                       6742, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25644, 3, 5596,
                                                                       6820, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 25878, 3, 5662,
                                                                       6898, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26112, 3, 5728,
                                                                       6976, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26346, 3, 5794,
                                                                       7054, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26580, 3, 5860,
                                                                       7132, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 26814, 3, 6124,
                                                                       7366, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27048, 3, 6190,
                                                                       7444, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27282, 3, 6256,
                                                                       7522, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27516, 3, 6322,
                                                                       7600, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27750, 3, 6388,
                                                                       7678, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 27984, 3, 6454,
                                                                       7756, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 28218, 3, 6742,
                                                                       8016, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 28491, 3, 6820,
                                                                       8107, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 28764, 3, 6898,
                                                                       8198, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 29037, 3, 6976,
                                                                       8289, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 29310, 3, 7054,
                                                                       8380, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 29583, 3, 7366,
                                                                       8653, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 29856, 3, 7444,
                                                                       8744, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 30129, 3, 7522,
                                                                       8835, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 30402, 3, 7600,
                                                                       8926, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 30675, 3, 7678,
                                                                       9017, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30948, 3, 7, 8,
                                                                       9108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30954, 3, 8, 9,
                                                                       9111, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30960, 3, 9, 10,
                                                                       9114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30966, 3, 10, 11,
                                                                       9117, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30972, 3, 11, 12,
                                                                       9120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30978, 3, 12, 13,
                                                                       9123, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30984, 3, 13, 14,
                                                                       9126, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30990, 3, 14, 15,
                                                                       9129, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 30996, 3, 15, 16,
                                                                       9132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31002, 3, 16, 17,
                                                                       9135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31008, 3, 17, 18,
                                                                       9138, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31014, 3, 18, 19,
                                                                       9141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31020, 3, 19, 20,
                                                                       9144, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31026, 3, 20, 21,
                                                                       9147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31032, 3, 21, 22,
                                                                       9150, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31038, 3, 22, 23,
                                                                       9153, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31044, 3, 23, 24,
                                                                       9156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31050, 3, 27, 28,
                                                                       9159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31056, 3, 28, 29,
                                                                       9162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31062, 3, 29, 30,
                                                                       9165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31068, 3, 30, 31,
                                                                       9168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31074, 3, 31, 32,
                                                                       9171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31080, 3, 32, 33,
                                                                       9174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31086, 3, 33, 34,
                                                                       9177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31092, 3, 34, 35,
                                                                       9180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31098, 3, 35, 36,
                                                                       9183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31104, 3, 36, 37,
                                                                       9186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31110, 3, 37, 38,
                                                                       9189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31116, 3, 38, 39,
                                                                       9192, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31122, 3, 39, 40,
                                                                       9195, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31128, 3, 40, 41,
                                                                       9198, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31134, 3, 41, 42,
                                                                       9201, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31140, 3, 42, 43,
                                                                       9204, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 31146, 3, 43, 44,
                                                                       9207, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31152, 0, 3,
                                                                       30948, 9108, 30954, 9210,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31170, 0, 3,
                                                                       30954, 9111, 30960, 9219,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31188, 0, 3,
                                                                       30960, 9114, 30966, 9228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31206, 0, 3,
                                                                       30966, 9117, 30972, 9237,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31224, 0, 3,
                                                                       30972, 9120, 30978, 9246,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31242, 0, 3,
                                                                       30978, 9123, 30984, 9255,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31260, 0, 3,
                                                                       30984, 9126, 30990, 9264,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31278, 0, 3,
                                                                       30990, 9129, 30996, 9273,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31296, 0, 3,
                                                                       30996, 9132, 31002, 9282,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31314, 0, 3,
                                                                       31002, 9135, 31008, 9291,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31332, 0, 3,
                                                                       31008, 9138, 31014, 9300,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31350, 0, 3,
                                                                       31014, 9141, 31020, 9309,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31368, 0, 3,
                                                                       31020, 9144, 31026, 9318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31386, 0, 3,
                                                                       31026, 9147, 31032, 9327,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31404, 0, 3,
                                                                       31032, 9150, 31038, 9336,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31422, 0, 3,
                                                                       31038, 9153, 31044, 9345,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31440, 0, 3,
                                                                       31050, 9159, 31056, 9354,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31458, 0, 3,
                                                                       31056, 9162, 31062, 9363,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31476, 0, 3,
                                                                       31062, 9165, 31068, 9372,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31494, 0, 3,
                                                                       31068, 9168, 31074, 9381,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31512, 0, 3,
                                                                       31074, 9171, 31080, 9390,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31530, 0, 3,
                                                                       31080, 9174, 31086, 9399,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31548, 0, 3,
                                                                       31086, 9177, 31092, 9408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31566, 0, 3,
                                                                       31092, 9180, 31098, 9417,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31584, 0, 3,
                                                                       31098, 9183, 31104, 9426,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31602, 0, 3,
                                                                       31104, 9186, 31110, 9435,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31620, 0, 3,
                                                                       31110, 9189, 31116, 9444,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31638, 0, 3,
                                                                       31116, 9192, 31122, 9453,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31656, 0, 3,
                                                                       31122, 9195, 31128, 9462,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31674, 0, 3,
                                                                       31128, 9198, 31134, 9471,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31692, 0, 3,
                                                                       31134, 9201, 31140, 9480,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 31710, 0, 3,
                                                                       31140, 9204, 31146, 9489,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 31728, 0, 3,
                                                                       31152, 9210, 31170, 154,
                                                                       160, 9498, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 31764, 0, 3,
                                                                       31170, 9219, 31188, 160,
                                                                       166, 9516, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 31800, 0, 3,
                                                                       31188, 9228, 31206, 166,
                                                                       172, 9534, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 31836, 0, 3,
                                                                       31206, 9237, 31224, 172,
                                                                       178, 9552, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 31872, 0, 3,
                                                                       31224, 9246, 31242, 178,
                                                                       184, 9570, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 31908, 0, 3,
                                                                       31242, 9255, 31260, 184,
                                                                       190, 9588, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 31944, 0, 3,
                                                                       31260, 9264, 31278, 190,
                                                                       196, 9606, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 31980, 0, 3,
                                                                       31278, 9273, 31296, 196,
                                                                       202, 9624, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32016, 0, 3,
                                                                       31296, 9282, 31314, 202,
                                                                       208, 9642, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32052, 0, 3,
                                                                       31314, 9291, 31332, 208,
                                                                       214, 9660, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32088, 0, 3,
                                                                       31332, 9300, 31350, 214,
                                                                       220, 9678, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32124, 0, 3,
                                                                       31350, 9309, 31368, 220,
                                                                       226, 9696, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32160, 0, 3,
                                                                       31368, 9318, 31386, 226,
                                                                       232, 9714, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32196, 0, 3,
                                                                       31386, 9327, 31404, 232,
                                                                       238, 9732, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32232, 0, 3,
                                                                       31404, 9336, 31422, 238,
                                                                       244, 9750, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32268, 0, 3,
                                                                       31440, 9354, 31458, 256,
                                                                       262, 9768, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32304, 0, 3,
                                                                       31458, 9363, 31476, 262,
                                                                       268, 9786, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32340, 0, 3,
                                                                       31476, 9372, 31494, 268,
                                                                       274, 9804, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32376, 0, 3,
                                                                       31494, 9381, 31512, 274,
                                                                       280, 9822, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32412, 0, 3,
                                                                       31512, 9390, 31530, 280,
                                                                       286, 9840, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32448, 0, 3,
                                                                       31530, 9399, 31548, 286,
                                                                       292, 9858, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32484, 0, 3,
                                                                       31548, 9408, 31566, 292,
                                                                       298, 9876, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32520, 0, 3,
                                                                       31566, 9417, 31584, 298,
                                                                       304, 9894, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32556, 0, 3,
                                                                       31584, 9426, 31602, 304,
                                                                       310, 9912, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32592, 0, 3,
                                                                       31602, 9435, 31620, 310,
                                                                       316, 9930, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32628, 0, 3,
                                                                       31620, 9444, 31638, 316,
                                                                       322, 9948, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32664, 0, 3,
                                                                       31638, 9453, 31656, 322,
                                                                       328, 9966, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32700, 0, 3,
                                                                       31656, 9462, 31674, 328,
                                                                       334, 9984, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32736, 0, 3,
                                                                       31674, 9471, 31692, 334,
                                                                       340, 10002, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 32772, 0, 3,
                                                                       31692, 9480, 31710, 340,
                                                                       346, 10020, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 32808, 0, 3,
                                                                       31728, 9498, 31764, 358,
                                                                       368, 10038, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 32868, 0, 3,
                                                                       31764, 9516, 31800, 368,
                                                                       378, 10068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 32928, 0, 3,
                                                                       31800, 9534, 31836, 378,
                                                                       388, 10098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 32988, 0, 3,
                                                                       31836, 9552, 31872, 388,
                                                                       398, 10128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33048, 0, 3,
                                                                       31872, 9570, 31908, 398,
                                                                       408, 10158, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33108, 0, 3,
                                                                       31908, 9588, 31944, 408,
                                                                       418, 10188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33168, 0, 3,
                                                                       31944, 9606, 31980, 418,
                                                                       428, 10218, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33228, 0, 3,
                                                                       31980, 9624, 32016, 428,
                                                                       438, 10248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33288, 0, 3,
                                                                       32016, 9642, 32052, 438,
                                                                       448, 10278, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33348, 0, 3,
                                                                       32052, 9660, 32088, 448,
                                                                       458, 10308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33408, 0, 3,
                                                                       32088, 9678, 32124, 458,
                                                                       468, 10338, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33468, 0, 3,
                                                                       32124, 9696, 32160, 468,
                                                                       478, 10368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33528, 0, 3,
                                                                       32160, 9714, 32196, 478,
                                                                       488, 10398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33588, 0, 3,
                                                                       32196, 9732, 32232, 488,
                                                                       498, 10428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33648, 0, 3,
                                                                       32268, 9768, 32304, 518,
                                                                       528, 10458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33708, 0, 3,
                                                                       32304, 9786, 32340, 528,
                                                                       538, 10488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33768, 0, 3,
                                                                       32340, 9804, 32376, 538,
                                                                       548, 10518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33828, 0, 3,
                                                                       32376, 9822, 32412, 548,
                                                                       558, 10548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33888, 0, 3,
                                                                       32412, 9840, 32448, 558,
                                                                       568, 10578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 33948, 0, 3,
                                                                       32448, 9858, 32484, 568,
                                                                       578, 10608, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 34008, 0, 3,
                                                                       32484, 9876, 32520, 578,
                                                                       588, 10638, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 34068, 0, 3,
                                                                       32520, 9894, 32556, 588,
                                                                       598, 10668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 34128, 0, 3,
                                                                       32556, 9912, 32592, 598,
                                                                       608, 10698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 34188, 0, 3,
                                                                       32592, 9930, 32628, 608,
                                                                       618, 10728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 34248, 0, 3,
                                                                       32628, 9948, 32664, 618,
                                                                       628, 10758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 34308, 0, 3,
                                                                       32664, 9966, 32700, 628,
                                                                       638, 10788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 34368, 0, 3,
                                                                       32700, 9984, 32736, 638,
                                                                       648, 10818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 34428, 0, 3,
                                                                       32736, 10002, 32772, 648,
                                                                       658, 10848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 34488, 0, 3,
                                                                       32808, 10038, 32868, 678,
                                                                       693, 10878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 34578, 0, 3,
                                                                       32868, 10068, 32928, 693,
                                                                       708, 10923, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 34668, 0, 3,
                                                                       32928, 10098, 32988, 708,
                                                                       723, 10968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 34758, 0, 3,
                                                                       32988, 10128, 33048, 723,
                                                                       738, 11013, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 34848, 0, 3,
                                                                       33048, 10158, 33108, 738,
                                                                       753, 11058, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 34938, 0, 3,
                                                                       33108, 10188, 33168, 753,
                                                                       768, 11103, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35028, 0, 3,
                                                                       33168, 10218, 33228, 768,
                                                                       783, 11148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35118, 0, 3,
                                                                       33228, 10248, 33288, 783,
                                                                       798, 11193, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35208, 0, 3,
                                                                       33288, 10278, 33348, 798,
                                                                       813, 11238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35298, 0, 3,
                                                                       33348, 10308, 33408, 813,
                                                                       828, 11283, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35388, 0, 3,
                                                                       33408, 10338, 33468, 828,
                                                                       843, 11328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35478, 0, 3,
                                                                       33468, 10368, 33528, 843,
                                                                       858, 11373, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35568, 0, 3,
                                                                       33528, 10398, 33588, 858,
                                                                       873, 11418, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35658, 0, 3,
                                                                       33648, 10458, 33708, 903,
                                                                       918, 11463, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35748, 0, 3,
                                                                       33708, 10488, 33768, 918,
                                                                       933, 11508, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35838, 0, 3,
                                                                       33768, 10518, 33828, 933,
                                                                       948, 11553, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 35928, 0, 3,
                                                                       33828, 10548, 33888, 948,
                                                                       963, 11598, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36018, 0, 3,
                                                                       33888, 10578, 33948, 963,
                                                                       978, 11643, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36108, 0, 3,
                                                                       33948, 10608, 34008, 978,
                                                                       993, 11688, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36198, 0, 3,
                                                                       34008, 10638, 34068, 993,
                                                                       1008, 11733, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36288, 0, 3,
                                                                       34068, 10668, 34128, 1008,
                                                                       1023, 11778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36378, 0, 3,
                                                                       34128, 10698, 34188, 1023,
                                                                       1038, 11823, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36468, 0, 3,
                                                                       34188, 10728, 34248, 1038,
                                                                       1053, 11868, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36558, 0, 3,
                                                                       34248, 10758, 34308, 1053,
                                                                       1068, 11913, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36648, 0, 3,
                                                                       34308, 10788, 34368, 1068,
                                                                       1083, 11958, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 36738, 0, 3,
                                                                       34368, 10818, 34428, 1083,
                                                                       1098, 12003, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36828, 0, 3,
                                                                       34488, 10878, 34578, 1128,
                                                                       1149, 12048, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 36954, 0, 3,
                                                                       34578, 10923, 34668, 1149,
                                                                       1170, 12111, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 37080, 0, 3,
                                                                       34668, 10968, 34758, 1170,
                                                                       1191, 12174, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 37206, 0, 3,
                                                                       34758, 11013, 34848, 1191,
                                                                       1212, 12237, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 37332, 0, 3,
                                                                       34848, 11058, 34938, 1212,
                                                                       1233, 12300, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 37458, 0, 3,
                                                                       34938, 11103, 35028, 1233,
                                                                       1254, 12363, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 37584, 0, 3,
                                                                       35028, 11148, 35118, 1254,
                                                                       1275, 12426, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 37710, 0, 3,
                                                                       35118, 11193, 35208, 1275,
                                                                       1296, 12489, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 37836, 0, 3,
                                                                       35208, 11238, 35298, 1296,
                                                                       1317, 12552, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 37962, 0, 3,
                                                                       35298, 11283, 35388, 1317,
                                                                       1338, 12615, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 38088, 0, 3,
                                                                       35388, 11328, 35478, 1338,
                                                                       1359, 12678, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 38214, 0, 3,
                                                                       35478, 11373, 35568, 1359,
                                                                       1380, 12741, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 38340, 0, 3,
                                                                       35658, 11463, 35748, 1422,
                                                                       1443, 12804, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 38466, 0, 3,
                                                                       35748, 11508, 35838, 1443,
                                                                       1464, 12867, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 38592, 0, 3,
                                                                       35838, 11553, 35928, 1464,
                                                                       1485, 12930, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 38718, 0, 3,
                                                                       35928, 11598, 36018, 1485,
                                                                       1506, 12993, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 38844, 0, 3,
                                                                       36018, 11643, 36108, 1506,
                                                                       1527, 13056, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 38970, 0, 3,
                                                                       36108, 11688, 36198, 1527,
                                                                       1548, 13119, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 39096, 0, 3,
                                                                       36198, 11733, 36288, 1548,
                                                                       1569, 13182, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 39222, 0, 3,
                                                                       36288, 11778, 36378, 1569,
                                                                       1590, 13245, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 39348, 0, 3,
                                                                       36378, 11823, 36468, 1590,
                                                                       1611, 13308, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 39474, 0, 3,
                                                                       36468, 11868, 36558, 1611,
                                                                       1632, 13371, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 39600, 0, 3,
                                                                       36558, 11913, 36648, 1632,
                                                                       1653, 13434, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 39726, 0, 3,
                                                                       36648, 11958, 36738, 1653,
                                                                       1674, 13497, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 39852, 0, 3,
                                                                       36828, 12048, 36954, 1716,
                                                                       1744, 13560, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 40020, 0, 3,
                                                                       36954, 12111, 37080, 1744,
                                                                       1772, 13644, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 40188, 0, 3,
                                                                       37080, 12174, 37206, 1772,
                                                                       1800, 13728, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 40356, 0, 3,
                                                                       37206, 12237, 37332, 1800,
                                                                       1828, 13812, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 40524, 0, 3,
                                                                       37332, 12300, 37458, 1828,
                                                                       1856, 13896, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 40692, 0, 3,
                                                                       37458, 12363, 37584, 1856,
                                                                       1884, 13980, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 40860, 0, 3,
                                                                       37584, 12426, 37710, 1884,
                                                                       1912, 14064, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 41028, 0, 3,
                                                                       37710, 12489, 37836, 1912,
                                                                       1940, 14148, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 41196, 0, 3,
                                                                       37836, 12552, 37962, 1940,
                                                                       1968, 14232, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 41364, 0, 3,
                                                                       37962, 12615, 38088, 1968,
                                                                       1996, 14316, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 41532, 0, 3,
                                                                       38088, 12678, 38214, 1996,
                                                                       2024, 14400, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 41700, 0, 3,
                                                                       38340, 12804, 38466, 2080,
                                                                       2108, 14484, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 41868, 0, 3,
                                                                       38466, 12867, 38592, 2108,
                                                                       2136, 14568, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 42036, 0, 3,
                                                                       38592, 12930, 38718, 2136,
                                                                       2164, 14652, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 42204, 0, 3,
                                                                       38718, 12993, 38844, 2164,
                                                                       2192, 14736, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 42372, 0, 3,
                                                                       38844, 13056, 38970, 2192,
                                                                       2220, 14820, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 42540, 0, 3,
                                                                       38970, 13119, 39096, 2220,
                                                                       2248, 14904, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 42708, 0, 3,
                                                                       39096, 13182, 39222, 2248,
                                                                       2276, 14988, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 42876, 0, 3,
                                                                       39222, 13245, 39348, 2276,
                                                                       2304, 15072, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 43044, 0, 3,
                                                                       39348, 13308, 39474, 2304,
                                                                       2332, 15156, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 43212, 0, 3,
                                                                       39474, 13371, 39600, 2332,
                                                                       2360, 15240, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 43380, 0, 3,
                                                                       39600, 13434, 39726, 2360,
                                                                       2388, 15324, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 43548, 0, 3,
                                                                       39852, 13560, 40020, 2444,
                                                                       2480, 15408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 43764, 0, 3,
                                                                       40020, 13644, 40188, 2480,
                                                                       2516, 15516, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 43980, 0, 3,
                                                                       40188, 13728, 40356, 2516,
                                                                       2552, 15624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 44196, 0, 3,
                                                                       40356, 13812, 40524, 2552,
                                                                       2588, 15732, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 44412, 0, 3,
                                                                       40524, 13896, 40692, 2588,
                                                                       2624, 15840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 44628, 0, 3,
                                                                       40692, 13980, 40860, 2624,
                                                                       2660, 15948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 44844, 0, 3,
                                                                       40860, 14064, 41028, 2660,
                                                                       2696, 16056, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 45060, 0, 3,
                                                                       41028, 14148, 41196, 2696,
                                                                       2732, 16164, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 45276, 0, 3,
                                                                       41196, 14232, 41364, 2732,
                                                                       2768, 16272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 45492, 0, 3,
                                                                       41364, 14316, 41532, 2768,
                                                                       2804, 16380, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 45708, 0, 3,
                                                                       41700, 14484, 41868, 2876,
                                                                       2912, 16488, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 45924, 0, 3,
                                                                       41868, 14568, 42036, 2912,
                                                                       2948, 16596, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 46140, 0, 3,
                                                                       42036, 14652, 42204, 2948,
                                                                       2984, 16704, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 46356, 0, 3,
                                                                       42204, 14736, 42372, 2984,
                                                                       3020, 16812, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 46572, 0, 3,
                                                                       42372, 14820, 42540, 3020,
                                                                       3056, 16920, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 46788, 0, 3,
                                                                       42540, 14904, 42708, 3056,
                                                                       3092, 17028, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 47004, 0, 3,
                                                                       42708, 14988, 42876, 3092,
                                                                       3128, 17136, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 47220, 0, 3,
                                                                       42876, 15072, 43044, 3128,
                                                                       3164, 17244, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 47436, 0, 3,
                                                                       43044, 15156, 43212, 3164,
                                                                       3200, 17352, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 47652, 0, 3,
                                                                       43212, 15240, 43380, 3200,
                                                                       3236, 17460, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 47868, 0, 3,
                                                                       43548, 15408, 43764, 3308,
                                                                       3353, 17568, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 48138, 0, 3,
                                                                       43764, 15516, 43980, 3353,
                                                                       3398, 17703, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 48408, 0, 3,
                                                                       43980, 15624, 44196, 3398,
                                                                       3443, 17838, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 48678, 0, 3,
                                                                       44196, 15732, 44412, 3443,
                                                                       3488, 17973, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 48948, 0, 3,
                                                                       44412, 15840, 44628, 3488,
                                                                       3533, 18108, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 49218, 0, 3,
                                                                       44628, 15948, 44844, 3533,
                                                                       3578, 18243, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 49488, 0, 3,
                                                                       44844, 16056, 45060, 3578,
                                                                       3623, 18378, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 49758, 0, 3,
                                                                       45060, 16164, 45276, 3623,
                                                                       3668, 18513, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 50028, 0, 3,
                                                                       45276, 16272, 45492, 3668,
                                                                       3713, 18648, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 50298, 0, 3,
                                                                       45708, 16488, 45924, 3803,
                                                                       3848, 18783, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 50568, 0, 3,
                                                                       45924, 16596, 46140, 3848,
                                                                       3893, 18918, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 50838, 0, 3,
                                                                       46140, 16704, 46356, 3893,
                                                                       3938, 19053, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 51108, 0, 3,
                                                                       46356, 16812, 46572, 3938,
                                                                       3983, 19188, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 51378, 0, 3,
                                                                       46572, 16920, 46788, 3983,
                                                                       4028, 19323, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 51648, 0, 3,
                                                                       46788, 17028, 47004, 4028,
                                                                       4073, 19458, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 51918, 0, 3,
                                                                       47004, 17136, 47220, 4073,
                                                                       4118, 19593, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 52188, 0, 3,
                                                                       47220, 17244, 47436, 4118,
                                                                       4163, 19728, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 52458, 0, 3,
                                                                       47436, 17352, 47652, 4163,
                                                                       4208, 19863, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 52728, 0, 3,
                                                                       47868, 17568, 48138, 4298,
                                                                       4353, 19998, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 53058, 0, 3,
                                                                       48138, 17703, 48408, 4353,
                                                                       4408, 20163, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 53388, 0, 3,
                                                                       48408, 17838, 48678, 4408,
                                                                       4463, 20328, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 53718, 0, 3,
                                                                       48678, 17973, 48948, 4463,
                                                                       4518, 20493, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 54048, 0, 3,
                                                                       48948, 18108, 49218, 4518,
                                                                       4573, 20658, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 54378, 0, 3,
                                                                       49218, 18243, 49488, 4573,
                                                                       4628, 20823, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 54708, 0, 3,
                                                                       49488, 18378, 49758, 4628,
                                                                       4683, 20988, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 55038, 0, 3,
                                                                       49758, 18513, 50028, 4683,
                                                                       4738, 21153, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 55368, 0, 3,
                                                                       50298, 18783, 50568, 4848,
                                                                       4903, 21318, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 55698, 0, 3,
                                                                       50568, 18918, 50838, 4903,
                                                                       4958, 21483, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 56028, 0, 3,
                                                                       50838, 19053, 51108, 4958,
                                                                       5013, 21648, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 56358, 0, 3,
                                                                       51108, 19188, 51378, 5013,
                                                                       5068, 21813, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 56688, 0, 3,
                                                                       51378, 19323, 51648, 5068,
                                                                       5123, 21978, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 57018, 0, 3,
                                                                       51648, 19458, 51918, 5123,
                                                                       5178, 22143, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 57348, 0, 3,
                                                                       51918, 19593, 52188, 5178,
                                                                       5233, 22308, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 57678, 0, 3,
                                                                       52188, 19728, 52458, 5233,
                                                                       5288, 22473, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 58008, 0, 3,
                                                                       52728, 19998, 53058, 5398,
                                                                       5464, 22638, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 58404, 0, 3,
                                                                       53058, 20163, 53388, 5464,
                                                                       5530, 22836, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 58800, 0, 3,
                                                                       53388, 20328, 53718, 5530,
                                                                       5596, 23034, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 59196, 0, 3,
                                                                       53718, 20493, 54048, 5596,
                                                                       5662, 23232, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 59592, 0, 3,
                                                                       54048, 20658, 54378, 5662,
                                                                       5728, 23430, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 59988, 0, 3,
                                                                       54378, 20823, 54708, 5728,
                                                                       5794, 23628, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 60384, 0, 3,
                                                                       54708, 20988, 55038, 5794,
                                                                       5860, 23826, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 60780, 0, 3,
                                                                       55368, 21318, 55698, 5992,
                                                                       6058, 24024, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 61176, 0, 3,
                                                                       55698, 21483, 56028, 6058,
                                                                       6124, 24222, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 61572, 0, 3,
                                                                       56028, 21648, 56358, 6124,
                                                                       6190, 24420, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 61968, 0, 3,
                                                                       56358, 21813, 56688, 6190,
                                                                       6256, 24618, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 62364, 0, 3,
                                                                       56688, 21978, 57018, 6256,
                                                                       6322, 24816, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 62760, 0, 3,
                                                                       57018, 22143, 57348, 6322,
                                                                       6388, 25014, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 63156, 0, 3,
                                                                       57348, 22308, 57678, 6388,
                                                                       6454, 25212, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 63552, 0, 3,
                                                                       58008, 22638, 58404, 6586,
                                                                       6664, 25410, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 64020, 0, 3,
                                                                       58404, 22836, 58800, 6664,
                                                                       6742, 25644, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 64488, 0, 3,
                                                                       58800, 23034, 59196, 6742,
                                                                       6820, 25878, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 64956, 0, 3,
                                                                       59196, 23232, 59592, 6820,
                                                                       6898, 26112, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 65424, 0, 3,
                                                                       59592, 23430, 59988, 6898,
                                                                       6976, 26346, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 65892, 0, 3,
                                                                       59988, 23628, 60384, 6976,
                                                                       7054, 26580, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 66360, 0, 3,
                                                                       60780, 24024, 61176, 7210,
                                                                       7288, 26814, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 66828, 0, 3,
                                                                       61176, 24222, 61572, 7288,
                                                                       7366, 27048, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 67296, 0, 3,
                                                                       61572, 24420, 61968, 7366,
                                                                       7444, 27282, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 67764, 0, 3,
                                                                       61968, 24618, 62364, 7444,
                                                                       7522, 27516, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 68232, 0, 3,
                                                                       62364, 24816, 62760, 7522,
                                                                       7600, 27750, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 68700, 0, 3,
                                                                       62760, 25014, 63156, 7600,
                                                                       7678, 27984, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 69168, 0, 3,
                                                                       63552, 25410, 64020, 7834,
                                                                       7925, 28218, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 69714, 0, 3,
                                                                       64020, 25644, 64488, 7925,
                                                                       8016, 28491, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 70260, 0, 3,
                                                                       64488, 25878, 64956, 8016,
                                                                       8107, 28764, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 70806, 0, 3,
                                                                       64956, 26112, 65424, 8107,
                                                                       8198, 29037, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 71352, 0, 3,
                                                                       65424, 26346, 65892, 8198,
                                                                       8289, 29310, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 71898, 0, 3,
                                                                       66360, 26814, 66828, 8471,
                                                                       8562, 29583, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 72444, 0, 3,
                                                                       66828, 27048, 67296, 8562,
                                                                       8653, 29856, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 72990, 0, 3,
                                                                       67296, 27282, 67764, 8653,
                                                                       8744, 30129, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 73536, 0, 3,
                                                                       67764, 27516, 68232, 8744,
                                                                       8835, 30402, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 74082, 0, 3,
                                                                       68232, 27750, 68700, 8835,
                                                                       8926, 30675, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74628, 3, 9108,
                                                                       9111, 30960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74638, 3, 9111,
                                                                       9114, 30966, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74648, 3, 9114,
                                                                       9117, 30972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74658, 3, 9117,
                                                                       9120, 30978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74668, 3, 9120,
                                                                       9123, 30984, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74678, 3, 9123,
                                                                       9126, 30990, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74688, 3, 9126,
                                                                       9129, 30996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74698, 3, 9129,
                                                                       9132, 31002, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74708, 3, 9132,
                                                                       9135, 31008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74718, 3, 9135,
                                                                       9138, 31014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74728, 3, 9138,
                                                                       9141, 31020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74738, 3, 9141,
                                                                       9144, 31026, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74748, 3, 9144,
                                                                       9147, 31032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74758, 3, 9147,
                                                                       9150, 31038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74768, 3, 9150,
                                                                       9153, 31044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74778, 3, 9159,
                                                                       9162, 31062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74788, 3, 9162,
                                                                       9165, 31068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74798, 3, 9165,
                                                                       9168, 31074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74808, 3, 9168,
                                                                       9171, 31080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74818, 3, 9171,
                                                                       9174, 31086, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74828, 3, 9174,
                                                                       9177, 31092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74838, 3, 9177,
                                                                       9180, 31098, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74848, 3, 9180,
                                                                       9183, 31104, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74858, 3, 9183,
                                                                       9186, 31110, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74868, 3, 9186,
                                                                       9189, 31116, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74878, 3, 9189,
                                                                       9192, 31122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74888, 3, 9192,
                                                                       9195, 31128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74898, 3, 9195,
                                                                       9198, 31134, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74908, 3, 9198,
                                                                       9201, 31140, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 74918, 3, 9201,
                                                                       9204, 31146, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 74928, 0, 3,
                                                                       74628, 30960, 74638, 9210,
                                                                       9219, 31188, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 74958, 0, 3,
                                                                       74638, 30966, 74648, 9219,
                                                                       9228, 31206, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 74988, 0, 3,
                                                                       74648, 30972, 74658, 9228,
                                                                       9237, 31224, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75018, 0, 3,
                                                                       74658, 30978, 74668, 9237,
                                                                       9246, 31242, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75048, 0, 3,
                                                                       74668, 30984, 74678, 9246,
                                                                       9255, 31260, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75078, 0, 3,
                                                                       74678, 30990, 74688, 9255,
                                                                       9264, 31278, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75108, 0, 3,
                                                                       74688, 30996, 74698, 9264,
                                                                       9273, 31296, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75138, 0, 3,
                                                                       74698, 31002, 74708, 9273,
                                                                       9282, 31314, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75168, 0, 3,
                                                                       74708, 31008, 74718, 9282,
                                                                       9291, 31332, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75198, 0, 3,
                                                                       74718, 31014, 74728, 9291,
                                                                       9300, 31350, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75228, 0, 3,
                                                                       74728, 31020, 74738, 9300,
                                                                       9309, 31368, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75258, 0, 3,
                                                                       74738, 31026, 74748, 9309,
                                                                       9318, 31386, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75288, 0, 3,
                                                                       74748, 31032, 74758, 9318,
                                                                       9327, 31404, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75318, 0, 3,
                                                                       74758, 31038, 74768, 9327,
                                                                       9336, 31422, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75348, 0, 3,
                                                                       74778, 31062, 74788, 9354,
                                                                       9363, 31476, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75378, 0, 3,
                                                                       74788, 31068, 74798, 9363,
                                                                       9372, 31494, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75408, 0, 3,
                                                                       74798, 31074, 74808, 9372,
                                                                       9381, 31512, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75438, 0, 3,
                                                                       74808, 31080, 74818, 9381,
                                                                       9390, 31530, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75468, 0, 3,
                                                                       74818, 31086, 74828, 9390,
                                                                       9399, 31548, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75498, 0, 3,
                                                                       74828, 31092, 74838, 9399,
                                                                       9408, 31566, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75528, 0, 3,
                                                                       74838, 31098, 74848, 9408,
                                                                       9417, 31584, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75558, 0, 3,
                                                                       74848, 31104, 74858, 9417,
                                                                       9426, 31602, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75588, 0, 3,
                                                                       74858, 31110, 74868, 9426,
                                                                       9435, 31620, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75618, 0, 3,
                                                                       74868, 31116, 74878, 9435,
                                                                       9444, 31638, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75648, 0, 3,
                                                                       74878, 31122, 74888, 9444,
                                                                       9453, 31656, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75678, 0, 3,
                                                                       74888, 31128, 74898, 9453,
                                                                       9462, 31674, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75708, 0, 3,
                                                                       74898, 31134, 74908, 9462,
                                                                       9471, 31692, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 75738, 0, 3,
                                                                       74908, 31140, 74918, 9471,
                                                                       9480, 31710, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 75768, 0, 3,
                                                                       74928, 31188, 74958, 9498,
                                                                       9516, 31800, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 75828, 0, 3,
                                                                       74958, 31206, 74988, 9516,
                                                                       9534, 31836, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 75888, 0, 3,
                                                                       74988, 31224, 75018, 9534,
                                                                       9552, 31872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 75948, 0, 3,
                                                                       75018, 31242, 75048, 9552,
                                                                       9570, 31908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76008, 0, 3,
                                                                       75048, 31260, 75078, 9570,
                                                                       9588, 31944, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76068, 0, 3,
                                                                       75078, 31278, 75108, 9588,
                                                                       9606, 31980, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76128, 0, 3,
                                                                       75108, 31296, 75138, 9606,
                                                                       9624, 32016, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76188, 0, 3,
                                                                       75138, 31314, 75168, 9624,
                                                                       9642, 32052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76248, 0, 3,
                                                                       75168, 31332, 75198, 9642,
                                                                       9660, 32088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76308, 0, 3,
                                                                       75198, 31350, 75228, 9660,
                                                                       9678, 32124, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76368, 0, 3,
                                                                       75228, 31368, 75258, 9678,
                                                                       9696, 32160, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76428, 0, 3,
                                                                       75258, 31386, 75288, 9696,
                                                                       9714, 32196, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76488, 0, 3,
                                                                       75288, 31404, 75318, 9714,
                                                                       9732, 32232, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76548, 0, 3,
                                                                       75348, 31476, 75378, 9768,
                                                                       9786, 32340, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76608, 0, 3,
                                                                       75378, 31494, 75408, 9786,
                                                                       9804, 32376, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76668, 0, 3,
                                                                       75408, 31512, 75438, 9804,
                                                                       9822, 32412, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76728, 0, 3,
                                                                       75438, 31530, 75468, 9822,
                                                                       9840, 32448, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76788, 0, 3,
                                                                       75468, 31548, 75498, 9840,
                                                                       9858, 32484, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76848, 0, 3,
                                                                       75498, 31566, 75528, 9858,
                                                                       9876, 32520, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76908, 0, 3,
                                                                       75528, 31584, 75558, 9876,
                                                                       9894, 32556, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 76968, 0, 3,
                                                                       75558, 31602, 75588, 9894,
                                                                       9912, 32592, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 77028, 0, 3,
                                                                       75588, 31620, 75618, 9912,
                                                                       9930, 32628, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 77088, 0, 3,
                                                                       75618, 31638, 75648, 9930,
                                                                       9948, 32664, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 77148, 0, 3,
                                                                       75648, 31656, 75678, 9948,
                                                                       9966, 32700, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 77208, 0, 3,
                                                                       75678, 31674, 75708, 9966,
                                                                       9984, 32736, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 77268, 0, 3,
                                                                       75708, 31692, 75738, 9984,
                                                                       10002, 32772, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 77328, 0, 3,
                                                                       75768, 31800, 75828,
                                                                       10038, 10068, 32928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 77428, 0, 3,
                                                                       75828, 31836, 75888,
                                                                       10068, 10098, 32988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 77528, 0, 3,
                                                                       75888, 31872, 75948,
                                                                       10098, 10128, 33048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 77628, 0, 3,
                                                                       75948, 31908, 76008,
                                                                       10128, 10158, 33108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 77728, 0, 3,
                                                                       76008, 31944, 76068,
                                                                       10158, 10188, 33168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 77828, 0, 3,
                                                                       76068, 31980, 76128,
                                                                       10188, 10218, 33228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 77928, 0, 3,
                                                                       76128, 32016, 76188,
                                                                       10218, 10248, 33288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78028, 0, 3,
                                                                       76188, 32052, 76248,
                                                                       10248, 10278, 33348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78128, 0, 3,
                                                                       76248, 32088, 76308,
                                                                       10278, 10308, 33408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78228, 0, 3,
                                                                       76308, 32124, 76368,
                                                                       10308, 10338, 33468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78328, 0, 3,
                                                                       76368, 32160, 76428,
                                                                       10338, 10368, 33528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78428, 0, 3,
                                                                       76428, 32196, 76488,
                                                                       10368, 10398, 33588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78528, 0, 3,
                                                                       76548, 32340, 76608,
                                                                       10458, 10488, 33768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78628, 0, 3,
                                                                       76608, 32376, 76668,
                                                                       10488, 10518, 33828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78728, 0, 3,
                                                                       76668, 32412, 76728,
                                                                       10518, 10548, 33888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78828, 0, 3,
                                                                       76728, 32448, 76788,
                                                                       10548, 10578, 33948,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 78928, 0, 3,
                                                                       76788, 32484, 76848,
                                                                       10578, 10608, 34008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 79028, 0, 3,
                                                                       76848, 32520, 76908,
                                                                       10608, 10638, 34068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 79128, 0, 3,
                                                                       76908, 32556, 76968,
                                                                       10638, 10668, 34128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 79228, 0, 3,
                                                                       76968, 32592, 77028,
                                                                       10668, 10698, 34188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 79328, 0, 3,
                                                                       77028, 32628, 77088,
                                                                       10698, 10728, 34248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 79428, 0, 3,
                                                                       77088, 32664, 77148,
                                                                       10728, 10758, 34308,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 79528, 0, 3,
                                                                       77148, 32700, 77208,
                                                                       10758, 10788, 34368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 79628, 0, 3,
                                                                       77208, 32736, 77268,
                                                                       10788, 10818, 34428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 79728, 0, 3,
                                                                       77328, 32928, 77428,
                                                                       10878, 10923, 34668,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 79878, 0, 3,
                                                                       77428, 32988, 77528,
                                                                       10923, 10968, 34758,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 80028, 0, 3,
                                                                       77528, 33048, 77628,
                                                                       10968, 11013, 34848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 80178, 0, 3,
                                                                       77628, 33108, 77728,
                                                                       11013, 11058, 34938,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 80328, 0, 3,
                                                                       77728, 33168, 77828,
                                                                       11058, 11103, 35028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 80478, 0, 3,
                                                                       77828, 33228, 77928,
                                                                       11103, 11148, 35118,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 80628, 0, 3,
                                                                       77928, 33288, 78028,
                                                                       11148, 11193, 35208,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 80778, 0, 3,
                                                                       78028, 33348, 78128,
                                                                       11193, 11238, 35298,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 80928, 0, 3,
                                                                       78128, 33408, 78228,
                                                                       11238, 11283, 35388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 81078, 0, 3,
                                                                       78228, 33468, 78328,
                                                                       11283, 11328, 35478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 81228, 0, 3,
                                                                       78328, 33528, 78428,
                                                                       11328, 11373, 35568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 81378, 0, 3,
                                                                       78528, 33768, 78628,
                                                                       11463, 11508, 35838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 81528, 0, 3,
                                                                       78628, 33828, 78728,
                                                                       11508, 11553, 35928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 81678, 0, 3,
                                                                       78728, 33888, 78828,
                                                                       11553, 11598, 36018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 81828, 0, 3,
                                                                       78828, 33948, 78928,
                                                                       11598, 11643, 36108,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 81978, 0, 3,
                                                                       78928, 34008, 79028,
                                                                       11643, 11688, 36198,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 82128, 0, 3,
                                                                       79028, 34068, 79128,
                                                                       11688, 11733, 36288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 82278, 0, 3,
                                                                       79128, 34128, 79228,
                                                                       11733, 11778, 36378,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 82428, 0, 3,
                                                                       79228, 34188, 79328,
                                                                       11778, 11823, 36468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 82578, 0, 3,
                                                                       79328, 34248, 79428,
                                                                       11823, 11868, 36558,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 82728, 0, 3,
                                                                       79428, 34308, 79528,
                                                                       11868, 11913, 36648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 82878, 0, 3,
                                                                       79528, 34368, 79628,
                                                                       11913, 11958, 36738,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 83028, 0, 3,
                                                                       79728, 34668, 79878,
                                                                       12048, 12111, 37080,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 83238, 0, 3,
                                                                       79878, 34758, 80028,
                                                                       12111, 12174, 37206,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 83448, 0, 3,
                                                                       80028, 34848, 80178,
                                                                       12174, 12237, 37332,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 83658, 0, 3,
                                                                       80178, 34938, 80328,
                                                                       12237, 12300, 37458,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 83868, 0, 3,
                                                                       80328, 35028, 80478,
                                                                       12300, 12363, 37584,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 84078, 0, 3,
                                                                       80478, 35118, 80628,
                                                                       12363, 12426, 37710,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 84288, 0, 3,
                                                                       80628, 35208, 80778,
                                                                       12426, 12489, 37836,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 84498, 0, 3,
                                                                       80778, 35298, 80928,
                                                                       12489, 12552, 37962,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 84708, 0, 3,
                                                                       80928, 35388, 81078,
                                                                       12552, 12615, 38088,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 84918, 0, 3,
                                                                       81078, 35478, 81228,
                                                                       12615, 12678, 38214,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 85128, 0, 3,
                                                                       81378, 35838, 81528,
                                                                       12804, 12867, 38592,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 85338, 0, 3,
                                                                       81528, 35928, 81678,
                                                                       12867, 12930, 38718,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 85548, 0, 3,
                                                                       81678, 36018, 81828,
                                                                       12930, 12993, 38844,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 85758, 0, 3,
                                                                       81828, 36108, 81978,
                                                                       12993, 13056, 38970,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 85968, 0, 3,
                                                                       81978, 36198, 82128,
                                                                       13056, 13119, 39096,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 86178, 0, 3,
                                                                       82128, 36288, 82278,
                                                                       13119, 13182, 39222,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 86388, 0, 3,
                                                                       82278, 36378, 82428,
                                                                       13182, 13245, 39348,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 86598, 0, 3,
                                                                       82428, 36468, 82578,
                                                                       13245, 13308, 39474,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 86808, 0, 3,
                                                                       82578, 36558, 82728,
                                                                       13308, 13371, 39600,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 87018, 0, 3,
                                                                       82728, 36648, 82878,
                                                                       13371, 13434, 39726,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 87228, 0, 3,
                                                                       83028, 37080, 83238,
                                                                       13560, 13644, 40188,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 87508, 0, 3,
                                                                       83238, 37206, 83448,
                                                                       13644, 13728, 40356,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 87788, 0, 3,
                                                                       83448, 37332, 83658,
                                                                       13728, 13812, 40524,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 88068, 0, 3,
                                                                       83658, 37458, 83868,
                                                                       13812, 13896, 40692,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 88348, 0, 3,
                                                                       83868, 37584, 84078,
                                                                       13896, 13980, 40860,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 88628, 0, 3,
                                                                       84078, 37710, 84288,
                                                                       13980, 14064, 41028,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 88908, 0, 3,
                                                                       84288, 37836, 84498,
                                                                       14064, 14148, 41196,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 89188, 0, 3,
                                                                       84498, 37962, 84708,
                                                                       14148, 14232, 41364,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 89468, 0, 3,
                                                                       84708, 38088, 84918,
                                                                       14232, 14316, 41532,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 89748, 0, 3,
                                                                       85128, 38592, 85338,
                                                                       14484, 14568, 42036,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 90028, 0, 3,
                                                                       85338, 38718, 85548,
                                                                       14568, 14652, 42204,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 90308, 0, 3,
                                                                       85548, 38844, 85758,
                                                                       14652, 14736, 42372,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 90588, 0, 3,
                                                                       85758, 38970, 85968,
                                                                       14736, 14820, 42540,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 90868, 0, 3,
                                                                       85968, 39096, 86178,
                                                                       14820, 14904, 42708,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 91148, 0, 3,
                                                                       86178, 39222, 86388,
                                                                       14904, 14988, 42876,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 91428, 0, 3,
                                                                       86388, 39348, 86598,
                                                                       14988, 15072, 43044,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 91708, 0, 3,
                                                                       86598, 39474, 86808,
                                                                       15072, 15156, 43212,
                                                                       ncols, gamma, p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 91988, 0, 3,
                                                                       86808, 39600, 87018,
                                                                       15156, 15240, 43380,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 92268, 0, 3,
                                                                       87228, 40188, 87508,
                                                                       15408, 15516, 43980,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 92628, 0, 3,
                                                                       87508, 40356, 87788,
                                                                       15516, 15624, 44196,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 92988, 0, 3,
                                                                       87788, 40524, 88068,
                                                                       15624, 15732, 44412,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 93348, 0, 3,
                                                                       88068, 40692, 88348,
                                                                       15732, 15840, 44628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 93708, 0, 3,
                                                                       88348, 40860, 88628,
                                                                       15840, 15948, 44844,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 94068, 0, 3,
                                                                       88628, 41028, 88908,
                                                                       15948, 16056, 45060,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 94428, 0, 3,
                                                                       88908, 41196, 89188,
                                                                       16056, 16164, 45276,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 94788, 0, 3,
                                                                       89188, 41364, 89468,
                                                                       16164, 16272, 45492,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 95148, 0, 3,
                                                                       89748, 42036, 90028,
                                                                       16488, 16596, 46140,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 95508, 0, 3,
                                                                       90028, 42204, 90308,
                                                                       16596, 16704, 46356,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 95868, 0, 3,
                                                                       90308, 42372, 90588,
                                                                       16704, 16812, 46572,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 96228, 0, 3,
                                                                       90588, 42540, 90868,
                                                                       16812, 16920, 46788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 96588, 0, 3,
                                                                       90868, 42708, 91148,
                                                                       16920, 17028, 47004,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 96948, 0, 3,
                                                                       91148, 42876, 91428,
                                                                       17028, 17136, 47220,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 97308, 0, 3,
                                                                       91428, 43044, 91708,
                                                                       17136, 17244, 47436,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 97668, 0, 3,
                                                                       91708, 43212, 91988,
                                                                       17244, 17352, 47652,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 98028, 0, 3,
                                                                       92268, 43980, 92628,
                                                                       17568, 17703, 48408,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 98478, 0, 3,
                                                                       92628, 44196, 92988,
                                                                       17703, 17838, 48678,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 98928, 0, 3,
                                                                       92988, 44412, 93348,
                                                                       17838, 17973, 48948,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 99378, 0, 3,
                                                                       93348, 44628, 93708,
                                                                       17973, 18108, 49218,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 99828, 0, 3,
                                                                       93708, 44844, 94068,
                                                                       18108, 18243, 49488,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 100278, 0, 3,
                                                                       94068, 45060, 94428,
                                                                       18243, 18378, 49758,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 100728, 0, 3,
                                                                       94428, 45276, 94788,
                                                                       18378, 18513, 50028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 101178, 0, 3,
                                                                       95148, 46140, 95508,
                                                                       18783, 18918, 50838,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 101628, 0, 3,
                                                                       95508, 46356, 95868,
                                                                       18918, 19053, 51108,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 102078, 0, 3,
                                                                       95868, 46572, 96228,
                                                                       19053, 19188, 51378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 102528, 0, 3,
                                                                       96228, 46788, 96588,
                                                                       19188, 19323, 51648,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 102978, 0, 3,
                                                                       96588, 47004, 96948,
                                                                       19323, 19458, 51918,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 103428, 0, 3,
                                                                       96948, 47220, 97308,
                                                                       19458, 19593, 52188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 103878, 0, 3,
                                                                       97308, 47436, 97668,
                                                                       19593, 19728, 52458,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 104328, 0, 3,
                                                                       98028, 48408, 98478,
                                                                       19998, 20163, 53388,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 104878, 0, 3,
                                                                       98478, 48678, 98928,
                                                                       20163, 20328, 53718,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 105428, 0, 3,
                                                                       98928, 48948, 99378,
                                                                       20328, 20493, 54048,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 105978, 0, 3,
                                                                       99378, 49218, 99828,
                                                                       20493, 20658, 54378,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 106528, 0, 3,
                                                                       99828, 49488, 100278,
                                                                       20658, 20823, 54708,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 107078, 0, 3,
                                                                       100278, 49758, 100728,
                                                                       20823, 20988, 55038,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 107628, 0, 3,
                                                                       101178, 50838, 101628,
                                                                       21318, 21483, 56028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 108178, 0, 3,
                                                                       101628, 51108, 102078,
                                                                       21483, 21648, 56358,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 108728, 0, 3,
                                                                       102078, 51378, 102528,
                                                                       21648, 21813, 56688,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 109278, 0, 3,
                                                                       102528, 51648, 102978,
                                                                       21813, 21978, 57018,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 109828, 0, 3,
                                                                       102978, 51918, 103428,
                                                                       21978, 22143, 57348,
                                                                       ncols, gamma, p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 110378, 0, 3,
                                                                       103428, 52188, 103878,
                                                                       22143, 22308, 57678,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 110928, 0, 3,
                                                                       104328, 53388, 104878,
                                                                       22638, 22836, 58800,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 111588, 0, 3,
                                                                       104878, 53718, 105428,
                                                                       22836, 23034, 59196,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 112248, 0, 3,
                                                                       105428, 54048, 105978,
                                                                       23034, 23232, 59592,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 112908, 0, 3,
                                                                       105978, 54378, 106528,
                                                                       23232, 23430, 59988,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 113568, 0, 3,
                                                                       106528, 54708, 107078,
                                                                       23430, 23628, 60384,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 114228, 0, 3,
                                                                       107628, 56028, 108178,
                                                                       24024, 24222, 61572,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 114888, 0, 3,
                                                                       108178, 56358, 108728,
                                                                       24222, 24420, 61968,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 115548, 0, 3,
                                                                       108728, 56688, 109278,
                                                                       24420, 24618, 62364,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 116208, 0, 3,
                                                                       109278, 57018, 109828,
                                                                       24618, 24816, 62760,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 116868, 0, 3,
                                                                       109828, 57348, 110378,
                                                                       24816, 25014, 63156,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 117528, 0, 3,
                                                                       110928, 58800, 111588,
                                                                       25410, 25644, 64488,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 118308, 0, 3,
                                                                       111588, 59196, 112248,
                                                                       25644, 25878, 64956,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 119088, 0, 3,
                                                                       112248, 59592, 112908,
                                                                       25878, 26112, 65424,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 119868, 0, 3,
                                                                       112908, 59988, 113568,
                                                                       26112, 26346, 65892,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 120648, 0, 3,
                                                                       114228, 61572, 114888,
                                                                       26814, 27048, 67296,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 121428, 0, 3,
                                                                       114888, 61968, 115548,
                                                                       27048, 27282, 67764,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 122208, 0, 3,
                                                                       115548, 62364, 116208,
                                                                       27282, 27516, 68232,
                                                                       ncols, gamma, p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 122988, 0, 3,
                                                                       116208, 62760, 116868,
                                                                       27516, 27750, 68700,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 123768, 0, 3,
                                                                       117528, 64488, 118308,
                                                                       28218, 28491, 70260,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 124678, 0, 3,
                                                                       118308, 64956, 119088,
                                                                       28491, 28764, 70806,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 125588, 0, 3,
                                                                       119088, 65424, 119868,
                                                                       28764, 29037, 71352,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 126498, 0, 3,
                                                                       120648, 67296, 121428,
                                                                       29583, 29856, 72990,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 127408, 0, 3,
                                                                       121428, 67764, 122208,
                                                                       29856, 30129, 73536,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 128318, 0, 3,
                                                                       122208, 68232, 122988,
                                                                       30129, 30402, 74082,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129228, 3, 30948,
                                                                       30954, 74628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129243, 3, 30954,
                                                                       30960, 74638, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129258, 3, 30960,
                                                                       30966, 74648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129273, 3, 30966,
                                                                       30972, 74658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129288, 3, 30972,
                                                                       30978, 74668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129303, 3, 30978,
                                                                       30984, 74678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129318, 3, 30984,
                                                                       30990, 74688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129333, 3, 30990,
                                                                       30996, 74698, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129348, 3, 30996,
                                                                       31002, 74708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129363, 3, 31002,
                                                                       31008, 74718, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129378, 3, 31008,
                                                                       31014, 74728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129393, 3, 31014,
                                                                       31020, 74738, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129408, 3, 31020,
                                                                       31026, 74748, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129423, 3, 31026,
                                                                       31032, 74758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129438, 3, 31032,
                                                                       31038, 74768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129453, 3, 31050,
                                                                       31056, 74778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129468, 3, 31056,
                                                                       31062, 74788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129483, 3, 31062,
                                                                       31068, 74798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129498, 3, 31068,
                                                                       31074, 74808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129513, 3, 31074,
                                                                       31080, 74818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129528, 3, 31080,
                                                                       31086, 74828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129543, 3, 31086,
                                                                       31092, 74838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129558, 3, 31092,
                                                                       31098, 74848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129573, 3, 31098,
                                                                       31104, 74858, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129588, 3, 31104,
                                                                       31110, 74868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129603, 3, 31110,
                                                                       31116, 74878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129618, 3, 31116,
                                                                       31122, 74888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129633, 3, 31122,
                                                                       31128, 74898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129648, 3, 31128,
                                                                       31134, 74908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 129663, 3, 31134,
                                                                       31140, 74918, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 129678, 0, 3,
                                                                       129228, 74628, 129243,
                                                                       31152, 31170, 74928,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 129723, 0, 3,
                                                                       129243, 74638, 129258,
                                                                       31170, 31188, 74958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 129768, 0, 3,
                                                                       129258, 74648, 129273,
                                                                       31188, 31206, 74988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 129813, 0, 3,
                                                                       129273, 74658, 129288,
                                                                       31206, 31224, 75018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 129858, 0, 3,
                                                                       129288, 74668, 129303,
                                                                       31224, 31242, 75048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 129903, 0, 3,
                                                                       129303, 74678, 129318,
                                                                       31242, 31260, 75078,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 129948, 0, 3,
                                                                       129318, 74688, 129333,
                                                                       31260, 31278, 75108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 129993, 0, 3,
                                                                       129333, 74698, 129348,
                                                                       31278, 31296, 75138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130038, 0, 3,
                                                                       129348, 74708, 129363,
                                                                       31296, 31314, 75168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130083, 0, 3,
                                                                       129363, 74718, 129378,
                                                                       31314, 31332, 75198,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130128, 0, 3,
                                                                       129378, 74728, 129393,
                                                                       31332, 31350, 75228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130173, 0, 3,
                                                                       129393, 74738, 129408,
                                                                       31350, 31368, 75258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130218, 0, 3,
                                                                       129408, 74748, 129423,
                                                                       31368, 31386, 75288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130263, 0, 3,
                                                                       129423, 74758, 129438,
                                                                       31386, 31404, 75318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130308, 0, 3,
                                                                       129453, 74778, 129468,
                                                                       31440, 31458, 75348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130353, 0, 3,
                                                                       129468, 74788, 129483,
                                                                       31458, 31476, 75378,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130398, 0, 3,
                                                                       129483, 74798, 129498,
                                                                       31476, 31494, 75408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130443, 0, 3,
                                                                       129498, 74808, 129513,
                                                                       31494, 31512, 75438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130488, 0, 3,
                                                                       129513, 74818, 129528,
                                                                       31512, 31530, 75468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130533, 0, 3,
                                                                       129528, 74828, 129543,
                                                                       31530, 31548, 75498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130578, 0, 3,
                                                                       129543, 74838, 129558,
                                                                       31548, 31566, 75528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130623, 0, 3,
                                                                       129558, 74848, 129573,
                                                                       31566, 31584, 75558,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130668, 0, 3,
                                                                       129573, 74858, 129588,
                                                                       31584, 31602, 75588,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130713, 0, 3,
                                                                       129588, 74868, 129603,
                                                                       31602, 31620, 75618,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130758, 0, 3,
                                                                       129603, 74878, 129618,
                                                                       31620, 31638, 75648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130803, 0, 3,
                                                                       129618, 74888, 129633,
                                                                       31638, 31656, 75678,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130848, 0, 3,
                                                                       129633, 74898, 129648,
                                                                       31656, 31674, 75708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 130893, 0, 3,
                                                                       129648, 74908, 129663,
                                                                       31674, 31692, 75738,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 130938, 0, 3,
                                                                       129678, 74928, 129723,
                                                                       31728, 31764, 75768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131028, 0, 3,
                                                                       129723, 74958, 129768,
                                                                       31764, 31800, 75828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131118, 0, 3,
                                                                       129768, 74988, 129813,
                                                                       31800, 31836, 75888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131208, 0, 3,
                                                                       129813, 75018, 129858,
                                                                       31836, 31872, 75948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131298, 0, 3,
                                                                       129858, 75048, 129903,
                                                                       31872, 31908, 76008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131388, 0, 3,
                                                                       129903, 75078, 129948,
                                                                       31908, 31944, 76068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131478, 0, 3,
                                                                       129948, 75108, 129993,
                                                                       31944, 31980, 76128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131568, 0, 3,
                                                                       129993, 75138, 130038,
                                                                       31980, 32016, 76188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131658, 0, 3,
                                                                       130038, 75168, 130083,
                                                                       32016, 32052, 76248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131748, 0, 3,
                                                                       130083, 75198, 130128,
                                                                       32052, 32088, 76308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131838, 0, 3,
                                                                       130128, 75228, 130173,
                                                                       32088, 32124, 76368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 131928, 0, 3,
                                                                       130173, 75258, 130218,
                                                                       32124, 32160, 76428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132018, 0, 3,
                                                                       130218, 75288, 130263,
                                                                       32160, 32196, 76488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132108, 0, 3,
                                                                       130308, 75348, 130353,
                                                                       32268, 32304, 76548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132198, 0, 3,
                                                                       130353, 75378, 130398,
                                                                       32304, 32340, 76608,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132288, 0, 3,
                                                                       130398, 75408, 130443,
                                                                       32340, 32376, 76668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132378, 0, 3,
                                                                       130443, 75438, 130488,
                                                                       32376, 32412, 76728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132468, 0, 3,
                                                                       130488, 75468, 130533,
                                                                       32412, 32448, 76788,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132558, 0, 3,
                                                                       130533, 75498, 130578,
                                                                       32448, 32484, 76848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132648, 0, 3,
                                                                       130578, 75528, 130623,
                                                                       32484, 32520, 76908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132738, 0, 3,
                                                                       130623, 75558, 130668,
                                                                       32520, 32556, 76968,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132828, 0, 3,
                                                                       130668, 75588, 130713,
                                                                       32556, 32592, 77028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 132918, 0, 3,
                                                                       130713, 75618, 130758,
                                                                       32592, 32628, 77088,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 133008, 0, 3,
                                                                       130758, 75648, 130803,
                                                                       32628, 32664, 77148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 133098, 0, 3,
                                                                       130803, 75678, 130848,
                                                                       32664, 32700, 77208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 133188, 0, 3,
                                                                       130848, 75708, 130893,
                                                                       32700, 32736, 77268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 133278, 0, 3,
                                                                       130938, 75768, 131028,
                                                                       32808, 32868, 77328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 133428, 0, 3,
                                                                       131028, 75828, 131118,
                                                                       32868, 32928, 77428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 133578, 0, 3,
                                                                       131118, 75888, 131208,
                                                                       32928, 32988, 77528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 133728, 0, 3,
                                                                       131208, 75948, 131298,
                                                                       32988, 33048, 77628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 133878, 0, 3,
                                                                       131298, 76008, 131388,
                                                                       33048, 33108, 77728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 134028, 0, 3,
                                                                       131388, 76068, 131478,
                                                                       33108, 33168, 77828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 134178, 0, 3,
                                                                       131478, 76128, 131568,
                                                                       33168, 33228, 77928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 134328, 0, 3,
                                                                       131568, 76188, 131658,
                                                                       33228, 33288, 78028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 134478, 0, 3,
                                                                       131658, 76248, 131748,
                                                                       33288, 33348, 78128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 134628, 0, 3,
                                                                       131748, 76308, 131838,
                                                                       33348, 33408, 78228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 134778, 0, 3,
                                                                       131838, 76368, 131928,
                                                                       33408, 33468, 78328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 134928, 0, 3,
                                                                       131928, 76428, 132018,
                                                                       33468, 33528, 78428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 135078, 0, 3,
                                                                       132108, 76548, 132198,
                                                                       33648, 33708, 78528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 135228, 0, 3,
                                                                       132198, 76608, 132288,
                                                                       33708, 33768, 78628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 135378, 0, 3,
                                                                       132288, 76668, 132378,
                                                                       33768, 33828, 78728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 135528, 0, 3,
                                                                       132378, 76728, 132468,
                                                                       33828, 33888, 78828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 135678, 0, 3,
                                                                       132468, 76788, 132558,
                                                                       33888, 33948, 78928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 135828, 0, 3,
                                                                       132558, 76848, 132648,
                                                                       33948, 34008, 79028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 135978, 0, 3,
                                                                       132648, 76908, 132738,
                                                                       34008, 34068, 79128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 136128, 0, 3,
                                                                       132738, 76968, 132828,
                                                                       34068, 34128, 79228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 136278, 0, 3,
                                                                       132828, 77028, 132918,
                                                                       34128, 34188, 79328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 136428, 0, 3,
                                                                       132918, 77088, 133008,
                                                                       34188, 34248, 79428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 136578, 0, 3,
                                                                       133008, 77148, 133098,
                                                                       34248, 34308, 79528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 136728, 0, 3,
                                                                       133098, 77208, 133188,
                                                                       34308, 34368, 79628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 136878, 0, 3,
                                                                       133278, 77328, 133428,
                                                                       34488, 34578, 79728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 137103, 0, 3,
                                                                       133428, 77428, 133578,
                                                                       34578, 34668, 79878,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 137328, 0, 3,
                                                                       133578, 77528, 133728,
                                                                       34668, 34758, 80028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 137553, 0, 3,
                                                                       133728, 77628, 133878,
                                                                       34758, 34848, 80178,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 137778, 0, 3,
                                                                       133878, 77728, 134028,
                                                                       34848, 34938, 80328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 138003, 0, 3,
                                                                       134028, 77828, 134178,
                                                                       34938, 35028, 80478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 138228, 0, 3,
                                                                       134178, 77928, 134328,
                                                                       35028, 35118, 80628,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 138453, 0, 3,
                                                                       134328, 78028, 134478,
                                                                       35118, 35208, 80778,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 138678, 0, 3,
                                                                       134478, 78128, 134628,
                                                                       35208, 35298, 80928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 138903, 0, 3,
                                                                       134628, 78228, 134778,
                                                                       35298, 35388, 81078,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 139128, 0, 3,
                                                                       134778, 78328, 134928,
                                                                       35388, 35478, 81228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 139353, 0, 3,
                                                                       135078, 78528, 135228,
                                                                       35658, 35748, 81378,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 139578, 0, 3,
                                                                       135228, 78628, 135378,
                                                                       35748, 35838, 81528,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 139803, 0, 3,
                                                                       135378, 78728, 135528,
                                                                       35838, 35928, 81678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 140028, 0, 3,
                                                                       135528, 78828, 135678,
                                                                       35928, 36018, 81828,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 140253, 0, 3,
                                                                       135678, 78928, 135828,
                                                                       36018, 36108, 81978,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 140478, 0, 3,
                                                                       135828, 79028, 135978,
                                                                       36108, 36198, 82128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 140703, 0, 3,
                                                                       135978, 79128, 136128,
                                                                       36198, 36288, 82278,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 140928, 0, 3,
                                                                       136128, 79228, 136278,
                                                                       36288, 36378, 82428,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 141153, 0, 3,
                                                                       136278, 79328, 136428,
                                                                       36378, 36468, 82578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 141378, 0, 3,
                                                                       136428, 79428, 136578,
                                                                       36468, 36558, 82728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 141603, 0, 3,
                                                                       136578, 79528, 136728,
                                                                       36558, 36648, 82878,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 141828, 0, 3,
                                                                       136878, 79728, 137103,
                                                                       36828, 36954, 83028,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 142143, 0, 3,
                                                                       137103, 79878, 137328,
                                                                       36954, 37080, 83238,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 142458, 0, 3,
                                                                       137328, 80028, 137553,
                                                                       37080, 37206, 83448,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 142773, 0, 3,
                                                                       137553, 80178, 137778,
                                                                       37206, 37332, 83658,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 143088, 0, 3,
                                                                       137778, 80328, 138003,
                                                                       37332, 37458, 83868,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 143403, 0, 3,
                                                                       138003, 80478, 138228,
                                                                       37458, 37584, 84078,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 143718, 0, 3,
                                                                       138228, 80628, 138453,
                                                                       37584, 37710, 84288,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 144033, 0, 3,
                                                                       138453, 80778, 138678,
                                                                       37710, 37836, 84498,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 144348, 0, 3,
                                                                       138678, 80928, 138903,
                                                                       37836, 37962, 84708,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 144663, 0, 3,
                                                                       138903, 81078, 139128,
                                                                       37962, 38088, 84918,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 144978, 0, 3,
                                                                       139353, 81378, 139578,
                                                                       38340, 38466, 85128,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 145293, 0, 3,
                                                                       139578, 81528, 139803,
                                                                       38466, 38592, 85338,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 145608, 0, 3,
                                                                       139803, 81678, 140028,
                                                                       38592, 38718, 85548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 145923, 0, 3,
                                                                       140028, 81828, 140253,
                                                                       38718, 38844, 85758,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 146238, 0, 3,
                                                                       140253, 81978, 140478,
                                                                       38844, 38970, 85968,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 146553, 0, 3,
                                                                       140478, 82128, 140703,
                                                                       38970, 39096, 86178,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 146868, 0, 3,
                                                                       140703, 82278, 140928,
                                                                       39096, 39222, 86388,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 147183, 0, 3,
                                                                       140928, 82428, 141153,
                                                                       39222, 39348, 86598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 147498, 0, 3,
                                                                       141153, 82578, 141378,
                                                                       39348, 39474, 86808,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 147813, 0, 3,
                                                                       141378, 82728, 141603,
                                                                       39474, 39600, 87018,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 148128, 0, 3,
                                                                       141828, 83028, 142143,
                                                                       39852, 40020, 87228,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 148548, 0, 3,
                                                                       142143, 83238, 142458,
                                                                       40020, 40188, 87508,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 148968, 0, 3,
                                                                       142458, 83448, 142773,
                                                                       40188, 40356, 87788,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 149388, 0, 3,
                                                                       142773, 83658, 143088,
                                                                       40356, 40524, 88068,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 149808, 0, 3,
                                                                       143088, 83868, 143403,
                                                                       40524, 40692, 88348,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 150228, 0, 3,
                                                                       143403, 84078, 143718,
                                                                       40692, 40860, 88628,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 150648, 0, 3,
                                                                       143718, 84288, 144033,
                                                                       40860, 41028, 88908,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 151068, 0, 3,
                                                                       144033, 84498, 144348,
                                                                       41028, 41196, 89188,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 151488, 0, 3,
                                                                       144348, 84708, 144663,
                                                                       41196, 41364, 89468,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 151908, 0, 3,
                                                                       144978, 85128, 145293,
                                                                       41700, 41868, 89748,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 152328, 0, 3,
                                                                       145293, 85338, 145608,
                                                                       41868, 42036, 90028,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 152748, 0, 3,
                                                                       145608, 85548, 145923,
                                                                       42036, 42204, 90308,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 153168, 0, 3,
                                                                       145923, 85758, 146238,
                                                                       42204, 42372, 90588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 153588, 0, 3,
                                                                       146238, 85968, 146553,
                                                                       42372, 42540, 90868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 154008, 0, 3,
                                                                       146553, 86178, 146868,
                                                                       42540, 42708, 91148,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 154428, 0, 3,
                                                                       146868, 86388, 147183,
                                                                       42708, 42876, 91428,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 154848, 0, 3,
                                                                       147183, 86598, 147498,
                                                                       42876, 43044, 91708,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 155268, 0, 3,
                                                                       147498, 86808, 147813,
                                                                       43044, 43212, 91988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 155688, 0, 3,
                                                                       148128, 87228, 148548,
                                                                       43548, 43764, 92268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 156228, 0, 3,
                                                                       148548, 87508, 148968,
                                                                       43764, 43980, 92628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 156768, 0, 3,
                                                                       148968, 87788, 149388,
                                                                       43980, 44196, 92988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 157308, 0, 3,
                                                                       149388, 88068, 149808,
                                                                       44196, 44412, 93348,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 157848, 0, 3,
                                                                       149808, 88348, 150228,
                                                                       44412, 44628, 93708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 158388, 0, 3,
                                                                       150228, 88628, 150648,
                                                                       44628, 44844, 94068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 158928, 0, 3,
                                                                       150648, 88908, 151068,
                                                                       44844, 45060, 94428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 159468, 0, 3,
                                                                       151068, 89188, 151488,
                                                                       45060, 45276, 94788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 160008, 0, 3,
                                                                       151908, 89748, 152328,
                                                                       45708, 45924, 95148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 160548, 0, 3,
                                                                       152328, 90028, 152748,
                                                                       45924, 46140, 95508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 161088, 0, 3,
                                                                       152748, 90308, 153168,
                                                                       46140, 46356, 95868,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 161628, 0, 3,
                                                                       153168, 90588, 153588,
                                                                       46356, 46572, 96228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 162168, 0, 3,
                                                                       153588, 90868, 154008,
                                                                       46572, 46788, 96588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 162708, 0, 3,
                                                                       154008, 91148, 154428,
                                                                       46788, 47004, 96948,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 163248, 0, 3,
                                                                       154428, 91428, 154848,
                                                                       47004, 47220, 97308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 163788, 0, 3,
                                                                       154848, 91708, 155268,
                                                                       47220, 47436, 97668,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 164328, 0, 3,
                                                                       155688, 92268, 156228,
                                                                       47868, 48138, 98028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 165003, 0, 3,
                                                                       156228, 92628, 156768,
                                                                       48138, 48408, 98478,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 165678, 0, 3,
                                                                       156768, 92988, 157308,
                                                                       48408, 48678, 98928,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 166353, 0, 3,
                                                                       157308, 93348, 157848,
                                                                       48678, 48948, 99378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 167028, 0, 3,
                                                                       157848, 93708, 158388,
                                                                       48948, 49218, 99828,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 167703, 0, 3,
                                                                       158388, 94068, 158928,
                                                                       49218, 49488, 100278,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 168378, 0, 3,
                                                                       158928, 94428, 159468,
                                                                       49488, 49758, 100728,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 169053, 0, 3,
                                                                       160008, 95148, 160548,
                                                                       50298, 50568, 101178,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 169728, 0, 3,
                                                                       160548, 95508, 161088,
                                                                       50568, 50838, 101628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 170403, 0, 3,
                                                                       161088, 95868, 161628,
                                                                       50838, 51108, 102078,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 171078, 0, 3,
                                                                       161628, 96228, 162168,
                                                                       51108, 51378, 102528,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 171753, 0, 3,
                                                                       162168, 96588, 162708,
                                                                       51378, 51648, 102978,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 172428, 0, 3,
                                                                       162708, 96948, 163248,
                                                                       51648, 51918, 103428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 173103, 0, 3,
                                                                       163248, 97308, 163788,
                                                                       51918, 52188, 103878,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 173778, 0, 3,
                                                                       164328, 98028, 165003,
                                                                       52728, 53058, 104328,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 174603, 0, 3,
                                                                       165003, 98478, 165678,
                                                                       53058, 53388, 104878,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 175428, 0, 3,
                                                                       165678, 98928, 166353,
                                                                       53388, 53718, 105428,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 176253, 0, 3,
                                                                       166353, 99378, 167028,
                                                                       53718, 54048, 105978,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 177078, 0, 3,
                                                                       167028, 99828, 167703,
                                                                       54048, 54378, 106528,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 177903, 0, 3,
                                                                       167703, 100278, 168378,
                                                                       54378, 54708, 107078,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 178728, 0, 3,
                                                                       169053, 101178, 169728,
                                                                       55368, 55698, 107628,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 179553, 0, 3,
                                                                       169728, 101628, 170403,
                                                                       55698, 56028, 108178,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 180378, 0, 3,
                                                                       170403, 102078, 171078,
                                                                       56028, 56358, 108728,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 181203, 0, 3,
                                                                       171078, 102528, 171753,
                                                                       56358, 56688, 109278,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 182028, 0, 3,
                                                                       171753, 102978, 172428,
                                                                       56688, 57018, 109828,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 182853, 0, 3,
                                                                       172428, 103428, 173103,
                                                                       57018, 57348, 110378,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 183678, 0, 3,
                                                                       173778, 104328, 174603,
                                                                       58008, 58404, 110928,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 184668, 0, 3,
                                                                       174603, 104878, 175428,
                                                                       58404, 58800, 111588,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 185658, 0, 3,
                                                                       175428, 105428, 176253,
                                                                       58800, 59196, 112248,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 186648, 0, 3,
                                                                       176253, 105978, 177078,
                                                                       59196, 59592, 112908,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 187638, 0, 3,
                                                                       177078, 106528, 177903,
                                                                       59592, 59988, 113568,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 188628, 0, 3,
                                                                       178728, 107628, 179553,
                                                                       60780, 61176, 114228,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 189618, 0, 3,
                                                                       179553, 108178, 180378,
                                                                       61176, 61572, 114888,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 190608, 0, 3,
                                                                       180378, 108728, 181203,
                                                                       61572, 61968, 115548,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 191598, 0, 3,
                                                                       181203, 109278, 182028,
                                                                       61968, 62364, 116208,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 192588, 0, 3,
                                                                       182028, 109828, 182853,
                                                                       62364, 62760, 116868,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 193578, 0, 3,
                                                                       183678, 110928, 184668,
                                                                       63552, 64020, 117528,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 194748, 0, 3,
                                                                       184668, 111588, 185658,
                                                                       64020, 64488, 118308,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 195918, 0, 3,
                                                                       185658, 112248, 186648,
                                                                       64488, 64956, 119088,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 197088, 0, 3,
                                                                       186648, 112908, 187638,
                                                                       64956, 65424, 119868,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 198258, 0, 3,
                                                                       188628, 114228, 189618,
                                                                       66360, 66828, 120648,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 199428, 0, 3,
                                                                       189618, 114888, 190608,
                                                                       66828, 67296, 121428,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 200598, 0, 3,
                                                                       190608, 115548, 191598,
                                                                       67296, 67764, 122208,
                                                                       ncols, gamma, p, q);

                    compute_prim_osg_three_center_electron_repulsion_0(buffer, 201768, 0, 3,
                                                                       191598, 116208, 192588,
                                                                       67764, 68232, 122988,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 202938, 0, 3,
                                                                       193578, 117528, 194748,
                                                                       69168, 69714, 123768,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 204303, 0, 3,
                                                                       194748, 118308, 195918,
                                                                       69714, 70260, 124678,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 205668, 0, 3,
                                                                       195918, 119088, 197088,
                                                                       70260, 70806, 125588,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 207033, 0, 3,
                                                                       198258, 120648, 199428,
                                                                       71898, 72444, 126498,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 208398, 0, 3,
                                                                       199428, 121428, 200598,
                                                                       72444, 72990, 127408,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsg_three_center_electron_repulsion_0(buffer, 209763, 0, 3,
                                                                       200598, 122208, 201768,
                                                                       72990, 73536, 128318,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211128, 3, 74628,
                                                                       74638, 129258, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211149, 3, 74638,
                                                                       74648, 129273, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211170, 3, 74648,
                                                                       74658, 129288, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211191, 3, 74658,
                                                                       74668, 129303, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211212, 3, 74668,
                                                                       74678, 129318, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211233, 3, 74678,
                                                                       74688, 129333, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211254, 3, 74688,
                                                                       74698, 129348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211275, 3, 74698,
                                                                       74708, 129363, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211296, 3, 74708,
                                                                       74718, 129378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211317, 3, 74718,
                                                                       74728, 129393, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211338, 3, 74728,
                                                                       74738, 129408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211359, 3, 74738,
                                                                       74748, 129423, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211380, 3, 74748,
                                                                       74758, 129438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211401, 3, 74778,
                                                                       74788, 129483, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211422, 3, 74788,
                                                                       74798, 129498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211443, 3, 74798,
                                                                       74808, 129513, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211464, 3, 74808,
                                                                       74818, 129528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211485, 3, 74818,
                                                                       74828, 129543, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211506, 3, 74828,
                                                                       74838, 129558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211527, 3, 74838,
                                                                       74848, 129573, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211548, 3, 74848,
                                                                       74858, 129588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211569, 3, 74858,
                                                                       74868, 129603, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211590, 3, 74868,
                                                                       74878, 129618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211611, 3, 74878,
                                                                       74888, 129633, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211632, 3, 74888,
                                                                       74898, 129648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 211653, 3, 74898,
                                                                       74908, 129663, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 211674, 0, 3,
                                                                       211128, 129258, 211149,
                                                                       74928, 74958, 129768,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 211737, 0, 3,
                                                                       211149, 129273, 211170,
                                                                       74958, 74988, 129813,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 211800, 0, 3,
                                                                       211170, 129288, 211191,
                                                                       74988, 75018, 129858,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 211863, 0, 3,
                                                                       211191, 129303, 211212,
                                                                       75018, 75048, 129903,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 211926, 0, 3,
                                                                       211212, 129318, 211233,
                                                                       75048, 75078, 129948,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 211989, 0, 3,
                                                                       211233, 129333, 211254,
                                                                       75078, 75108, 129993,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212052, 0, 3,
                                                                       211254, 129348, 211275,
                                                                       75108, 75138, 130038,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212115, 0, 3,
                                                                       211275, 129363, 211296,
                                                                       75138, 75168, 130083,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212178, 0, 3,
                                                                       211296, 129378, 211317,
                                                                       75168, 75198, 130128,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212241, 0, 3,
                                                                       211317, 129393, 211338,
                                                                       75198, 75228, 130173,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212304, 0, 3,
                                                                       211338, 129408, 211359,
                                                                       75228, 75258, 130218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212367, 0, 3,
                                                                       211359, 129423, 211380,
                                                                       75258, 75288, 130263,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212430, 0, 3,
                                                                       211401, 129483, 211422,
                                                                       75348, 75378, 130398,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212493, 0, 3,
                                                                       211422, 129498, 211443,
                                                                       75378, 75408, 130443,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212556, 0, 3,
                                                                       211443, 129513, 211464,
                                                                       75408, 75438, 130488,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212619, 0, 3,
                                                                       211464, 129528, 211485,
                                                                       75438, 75468, 130533,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212682, 0, 3,
                                                                       211485, 129543, 211506,
                                                                       75468, 75498, 130578,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212745, 0, 3,
                                                                       211506, 129558, 211527,
                                                                       75498, 75528, 130623,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212808, 0, 3,
                                                                       211527, 129573, 211548,
                                                                       75528, 75558, 130668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212871, 0, 3,
                                                                       211548, 129588, 211569,
                                                                       75558, 75588, 130713,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212934, 0, 3,
                                                                       211569, 129603, 211590,
                                                                       75588, 75618, 130758,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 212997, 0, 3,
                                                                       211590, 129618, 211611,
                                                                       75618, 75648, 130803,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 213060, 0, 3,
                                                                       211611, 129633, 211632,
                                                                       75648, 75678, 130848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 213123, 0, 3,
                                                                       211632, 129648, 211653,
                                                                       75678, 75708, 130893,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 213186, 0, 3,
                                                                       211674, 129768, 211737,
                                                                       75768, 75828, 131118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 213312, 0, 3,
                                                                       211737, 129813, 211800,
                                                                       75828, 75888, 131208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 213438, 0, 3,
                                                                       211800, 129858, 211863,
                                                                       75888, 75948, 131298,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 213564, 0, 3,
                                                                       211863, 129903, 211926,
                                                                       75948, 76008, 131388,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 213690, 0, 3,
                                                                       211926, 129948, 211989,
                                                                       76008, 76068, 131478,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 213816, 0, 3,
                                                                       211989, 129993, 212052,
                                                                       76068, 76128, 131568,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 213942, 0, 3,
                                                                       212052, 130038, 212115,
                                                                       76128, 76188, 131658,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 214068, 0, 3,
                                                                       212115, 130083, 212178,
                                                                       76188, 76248, 131748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 214194, 0, 3,
                                                                       212178, 130128, 212241,
                                                                       76248, 76308, 131838,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 214320, 0, 3,
                                                                       212241, 130173, 212304,
                                                                       76308, 76368, 131928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 214446, 0, 3,
                                                                       212304, 130218, 212367,
                                                                       76368, 76428, 132018,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 214572, 0, 3,
                                                                       212430, 130398, 212493,
                                                                       76548, 76608, 132288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 214698, 0, 3,
                                                                       212493, 130443, 212556,
                                                                       76608, 76668, 132378,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 214824, 0, 3,
                                                                       212556, 130488, 212619,
                                                                       76668, 76728, 132468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 214950, 0, 3,
                                                                       212619, 130533, 212682,
                                                                       76728, 76788, 132558,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 215076, 0, 3,
                                                                       212682, 130578, 212745,
                                                                       76788, 76848, 132648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 215202, 0, 3,
                                                                       212745, 130623, 212808,
                                                                       76848, 76908, 132738,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 215328, 0, 3,
                                                                       212808, 130668, 212871,
                                                                       76908, 76968, 132828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 215454, 0, 3,
                                                                       212871, 130713, 212934,
                                                                       76968, 77028, 132918,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 215580, 0, 3,
                                                                       212934, 130758, 212997,
                                                                       77028, 77088, 133008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 215706, 0, 3,
                                                                       212997, 130803, 213060,
                                                                       77088, 77148, 133098,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 215832, 0, 3,
                                                                       213060, 130848, 213123,
                                                                       77148, 77208, 133188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 215958, 0, 3,
                                                                       213186, 131118, 213312,
                                                                       77328, 77428, 133578,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 216168, 0, 3,
                                                                       213312, 131208, 213438,
                                                                       77428, 77528, 133728,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 216378, 0, 3,
                                                                       213438, 131298, 213564,
                                                                       77528, 77628, 133878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 216588, 0, 3,
                                                                       213564, 131388, 213690,
                                                                       77628, 77728, 134028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 216798, 0, 3,
                                                                       213690, 131478, 213816,
                                                                       77728, 77828, 134178,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 217008, 0, 3,
                                                                       213816, 131568, 213942,
                                                                       77828, 77928, 134328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 217218, 0, 3,
                                                                       213942, 131658, 214068,
                                                                       77928, 78028, 134478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 217428, 0, 3,
                                                                       214068, 131748, 214194,
                                                                       78028, 78128, 134628,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 217638, 0, 3,
                                                                       214194, 131838, 214320,
                                                                       78128, 78228, 134778,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 217848, 0, 3,
                                                                       214320, 131928, 214446,
                                                                       78228, 78328, 134928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 218058, 0, 3,
                                                                       214572, 132288, 214698,
                                                                       78528, 78628, 135378,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 218268, 0, 3,
                                                                       214698, 132378, 214824,
                                                                       78628, 78728, 135528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 218478, 0, 3,
                                                                       214824, 132468, 214950,
                                                                       78728, 78828, 135678,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 218688, 0, 3,
                                                                       214950, 132558, 215076,
                                                                       78828, 78928, 135828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 218898, 0, 3,
                                                                       215076, 132648, 215202,
                                                                       78928, 79028, 135978,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 219108, 0, 3,
                                                                       215202, 132738, 215328,
                                                                       79028, 79128, 136128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 219318, 0, 3,
                                                                       215328, 132828, 215454,
                                                                       79128, 79228, 136278,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 219528, 0, 3,
                                                                       215454, 132918, 215580,
                                                                       79228, 79328, 136428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 219738, 0, 3,
                                                                       215580, 133008, 215706,
                                                                       79328, 79428, 136578,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 219948, 0, 3,
                                                                       215706, 133098, 215832,
                                                                       79428, 79528, 136728,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 220158, 0, 3,
                                                                       215958, 133578, 216168,
                                                                       79728, 79878, 137328,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 220473, 0, 3,
                                                                       216168, 133728, 216378,
                                                                       79878, 80028, 137553,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 220788, 0, 3,
                                                                       216378, 133878, 216588,
                                                                       80028, 80178, 137778,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 221103, 0, 3,
                                                                       216588, 134028, 216798,
                                                                       80178, 80328, 138003,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 221418, 0, 3,
                                                                       216798, 134178, 217008,
                                                                       80328, 80478, 138228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 221733, 0, 3,
                                                                       217008, 134328, 217218,
                                                                       80478, 80628, 138453,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 222048, 0, 3,
                                                                       217218, 134478, 217428,
                                                                       80628, 80778, 138678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 222363, 0, 3,
                                                                       217428, 134628, 217638,
                                                                       80778, 80928, 138903,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 222678, 0, 3,
                                                                       217638, 134778, 217848,
                                                                       80928, 81078, 139128,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 222993, 0, 3,
                                                                       218058, 135378, 218268,
                                                                       81378, 81528, 139803,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 223308, 0, 3,
                                                                       218268, 135528, 218478,
                                                                       81528, 81678, 140028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 223623, 0, 3,
                                                                       218478, 135678, 218688,
                                                                       81678, 81828, 140253,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 223938, 0, 3,
                                                                       218688, 135828, 218898,
                                                                       81828, 81978, 140478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 224253, 0, 3,
                                                                       218898, 135978, 219108,
                                                                       81978, 82128, 140703,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 224568, 0, 3,
                                                                       219108, 136128, 219318,
                                                                       82128, 82278, 140928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 224883, 0, 3,
                                                                       219318, 136278, 219528,
                                                                       82278, 82428, 141153,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 225198, 0, 3,
                                                                       219528, 136428, 219738,
                                                                       82428, 82578, 141378,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 225513, 0, 3,
                                                                       219738, 136578, 219948,
                                                                       82578, 82728, 141603,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 225828, 0, 3,
                                                                       220158, 137328, 220473,
                                                                       83028, 83238, 142458,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 226269, 0, 3,
                                                                       220473, 137553, 220788,
                                                                       83238, 83448, 142773,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 226710, 0, 3,
                                                                       220788, 137778, 221103,
                                                                       83448, 83658, 143088,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 227151, 0, 3,
                                                                       221103, 138003, 221418,
                                                                       83658, 83868, 143403,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 227592, 0, 3,
                                                                       221418, 138228, 221733,
                                                                       83868, 84078, 143718,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 228033, 0, 3,
                                                                       221733, 138453, 222048,
                                                                       84078, 84288, 144033,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 228474, 0, 3,
                                                                       222048, 138678, 222363,
                                                                       84288, 84498, 144348,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 228915, 0, 3,
                                                                       222363, 138903, 222678,
                                                                       84498, 84708, 144663,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 229356, 0, 3,
                                                                       222993, 139803, 223308,
                                                                       85128, 85338, 145608,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 229797, 0, 3,
                                                                       223308, 140028, 223623,
                                                                       85338, 85548, 145923,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 230238, 0, 3,
                                                                       223623, 140253, 223938,
                                                                       85548, 85758, 146238,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 230679, 0, 3,
                                                                       223938, 140478, 224253,
                                                                       85758, 85968, 146553,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 231120, 0, 3,
                                                                       224253, 140703, 224568,
                                                                       85968, 86178, 146868,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 231561, 0, 3,
                                                                       224568, 140928, 224883,
                                                                       86178, 86388, 147183,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 232002, 0, 3,
                                                                       224883, 141153, 225198,
                                                                       86388, 86598, 147498,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 232443, 0, 3,
                                                                       225198, 141378, 225513,
                                                                       86598, 86808, 147813,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 232884, 0, 3,
                                                                       225828, 142458, 226269,
                                                                       87228, 87508, 148968,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 233472, 0, 3,
                                                                       226269, 142773, 226710,
                                                                       87508, 87788, 149388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 234060, 0, 3,
                                                                       226710, 143088, 227151,
                                                                       87788, 88068, 149808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 234648, 0, 3,
                                                                       227151, 143403, 227592,
                                                                       88068, 88348, 150228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 235236, 0, 3,
                                                                       227592, 143718, 228033,
                                                                       88348, 88628, 150648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 235824, 0, 3,
                                                                       228033, 144033, 228474,
                                                                       88628, 88908, 151068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 236412, 0, 3,
                                                                       228474, 144348, 228915,
                                                                       88908, 89188, 151488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 237000, 0, 3,
                                                                       229356, 145608, 229797,
                                                                       89748, 90028, 152748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 237588, 0, 3,
                                                                       229797, 145923, 230238,
                                                                       90028, 90308, 153168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 238176, 0, 3,
                                                                       230238, 146238, 230679,
                                                                       90308, 90588, 153588,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 238764, 0, 3,
                                                                       230679, 146553, 231120,
                                                                       90588, 90868, 154008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 239352, 0, 3,
                                                                       231120, 146868, 231561,
                                                                       90868, 91148, 154428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 239940, 0, 3,
                                                                       231561, 147183, 232002,
                                                                       91148, 91428, 154848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 240528, 0, 3,
                                                                       232002, 147498, 232443,
                                                                       91428, 91708, 155268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 241116, 0, 3,
                                                                       232884, 148968, 233472,
                                                                       92268, 92628, 156768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 241872, 0, 3,
                                                                       233472, 149388, 234060,
                                                                       92628, 92988, 157308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 242628, 0, 3,
                                                                       234060, 149808, 234648,
                                                                       92988, 93348, 157848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 243384, 0, 3,
                                                                       234648, 150228, 235236,
                                                                       93348, 93708, 158388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 244140, 0, 3,
                                                                       235236, 150648, 235824,
                                                                       93708, 94068, 158928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 244896, 0, 3,
                                                                       235824, 151068, 236412,
                                                                       94068, 94428, 159468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 245652, 0, 3,
                                                                       237000, 152748, 237588,
                                                                       95148, 95508, 161088,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 246408, 0, 3,
                                                                       237588, 153168, 238176,
                                                                       95508, 95868, 161628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 247164, 0, 3,
                                                                       238176, 153588, 238764,
                                                                       95868, 96228, 162168,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 247920, 0, 3,
                                                                       238764, 154008, 239352,
                                                                       96228, 96588, 162708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 248676, 0, 3,
                                                                       239352, 154428, 239940,
                                                                       96588, 96948, 163248,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 249432, 0, 3,
                                                                       239940, 154848, 240528,
                                                                       96948, 97308, 163788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 250188, 0, 3,
                                                                       241116, 156768, 241872,
                                                                       98028, 98478, 165678,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 251133, 0, 3,
                                                                       241872, 157308, 242628,
                                                                       98478, 98928, 166353,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 252078, 0, 3,
                                                                       242628, 157848, 243384,
                                                                       98928, 99378, 167028,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 253023, 0, 3,
                                                                       243384, 158388, 244140,
                                                                       99378, 99828, 167703,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 253968, 0, 3,
                                                                       244140, 158928, 244896,
                                                                       99828, 100278, 168378,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 254913, 0, 3,
                                                                       245652, 161088, 246408,
                                                                       101178, 101628, 170403,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 255858, 0, 3,
                                                                       246408, 161628, 247164,
                                                                       101628, 102078, 171078,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 256803, 0, 3,
                                                                       247164, 162168, 247920,
                                                                       102078, 102528, 171753,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 257748, 0, 3,
                                                                       247920, 162708, 248676,
                                                                       102528, 102978, 172428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 258693, 0, 3,
                                                                       248676, 163248, 249432,
                                                                       102978, 103428, 173103,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 259638, 0, 3,
                                                                       250188, 165678, 251133,
                                                                       104328, 104878, 175428,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 260793, 0, 3,
                                                                       251133, 166353, 252078,
                                                                       104878, 105428, 176253,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 261948, 0, 3,
                                                                       252078, 167028, 253023,
                                                                       105428, 105978, 177078,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 263103, 0, 3,
                                                                       253023, 167703, 253968,
                                                                       105978, 106528, 177903,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 264258, 0, 3,
                                                                       254913, 170403, 255858,
                                                                       107628, 108178, 180378,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 265413, 0, 3,
                                                                       255858, 171078, 256803,
                                                                       108178, 108728, 181203,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 266568, 0, 3,
                                                                       256803, 171753, 257748,
                                                                       108728, 109278, 182028,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 267723, 0, 3,
                                                                       257748, 172428, 258693,
                                                                       109278, 109828, 182853,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 268878, 0, 3,
                                                                       259638, 175428, 260793,
                                                                       110928, 111588, 185658,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 270264, 0, 3,
                                                                       260793, 176253, 261948,
                                                                       111588, 112248, 186648,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 271650, 0, 3,
                                                                       261948, 177078, 263103,
                                                                       112248, 112908, 187638,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 273036, 0, 3,
                                                                       264258, 180378, 265413,
                                                                       114228, 114888, 190608,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 274422, 0, 3,
                                                                       265413, 181203, 266568,
                                                                       114888, 115548, 191598,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 275808, 0, 3,
                                                                       266568, 182028, 267723,
                                                                       115548, 116208, 192588,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 277194, 0, 3,
                                                                       268878, 185658, 270264,
                                                                       117528, 118308, 195918,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 278832, 0, 3,
                                                                       270264, 186648, 271650,
                                                                       118308, 119088, 197088,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 280470, 0, 3,
                                                                       273036, 190608, 274422,
                                                                       120648, 121428, 200598,
                                                                       ncols, gamma, p, q);

                    compute_prim_osh_three_center_electron_repulsion_0(buffer, 282108, 0, 3,
                                                                       274422, 191598, 275808,
                                                                       121428, 122208, 201768,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 283746, 0, 3,
                                                                       277194, 195918, 278832,
                                                                       123768, 124678, 205668,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsh_three_center_electron_repulsion_0(buffer, 285657, 0, 3,
                                                                       280470, 200598, 282108,
                                                                       126498, 127408, 209763,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287568, 3, 129228,
                                                                       129243, 211128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287596, 3, 129243,
                                                                       129258, 211149, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287624, 3, 129258,
                                                                       129273, 211170, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287652, 3, 129273,
                                                                       129288, 211191, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287680, 3, 129288,
                                                                       129303, 211212, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287708, 3, 129303,
                                                                       129318, 211233, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287736, 3, 129318,
                                                                       129333, 211254, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287764, 3, 129333,
                                                                       129348, 211275, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287792, 3, 129348,
                                                                       129363, 211296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287820, 3, 129363,
                                                                       129378, 211317, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287848, 3, 129378,
                                                                       129393, 211338, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287876, 3, 129393,
                                                                       129408, 211359, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287904, 3, 129408,
                                                                       129423, 211380, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287932, 3, 129453,
                                                                       129468, 211401, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287960, 3, 129468,
                                                                       129483, 211422, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 287988, 3, 129483,
                                                                       129498, 211443, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288016, 3, 129498,
                                                                       129513, 211464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288044, 3, 129513,
                                                                       129528, 211485, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288072, 3, 129528,
                                                                       129543, 211506, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288100, 3, 129543,
                                                                       129558, 211527, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288128, 3, 129558,
                                                                       129573, 211548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288156, 3, 129573,
                                                                       129588, 211569, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288184, 3, 129588,
                                                                       129603, 211590, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288212, 3, 129603,
                                                                       129618, 211611, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288240, 3, 129618,
                                                                       129633, 211632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 288268, 3, 129633,
                                                                       129648, 211653, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288296, 0, 3,
                                                                       287568, 211128, 287596,
                                                                       129678, 129723, 211674,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288380, 0, 3,
                                                                       287596, 211149, 287624,
                                                                       129723, 129768, 211737,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288464, 0, 3,
                                                                       287624, 211170, 287652,
                                                                       129768, 129813, 211800,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288548, 0, 3,
                                                                       287652, 211191, 287680,
                                                                       129813, 129858, 211863,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288632, 0, 3,
                                                                       287680, 211212, 287708,
                                                                       129858, 129903, 211926,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288716, 0, 3,
                                                                       287708, 211233, 287736,
                                                                       129903, 129948, 211989,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288800, 0, 3,
                                                                       287736, 211254, 287764,
                                                                       129948, 129993, 212052,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288884, 0, 3,
                                                                       287764, 211275, 287792,
                                                                       129993, 130038, 212115,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 288968, 0, 3,
                                                                       287792, 211296, 287820,
                                                                       130038, 130083, 212178,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289052, 0, 3,
                                                                       287820, 211317, 287848,
                                                                       130083, 130128, 212241,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289136, 0, 3,
                                                                       287848, 211338, 287876,
                                                                       130128, 130173, 212304,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289220, 0, 3,
                                                                       287876, 211359, 287904,
                                                                       130173, 130218, 212367,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289304, 0, 3,
                                                                       287932, 211401, 287960,
                                                                       130308, 130353, 212430,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289388, 0, 3,
                                                                       287960, 211422, 287988,
                                                                       130353, 130398, 212493,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289472, 0, 3,
                                                                       287988, 211443, 288016,
                                                                       130398, 130443, 212556,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289556, 0, 3,
                                                                       288016, 211464, 288044,
                                                                       130443, 130488, 212619,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289640, 0, 3,
                                                                       288044, 211485, 288072,
                                                                       130488, 130533, 212682,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289724, 0, 3,
                                                                       288072, 211506, 288100,
                                                                       130533, 130578, 212745,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289808, 0, 3,
                                                                       288100, 211527, 288128,
                                                                       130578, 130623, 212808,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289892, 0, 3,
                                                                       288128, 211548, 288156,
                                                                       130623, 130668, 212871,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 289976, 0, 3,
                                                                       288156, 211569, 288184,
                                                                       130668, 130713, 212934,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 290060, 0, 3,
                                                                       288184, 211590, 288212,
                                                                       130713, 130758, 212997,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 290144, 0, 3,
                                                                       288212, 211611, 288240,
                                                                       130758, 130803, 213060,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 290228, 0, 3,
                                                                       288240, 211632, 288268,
                                                                       130803, 130848, 213123,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 290312, 0, 3,
                                                                       288296, 211674, 288380,
                                                                       130938, 131028, 213186,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 290480, 0, 3,
                                                                       288380, 211737, 288464,
                                                                       131028, 131118, 213312,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 290648, 0, 3,
                                                                       288464, 211800, 288548,
                                                                       131118, 131208, 213438,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 290816, 0, 3,
                                                                       288548, 211863, 288632,
                                                                       131208, 131298, 213564,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 290984, 0, 3,
                                                                       288632, 211926, 288716,
                                                                       131298, 131388, 213690,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 291152, 0, 3,
                                                                       288716, 211989, 288800,
                                                                       131388, 131478, 213816,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 291320, 0, 3,
                                                                       288800, 212052, 288884,
                                                                       131478, 131568, 213942,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 291488, 0, 3,
                                                                       288884, 212115, 288968,
                                                                       131568, 131658, 214068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 291656, 0, 3,
                                                                       288968, 212178, 289052,
                                                                       131658, 131748, 214194,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 291824, 0, 3,
                                                                       289052, 212241, 289136,
                                                                       131748, 131838, 214320,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 291992, 0, 3,
                                                                       289136, 212304, 289220,
                                                                       131838, 131928, 214446,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 292160, 0, 3,
                                                                       289304, 212430, 289388,
                                                                       132108, 132198, 214572,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 292328, 0, 3,
                                                                       289388, 212493, 289472,
                                                                       132198, 132288, 214698,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 292496, 0, 3,
                                                                       289472, 212556, 289556,
                                                                       132288, 132378, 214824,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 292664, 0, 3,
                                                                       289556, 212619, 289640,
                                                                       132378, 132468, 214950,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 292832, 0, 3,
                                                                       289640, 212682, 289724,
                                                                       132468, 132558, 215076,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 293000, 0, 3,
                                                                       289724, 212745, 289808,
                                                                       132558, 132648, 215202,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 293168, 0, 3,
                                                                       289808, 212808, 289892,
                                                                       132648, 132738, 215328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 293336, 0, 3,
                                                                       289892, 212871, 289976,
                                                                       132738, 132828, 215454,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 293504, 0, 3,
                                                                       289976, 212934, 290060,
                                                                       132828, 132918, 215580,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 293672, 0, 3,
                                                                       290060, 212997, 290144,
                                                                       132918, 133008, 215706,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 293840, 0, 3,
                                                                       290144, 213060, 290228,
                                                                       133008, 133098, 215832,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 294008, 0, 3,
                                                                       290312, 213186, 290480,
                                                                       133278, 133428, 215958,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 294288, 0, 3,
                                                                       290480, 213312, 290648,
                                                                       133428, 133578, 216168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 294568, 0, 3,
                                                                       290648, 213438, 290816,
                                                                       133578, 133728, 216378,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 294848, 0, 3,
                                                                       290816, 213564, 290984,
                                                                       133728, 133878, 216588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 295128, 0, 3,
                                                                       290984, 213690, 291152,
                                                                       133878, 134028, 216798,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 295408, 0, 3,
                                                                       291152, 213816, 291320,
                                                                       134028, 134178, 217008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 295688, 0, 3,
                                                                       291320, 213942, 291488,
                                                                       134178, 134328, 217218,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 295968, 0, 3,
                                                                       291488, 214068, 291656,
                                                                       134328, 134478, 217428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 296248, 0, 3,
                                                                       291656, 214194, 291824,
                                                                       134478, 134628, 217638,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 296528, 0, 3,
                                                                       291824, 214320, 291992,
                                                                       134628, 134778, 217848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 296808, 0, 3,
                                                                       292160, 214572, 292328,
                                                                       135078, 135228, 218058,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 297088, 0, 3,
                                                                       292328, 214698, 292496,
                                                                       135228, 135378, 218268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 297368, 0, 3,
                                                                       292496, 214824, 292664,
                                                                       135378, 135528, 218478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 297648, 0, 3,
                                                                       292664, 214950, 292832,
                                                                       135528, 135678, 218688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 297928, 0, 3,
                                                                       292832, 215076, 293000,
                                                                       135678, 135828, 218898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 298208, 0, 3,
                                                                       293000, 215202, 293168,
                                                                       135828, 135978, 219108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 298488, 0, 3,
                                                                       293168, 215328, 293336,
                                                                       135978, 136128, 219318,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 298768, 0, 3,
                                                                       293336, 215454, 293504,
                                                                       136128, 136278, 219528,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 299048, 0, 3,
                                                                       293504, 215580, 293672,
                                                                       136278, 136428, 219738,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 299328, 0, 3,
                                                                       293672, 215706, 293840,
                                                                       136428, 136578, 219948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 299608, 0, 3,
                                                                       294008, 215958, 294288,
                                                                       136878, 137103, 220158,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 300028, 0, 3,
                                                                       294288, 216168, 294568,
                                                                       137103, 137328, 220473,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 300448, 0, 3,
                                                                       294568, 216378, 294848,
                                                                       137328, 137553, 220788,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 300868, 0, 3,
                                                                       294848, 216588, 295128,
                                                                       137553, 137778, 221103,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 301288, 0, 3,
                                                                       295128, 216798, 295408,
                                                                       137778, 138003, 221418,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 301708, 0, 3,
                                                                       295408, 217008, 295688,
                                                                       138003, 138228, 221733,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 302128, 0, 3,
                                                                       295688, 217218, 295968,
                                                                       138228, 138453, 222048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 302548, 0, 3,
                                                                       295968, 217428, 296248,
                                                                       138453, 138678, 222363,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 302968, 0, 3,
                                                                       296248, 217638, 296528,
                                                                       138678, 138903, 222678,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 303388, 0, 3,
                                                                       296808, 218058, 297088,
                                                                       139353, 139578, 222993,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 303808, 0, 3,
                                                                       297088, 218268, 297368,
                                                                       139578, 139803, 223308,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 304228, 0, 3,
                                                                       297368, 218478, 297648,
                                                                       139803, 140028, 223623,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 304648, 0, 3,
                                                                       297648, 218688, 297928,
                                                                       140028, 140253, 223938,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 305068, 0, 3,
                                                                       297928, 218898, 298208,
                                                                       140253, 140478, 224253,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 305488, 0, 3,
                                                                       298208, 219108, 298488,
                                                                       140478, 140703, 224568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 305908, 0, 3,
                                                                       298488, 219318, 298768,
                                                                       140703, 140928, 224883,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 306328, 0, 3,
                                                                       298768, 219528, 299048,
                                                                       140928, 141153, 225198,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 306748, 0, 3,
                                                                       299048, 219738, 299328,
                                                                       141153, 141378, 225513,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 307168, 0, 3,
                                                                       299608, 220158, 300028,
                                                                       141828, 142143, 225828,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 307756, 0, 3,
                                                                       300028, 220473, 300448,
                                                                       142143, 142458, 226269,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 308344, 0, 3,
                                                                       300448, 220788, 300868,
                                                                       142458, 142773, 226710,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 308932, 0, 3,
                                                                       300868, 221103, 301288,
                                                                       142773, 143088, 227151,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 309520, 0, 3,
                                                                       301288, 221418, 301708,
                                                                       143088, 143403, 227592,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 310108, 0, 3,
                                                                       301708, 221733, 302128,
                                                                       143403, 143718, 228033,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 310696, 0, 3,
                                                                       302128, 222048, 302548,
                                                                       143718, 144033, 228474,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 311284, 0, 3,
                                                                       302548, 222363, 302968,
                                                                       144033, 144348, 228915,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 311872, 0, 3,
                                                                       303388, 222993, 303808,
                                                                       144978, 145293, 229356,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 312460, 0, 3,
                                                                       303808, 223308, 304228,
                                                                       145293, 145608, 229797,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 313048, 0, 3,
                                                                       304228, 223623, 304648,
                                                                       145608, 145923, 230238,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 313636, 0, 3,
                                                                       304648, 223938, 305068,
                                                                       145923, 146238, 230679,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 314224, 0, 3,
                                                                       305068, 224253, 305488,
                                                                       146238, 146553, 231120,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 314812, 0, 3,
                                                                       305488, 224568, 305908,
                                                                       146553, 146868, 231561,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 315400, 0, 3,
                                                                       305908, 224883, 306328,
                                                                       146868, 147183, 232002,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 315988, 0, 3,
                                                                       306328, 225198, 306748,
                                                                       147183, 147498, 232443,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 316576, 0, 3,
                                                                       307168, 225828, 307756,
                                                                       148128, 148548, 232884,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 317360, 0, 3,
                                                                       307756, 226269, 308344,
                                                                       148548, 148968, 233472,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 318144, 0, 3,
                                                                       308344, 226710, 308932,
                                                                       148968, 149388, 234060,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 318928, 0, 3,
                                                                       308932, 227151, 309520,
                                                                       149388, 149808, 234648,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 319712, 0, 3,
                                                                       309520, 227592, 310108,
                                                                       149808, 150228, 235236,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 320496, 0, 3,
                                                                       310108, 228033, 310696,
                                                                       150228, 150648, 235824,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 321280, 0, 3,
                                                                       310696, 228474, 311284,
                                                                       150648, 151068, 236412,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 322064, 0, 3,
                                                                       311872, 229356, 312460,
                                                                       151908, 152328, 237000,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 322848, 0, 3,
                                                                       312460, 229797, 313048,
                                                                       152328, 152748, 237588,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 323632, 0, 3,
                                                                       313048, 230238, 313636,
                                                                       152748, 153168, 238176,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 324416, 0, 3,
                                                                       313636, 230679, 314224,
                                                                       153168, 153588, 238764,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 325200, 0, 3,
                                                                       314224, 231120, 314812,
                                                                       153588, 154008, 239352,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 325984, 0, 3,
                                                                       314812, 231561, 315400,
                                                                       154008, 154428, 239940,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 326768, 0, 3,
                                                                       315400, 232002, 315988,
                                                                       154428, 154848, 240528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 327552, 0, 3,
                                                                       316576, 232884, 317360,
                                                                       155688, 156228, 241116,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 328560, 0, 3,
                                                                       317360, 233472, 318144,
                                                                       156228, 156768, 241872,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 329568, 0, 3,
                                                                       318144, 234060, 318928,
                                                                       156768, 157308, 242628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 330576, 0, 3,
                                                                       318928, 234648, 319712,
                                                                       157308, 157848, 243384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 331584, 0, 3,
                                                                       319712, 235236, 320496,
                                                                       157848, 158388, 244140,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 332592, 0, 3,
                                                                       320496, 235824, 321280,
                                                                       158388, 158928, 244896,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 333600, 0, 3,
                                                                       322064, 237000, 322848,
                                                                       160008, 160548, 245652,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 334608, 0, 3,
                                                                       322848, 237588, 323632,
                                                                       160548, 161088, 246408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 335616, 0, 3,
                                                                       323632, 238176, 324416,
                                                                       161088, 161628, 247164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 336624, 0, 3,
                                                                       324416, 238764, 325200,
                                                                       161628, 162168, 247920,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 337632, 0, 3,
                                                                       325200, 239352, 325984,
                                                                       162168, 162708, 248676,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 338640, 0, 3,
                                                                       325984, 239940, 326768,
                                                                       162708, 163248, 249432,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 339648, 0, 3,
                                                                       327552, 241116, 328560,
                                                                       164328, 165003, 250188,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 340908, 0, 3,
                                                                       328560, 241872, 329568,
                                                                       165003, 165678, 251133,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 342168, 0, 3,
                                                                       329568, 242628, 330576,
                                                                       165678, 166353, 252078,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 343428, 0, 3,
                                                                       330576, 243384, 331584,
                                                                       166353, 167028, 253023,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 344688, 0, 3,
                                                                       331584, 244140, 332592,
                                                                       167028, 167703, 253968,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 345948, 0, 3,
                                                                       333600, 245652, 334608,
                                                                       169053, 169728, 254913,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 347208, 0, 3,
                                                                       334608, 246408, 335616,
                                                                       169728, 170403, 255858,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 348468, 0, 3,
                                                                       335616, 247164, 336624,
                                                                       170403, 171078, 256803,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 349728, 0, 3,
                                                                       336624, 247920, 337632,
                                                                       171078, 171753, 257748,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 350988, 0, 3,
                                                                       337632, 248676, 338640,
                                                                       171753, 172428, 258693,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 352248, 0, 3,
                                                                       339648, 250188, 340908,
                                                                       173778, 174603, 259638,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 353788, 0, 3,
                                                                       340908, 251133, 342168,
                                                                       174603, 175428, 260793,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 355328, 0, 3,
                                                                       342168, 252078, 343428,
                                                                       175428, 176253, 261948,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 356868, 0, 3,
                                                                       343428, 253023, 344688,
                                                                       176253, 177078, 263103,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 358408, 0, 3,
                                                                       345948, 254913, 347208,
                                                                       178728, 179553, 264258,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 359948, 0, 3,
                                                                       347208, 255858, 348468,
                                                                       179553, 180378, 265413,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 361488, 0, 3,
                                                                       348468, 256803, 349728,
                                                                       180378, 181203, 266568,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 363028, 0, 3,
                                                                       349728, 257748, 350988,
                                                                       181203, 182028, 267723,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 364568, 0, 3,
                                                                       352248, 259638, 353788,
                                                                       183678, 184668, 268878,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 366416, 0, 3,
                                                                       353788, 260793, 355328,
                                                                       184668, 185658, 270264,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 368264, 0, 3,
                                                                       355328, 261948, 356868,
                                                                       185658, 186648, 271650,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 370112, 0, 3,
                                                                       358408, 264258, 359948,
                                                                       188628, 189618, 273036,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 371960, 0, 3,
                                                                       359948, 265413, 361488,
                                                                       189618, 190608, 274422,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsi_three_center_electron_repulsion_0(buffer, 373808, 0, 3,
                                                                       361488, 266568, 363028,
                                                                       190608, 191598, 275808,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 375656, 0, 3,
                                                                       364568, 268878, 366416,
                                                                       193578, 194748, 277194,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 377840, 0, 3,
                                                                       366416, 270264, 368264,
                                                                       194748, 195918, 278832,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 380024, 0, 3,
                                                                       370112, 273036, 371960,
                                                                       198258, 199428, 280470,
                                                                       ncols, gamma, p, q);

                    compute_prim_osi_three_center_electron_repulsion_0(buffer, 382208, 0, 3,
                                                                       371960, 274422, 373808,
                                                                       199428, 200598, 282108,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 384392, 0, 3,
                                                                       375656, 277194, 377840,
                                                                       202938, 204303, 283746,
                                                                       ncols, gamma, p, q);

                    compute_prim_qsi_three_center_electron_repulsion_0(buffer, 386940, 0, 3,
                                                                       380024, 280470, 382208,
                                                                       207033, 208398, 285657,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 389488, 316576, 784, ncols);

                    simdfunc::contract_primitives(buffer, 390636, 322064, 784, ncols);

                    simdfunc::contract_primitives(buffer, 391784, 327552, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 393260, 333600, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 394736, 339648, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 396581, 345948, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 398426, 352248, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 400681, 358408, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 402936, 364568, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 405642, 370112, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 408348, 375656, 2184, ncols);

                    simdfunc::contract_primitives(buffer, 411546, 380024, 2184, ncols);

                    simdfunc::contract_primitives(buffer, 414744, 384392, 2548, ncols);

                    simdfunc::contract_primitives(buffer, 418475, 386940, 2548, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 390272, 389488, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 391420, 390636, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 392792, 391784, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 394268, 393260, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 395996, 394736, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 397841, 396581, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 399966, 398426, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 402221, 400681, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 404784, 402936, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 407490, 405642, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 410532, 408348, 78, 1, nmax);

        simdtrf::transform_i_inner(buffer, 413730, 411546, 78, 1, nmax);

        simdtrf::transform_i_inner(buffer, 417292, 414744, 91, 1, nmax);

        simdtrf::transform_i_inner(buffer, 421023, 418475, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 422206, 390272, 392792, 13,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 423298, 391420, 394268, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 424390, 392792, 395996, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 425794, 394268, 397841, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 427198, 395996, 399966, 13,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 428953, 397841, 402221, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 430708, 399966, 404784, 13,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 432853, 402221, 407490, 13,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 434998, 404784, 410532, 13,
                                             nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 437572, 407490, 413730, 13,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 440146, 410532, 417292, 13,
                                             nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 443188, 413730, 421023, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 446230, 422206, 424390, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 448414, 423298, 425794, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 450598, 424390, 427198, 13,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 453406, 425794, 428953, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 456214, 427198, 430708, 13,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 459724, 428953, 432853, 13,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 463234, 430708, 434998, 13,
                                             nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 467524, 432853, 437572, 13,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 471814, 434998, 440146, 13,
                                             nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 476962, 437572, 443188, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 482110, 446230, 450598, 13,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 485750, 448414, 453406, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 489390, 450598, 456214, 13,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 494070, 453406, 459724, 13,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 498750, 456214, 463234, 13,
                                             nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 504600, 459724, 467524, 13,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 510450, 463234, 471814, 13,
                                             nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 517600, 467524, 476962, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 524750, 482110, 489390, 13,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 530210, 485750, 494070, 13,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 535670, 489390, 498750, 13,
                                             nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 542690, 494070, 504600, 13,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 549710, 498750, 510450, 13,
                                             nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 558485, 504600, 517600, 13,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 567260, 524750, 535670, 13,
                                             nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 574904, 530210, 542690, 13,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 582548, 535670, 549710, 13,
                                             nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 592376, 542690, 558485, 13,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 602204, 567260, 582548, 13,
                                             nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 612396, 574904, 592376, 13,
                                             nmax);

        simdtrf::transform_i_inner(buffer, 622588, 612396, 28, 13, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 622588, 169, nmax);

        simdtrf::transform_i_inner(buffer, 622588, 602204, 28, 13, nmax);

        simdtrf::transform_i_outer(values + 2197 * nvalues + n * npairs, nvalues, buffer, 622588,
                                   169, nmax);
    }

    for (size_t m = 0; m < 4394; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
